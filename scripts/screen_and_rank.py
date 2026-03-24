#!/usr/bin/env python3
"""
screen_and_rank.py — Phase 4: Molecular docking with AutoDock Vina

Usage:
    python screen_and_rank.py [--config CONFIG] [--workers N] [--resume]

Reads selected_targets and compound library from config.yaml, runs AutoDock Vina
for each (target, compound) pair, ranks by binding affinity (ΔG), and saves:
  - results/ranking.csv
  - results/top50_poses/  (top N poses as .pdbqt)

Progress is saved incrementally; use --resume to continue an interrupted run.
"""

import argparse
import concurrent.futures
import csv
import json
import os
import re
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import yaml
from Bio.PDB import PDBParser
from openbabel import openbabel as ob
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

warnings.filterwarnings("ignore")

# ── Path helpers ──────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent


def project_path(p: str) -> Path:
    return PROJECT_DIR / p


# ── Config ────────────────────────────────────────────────────────────────────

def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# ── Binding site detection ────────────────────────────────────────────────────

def parse_site_records(pdb_path: Path):
    """
    Parse SITE records from a PDB file (fixed-width PDB format).
    Returns list of (resname, chain, resnum) tuples.
    """
    site_residues = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("SITE"):
                continue
            # Fixed-width groups start at column 18, each 11 chars wide:
            # resname(3) space chain(1) resnum(4) icode(1) space(1)
            i = 18
            while i + 10 <= len(line.rstrip()):
                group = line[i : i + 11]
                resname = group[0:3].strip()
                chain = group[4].strip() if len(group) > 4 else ""
                resnum_str = group[5:9].strip() if len(group) > 5 else ""
                if resname and chain and resnum_str:
                    try:
                        site_residues.append((resname, chain, int(resnum_str)))
                    except ValueError:
                        pass
                i += 11
    return site_residues


def get_binding_site(pdb_path: Path, config: dict, pdb_id: str = ""):
    """
    Returns (center_x, center_y, center_z, size_x, size_y, size_z).

    Priority:
      1. Per-target coordinates from config.yaml (binding_site.targets.<pdb_id>)
         — written by prepare_protein.py
      2. Manual override in config.yaml under binding_site.manual
      3. Auto-detection from PDB SITE records
      4. Fallback: geometric centre of all CA atoms + 30 Å box
    """
    bs = config.get("binding_site", {})

    # 1. Per-target coordinates (written by prepare_protein.py)
    targets_map = (bs or {}).get("targets", {})
    if pdb_id and pdb_id in targets_map:
        t = targets_map[pdb_id]
        cx = float(t["center_x"])
        cy = float(t["center_y"])
        cz = float(t["center_z"])
        sx = float(t.get("size_x", 25))
        sy = float(t.get("size_y", 25))
        sz = float(t.get("size_z", 25))
        print(f"  Binding site: from prepare_protein.py ({cx:.2f}, {cy:.2f}, {cz:.2f}), "
              f"box ({sx:.1f}x{sy:.1f}x{sz:.1f})")
        return cx, cy, cz, sx, sy, sz

    # 2. Manual override
    manual = bs.get("manual") if bs else None
    if manual:
        cx = float(manual["center_x"])
        cy = float(manual["center_y"])
        cz = float(manual["center_z"])
        sx = float(manual.get("size_x", 25))
        sy = float(manual.get("size_y", 25))
        sz = float(manual.get("size_z", 25))
        print(f"  Binding site: manual ({cx:.2f}, {cy:.2f}, {cz:.2f}), "
              f"box ({sx:.1f}x{sy:.1f}x{sz:.1f})")
        return cx, cy, cz, sx, sy, sz

    # 3. Auto-detect from SITE records
    site_residues = parse_site_records(pdb_path)
    if site_residues:
        parser = PDBParser()
        struct = parser.get_structure("X", str(pdb_path))
        model = struct[0]
        coords = []
        for resname, chain, resnum in site_residues:
            try:
                res = model[chain][resnum]
                ca = res["CA"]
                coords.append(ca.get_vector().get_array())
            except (KeyError, Exception):
                pass
        if coords:
            coords = np.array(coords)
            center = coords.mean(axis=0)
            extents = coords.max(axis=0) - coords.min(axis=0)
            # Add 5 Å padding on each side → +10 total per dimension; min 20 Å
            box = np.maximum(extents + 10, 20.0)
            cx, cy, cz = center
            sx, sy, sz = box
            print(f"  Binding site: auto-detected from SITE records "
                  f"({cx:.2f}, {cy:.2f}, {cz:.2f}), box ({sx:.1f}x{sy:.1f}x{sz:.1f})")
            return float(cx), float(cy), float(cz), float(sx), float(sy), float(sz)

    # 4. Fallback: centre of mass of all CA atoms
    print("  Warning: no SITE records found. Using geometric center of CA atoms as fallback.")
    parser = PDBParser()
    struct = parser.get_structure("X", str(pdb_path))
    model = struct[0]
    coords = []
    for chain in model:
        for res in chain:
            if "CA" in res:
                coords.append(res["CA"].get_vector().get_array())
    if not coords:
        raise ValueError(f"No CA atoms found in {pdb_path}")
    coords = np.array(coords)
    center = coords.mean(axis=0)
    cx, cy, cz = center
    sx = sy = sz = 30.0
    print(f"  Binding site: fallback to CA centroid ({cx:.2f}, {cy:.2f}, {cz:.2f}), box 30x30x30")
    return float(cx), float(cy), float(cz), sx, sy, sz


# ── Format conversions ────────────────────────────────────────────────────────

def prepare_receptor_pdbqt(pdb_path: Path, out_dir: Path) -> Path:
    """Convert protein PDB → PDBQT (rigid receptor) using OpenBabel."""
    out_path = out_dir / f"{pdb_path.stem}_receptor.pdbqt"
    if out_path.exists():
        return out_path

    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdbqt")
    conv.AddOption("r", ob.OBConversion.OUTOPTIONS)  # rigid
    mol = ob.OBMol()
    ok = conv.ReadFile(mol, str(pdb_path))
    if not ok or mol.NumAtoms() == 0:
        raise RuntimeError(f"OpenBabel could not read {pdb_path}")
    conv.WriteFile(mol, str(out_path))
    return out_path


# AutoDock Vina atom types — elements outside this set cause PDBQT parse errors
_VINA_SUPPORTED_ELEMENTS = {
    "H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I",
    "Mn", "Fe", "Zn", "Ca", "Mg", "Na", "K",
}


def mol_has_unsupported_elements(mol) -> bool:
    """Return True if the molecule contains any element unsupported by AutoDock Vina."""
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _VINA_SUPPORTED_ELEMENTS:
            return True
    return False


def prepare_ligand_pdbqt(mol, name: str, out_dir: Path) -> Path | None:
    """Convert RDKit Mol (with 3D coords) → PDBQT using Meeko.

    Skips molecules with atom types unsupported by AutoDock Vina (e.g. boron).
    """
    out_path = out_dir / f"{name}.pdbqt"
    if out_path.exists():
        return out_path
    if mol_has_unsupported_elements(mol):
        return None
    try:
        prep = MoleculePreparation()
        molsetups = prep.prepare(mol)
        if not molsetups:
            return None
        pdbqt_str, ok, msg = PDBQTWriterLegacy.write_string(molsetups[0])
        if not ok:
            return None
        out_path.write_text(pdbqt_str)
        return out_path
    except Exception as e:
        print(f"    Meeko failed for {name}: {e}", file=sys.stderr)
        return None


# ── Docking ───────────────────────────────────────────────────────────────────

def find_vina_binary() -> str:
    """Locate the AutoDock Vina binary."""
    # Check ~/bin first, then PATH
    candidates = [
        str(Path.home() / "bin" / "vina"),
        "vina",
    ]
    for candidate in candidates:
        try:
            result = subprocess.run(
                [candidate, "--version"],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                return candidate
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    raise RuntimeError(
        "AutoDock Vina not found. Install it and put it in ~/bin or PATH."
    )


def parse_vina_output(output: str) -> float | None:
    """Extract best binding affinity (mode 1 ΔG kcal/mol) from vina stdout."""
    for line in output.splitlines():
        m = re.match(r"\s*1\s+([-\d.]+)\s+", line)
        if m:
            return float(m.group(1))
    return None


def run_vina(
    vina_bin: str,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    out_pdbqt: Path,
    cx: float, cy: float, cz: float,
    sx: float, sy: float, sz: float,
    exhaustiveness: int,
    num_modes: int,
) -> float | None:
    """Run AutoDock Vina for one ligand. Returns best ΔG or None on failure."""
    cmd = [
        vina_bin,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--out", str(out_pdbqt),
        "--center_x", f"{cx:.4f}",
        "--center_y", f"{cy:.4f}",
        "--center_z", f"{cz:.4f}",
        "--size_x", f"{sx:.4f}",
        "--size_y", f"{sy:.4f}",
        "--size_z", f"{sz:.4f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300
        )
        if result.returncode != 0:
            return None
        return parse_vina_output(result.stdout)
    except subprocess.TimeoutExpired:
        return None
    except Exception:
        return None


# ── Progress tracking ─────────────────────────────────────────────────────────

def load_progress(progress_file: Path) -> dict:
    """Load previously completed docking results {compound_key: dg}."""
    if progress_file.exists():
        with open(progress_file) as f:
            return json.load(f)
    return {}


def save_progress(progress_file: Path, progress: dict):
    with open(progress_file, "w") as f:
        json.dump(progress, f, indent=2)


# ── Main pipeline ─────────────────────────────────────────────────────────────

def dock_one(args):
    """Worker function for parallel docking. Returns (key, name, pdb_id, dg)."""
    (
        vina_bin, receptor_pdbqt, lig_path, lig_name,
        cx, cy, cz, sx, sy, sz,
        exhaustiveness, num_modes,
        poses_dir, pdb_id
    ) = args

    key = f"{pdb_id}__{lig_name}"
    out_pdbqt = poses_dir / f"{key}.pdbqt"

    dg = run_vina(
        vina_bin, receptor_pdbqt, lig_path, out_pdbqt,
        cx, cy, cz, sx, sy, sz,
        exhaustiveness, num_modes,
    )
    return key, lig_name, pdb_id, dg


def main():
    parser = argparse.ArgumentParser(
        description="screen_and_rank.py — Molecular docking with AutoDock Vina"
    )
    parser.add_argument(
        "--config", default=str(project_path("config.yaml")),
        help="Path to config.yaml (default: project root)"
    )
    parser.add_argument(
        "--workers", type=int, default=None,
        help="Number of parallel docking workers (overrides config)"
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume a previously interrupted docking run"
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    config = load_config(config_path)

    # ── Resolve paths ──────────────────────────────────────────────────────────
    output_dir = project_path(config.get("output_dir", "results/"))
    output_dir.mkdir(parents=True, exist_ok=True)

    poses_dir = output_dir / "top50_poses"
    poses_dir.mkdir(exist_ok=True)

    prep_dir = project_path("data/prepared")
    prep_dir.mkdir(exist_ok=True)

    compounds_sdf = project_path("data/compounds/compounds.sdf")
    structures_dir = project_path("data/structures")

    ranking_csv = output_dir / "ranking.csv"
    progress_file = output_dir / ".docking_progress.json"

    # ── Config values ──────────────────────────────────────────────────────────
    selected_targets = config.get("selected_targets", [])
    if not selected_targets:
        print("Error: no selected_targets in config.yaml. Run select_targets.py first.")
        sys.exit(1)

    top_n = int(config.get("top_n_docking", 50))
    exhaustiveness = int(config.get("vina_exhaustiveness", 8))
    num_modes = int(config.get("vina_num_modes", 9))
    workers = args.workers or int(config.get("docking_workers", 4))

    # ── Find Vina ──────────────────────────────────────────────────────────────
    try:
        vina_bin = find_vina_binary()
        print(f"Using Vina binary: {vina_bin}")
    except RuntimeError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # ── Load compounds ─────────────────────────────────────────────────────────
    if not compounds_sdf.exists():
        print(f"Error: {compounds_sdf} not found. Run fetch_compounds.py first.")
        sys.exit(1)

    print(f"Loading compounds from {compounds_sdf}...")
    suppl = Chem.SDMolSupplier(str(compounds_sdf), removeHs=False)
    compounds = []
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"compound_{i}"
        # Sanitise name for use as filename (replace spaces/slashes)
        safe_name = re.sub(r"[^\w\-]", "_", name)[:80]
        compounds.append((safe_name, mol))

    print(f"  Loaded {len(compounds)} valid compounds.")

    if not compounds:
        print("Error: no valid compounds loaded from SDF.")
        sys.exit(1)

    # ── Load progress ──────────────────────────────────────────────────────────
    progress = {}
    if args.resume and progress_file.exists():
        progress = load_progress(progress_file)
        print(f"Resuming: {len(progress)} docking results already done.")

    # ── Prepare ligand PDBQTs ──────────────────────────────────────────────────
    print(f"Preparing ligand PDBQT files...")
    lig_pdbqt_map = {}  # safe_name → Path
    failed_prep = 0
    for safe_name, mol in compounds:
        lig_path = prepare_ligand_pdbqt(mol, safe_name, prep_dir)
        if lig_path:
            lig_pdbqt_map[safe_name] = lig_path
        else:
            failed_prep += 1

    print(f"  {len(lig_pdbqt_map)} ligands prepared, {failed_prep} failed.")

    # ── Main docking loop (per target) ────────────────────────────────────────
    all_results = []  # list of dicts

    prepared_dir = project_path("data/prepared")

    for pdb_id in selected_targets:
        # Prefer prepared structure; fall back to raw structure with warning
        prep_pdb = prepared_dir / f"{pdb_id}_prepared.pdb"
        raw_pdb  = structures_dir / f"{pdb_id}.pdb"

        if prep_pdb.exists():
            pdb_path = prep_pdb
        elif raw_pdb.exists():
            print(f"Warning: prepared structure for {pdb_id} not found. "
                  f"Using raw PDB. Run prepare_protein.py for best results.")
            pdb_path = raw_pdb
        else:
            print(f"Warning: {pdb_id}.pdb not found in data/structures/ or data/prepared/, skipping.")
            continue

        print(f"\n{'='*60}")
        print(f"Target: {pdb_id}  (source: {pdb_path.parent.name}/{pdb_path.name})")

        # Prepare receptor
        print(f"  Preparing receptor PDBQT...")
        try:
            receptor_pdbqt = prepare_receptor_pdbqt(pdb_path, prep_dir)
        except RuntimeError as e:
            print(f"  Error preparing receptor: {e}")
            continue

        # Binding site
        try:
            cx, cy, cz, sx, sy, sz = get_binding_site(pdb_path, config, pdb_id=pdb_id)
        except Exception as e:
            print(f"  Error detecting binding site: {e}")
            continue

        # Build list of jobs to run (skip already done)
        jobs = []
        for safe_name, lig_path in lig_pdbqt_map.items():
            key = f"{pdb_id}__{safe_name}"
            if key in progress:
                # Already done — add to all_results from cache
                all_results.append({
                    "pdb_id": pdb_id,
                    "compound_id": safe_name,
                    "dg_kcal_mol": progress[key],
                })
                continue
            jobs.append((
                vina_bin, receptor_pdbqt, lig_path, safe_name,
                cx, cy, cz, sx, sy, sz,
                exhaustiveness, num_modes,
                poses_dir, pdb_id,
            ))

        print(f"  Docking {len(jobs)} compounds ({len(lig_pdbqt_map) - len(jobs)} cached)...")

        done = 0
        failed = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            future_map = {executor.submit(dock_one, job): job for job in jobs}
            for future in concurrent.futures.as_completed(future_map):
                key, lig_name, pid, dg = future.result()
                done += 1
                if dg is not None:
                    all_results.append({
                        "pdb_id": pid,
                        "compound_id": lig_name,
                        "dg_kcal_mol": dg,
                    })
                    progress[key] = dg
                else:
                    failed += 1
                    progress[key] = None

                # Save progress every 10 completions
                if done % 10 == 0:
                    save_progress(progress_file, progress)
                    completed = done - failed
                    print(f"    {done}/{len(jobs)} done | {completed} successful | {failed} failed")

        save_progress(progress_file, progress)
        successful = sum(1 for r in all_results if r["pdb_id"] == pdb_id and r["dg_kcal_mol"] is not None)
        print(f"  Done: {successful} successful docking results for {pdb_id}.")

    # ── Rank and write CSV ─────────────────────────────────────────────────────
    if not all_results:
        print("\nNo docking results to rank.")
        sys.exit(0)

    # Sort by ΔG ascending (most negative = best binding)
    ranked = sorted(
        [r for r in all_results if r["dg_kcal_mol"] is not None],
        key=lambda r: r["dg_kcal_mol"]
    )

    print(f"\nWriting ranking.csv ({len(ranked)} entries)...")
    with open(ranking_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "pdb_id", "compound_id", "dg_kcal_mol"])
        writer.writeheader()
        for i, row in enumerate(ranked, 1):
            writer.writerow({
                "rank": i,
                "pdb_id": row["pdb_id"],
                "compound_id": row["compound_id"],
                "dg_kcal_mol": row["dg_kcal_mol"],
            })

    # ── Keep only top N poses on disk ─────────────────────────────────────────
    top_keys = set()
    for row in ranked[:top_n]:
        top_keys.add(f"{row['pdb_id']}__{row['compound_id']}")

    # Remove poses that are not in top N to save disk space
    for pose_file in poses_dir.glob("*.pdbqt"):
        if pose_file.stem not in top_keys:
            pose_file.unlink(missing_ok=True)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"Docking complete.")
    print(f"  Total ranked: {len(ranked)}")
    print(f"  Ranking CSV:  {ranking_csv}")
    print(f"  Top {top_n} poses: {poses_dir}")
    print(f"\nTop 10 results:")
    print(f"  {'Rank':<6} {'PDB':<6} {'Compound':<40} {'ΔG (kcal/mol)'}")
    print(f"  {'-'*6} {'-'*6} {'-'*40} {'-'*15}")
    for row in ranked[:10]:
        rank = ranked.index(row) + 1
        print(f"  {rank:<6} {row['pdb_id']:<6} {row['compound_id']:<40} {row['dg_kcal_mol']:.3f}")


if __name__ == "__main__":
    main()
