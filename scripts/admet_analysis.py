#!/usr/bin/env python3
"""
admet_analysis.py — Phase 5: ADMET analysis of top docking candidates

Usage:
    python admet_analysis.py [--config CONFIG] [--top N]

Reads results/ranking.csv, extracts SMILES for the top N compounds,
submits each to the SwissADME API, and generates:
  - results/admet_report.csv  — full ADMET analysis for all evaluated compounds
  - results/final_candidates.csv — only Green and Yellow classified compounds

Classification:
  Green  — passes all ADMET criteria (priority candidates for wet lab)
  Yellow — passes most criteria (requires manual review)
  Red    — fails key criteria (discarded)
"""

import argparse
import csv
import io
import re
import sys
import time
from pathlib import Path

import requests
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# ── Path helpers ──────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent


def project_path(p: str) -> Path:
    return PROJECT_DIR / p


# ── Config ────────────────────────────────────────────────────────────────────

def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# ── SMILES extraction from SDF ────────────────────────────────────────────────

def build_smiles_map(sdf_path: Path) -> dict[str, str]:
    """
    Returns {compound_id: canonical_smiles} from the compounds SDF.
    Falls back to generating SMILES from 3D coordinates with RDKit.
    """
    smiles_map = {}
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=True)
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else None
        if not name:
            continue
        # Prefer the stored SMILES property if available
        if mol.HasProp("smiles"):
            smi = mol.GetProp("smiles").strip()
        elif mol.HasProp("SMILES"):
            smi = mol.GetProp("SMILES").strip()
        else:
            smi = Chem.MolToSmiles(mol)
        smiles_map[name] = smi
    return smiles_map


# ── SwissADME API ─────────────────────────────────────────────────────────────

SWISSADME_URL = "https://www.swissadme.ch/index.php"
SWISSADME_CSV_TMPL = "https://www.swissadme.ch/results/{session}/swissadme.csv"


def query_swissadme(smiles: str, retries: int = 3, delay: float = 2.0) -> dict | None:
    """
    Submit a single SMILES to SwissADME and return the parsed CSV row as a dict.
    Returns None on failure after retries.
    """
    for attempt in range(retries):
        try:
            r = requests.post(
                SWISSADME_URL,
                data={"smiles": smiles},
                timeout=60,
            )
            r.raise_for_status()

            m = re.search(r"results/(\d+)/swissadme\.csv", r.text)
            if not m:
                time.sleep(delay)
                continue

            csv_url = SWISSADME_CSV_TMPL.format(session=m.group(1))
            time.sleep(1.0)  # small pause before fetching CSV
            csv_r = requests.get(csv_url, timeout=30)
            csv_r.raise_for_status()

            reader = csv.DictReader(io.StringIO(csv_r.text))
            rows = list(reader)
            if rows:
                return rows[0]

        except (requests.RequestException, Exception) as e:
            print(f"    SwissADME attempt {attempt + 1}/{retries} failed: {e}",
                  file=sys.stderr)
            time.sleep(delay * (attempt + 1))

    return None


# ── ADMET classification ──────────────────────────────────────────────────────

def classify_admet(row: dict, smiles: str) -> tuple[str, list[str], list[str]]:
    """
    Classify a SwissADME result row as Green, Yellow, or Red.
    Returns (classification, passed_criteria, failed_criteria).

    Criteria (from design spec):
      1. Lipinski Ro5: #violations == 0
      2. Oral bioavailability: Bioavailability Score >= 0.55
      3. BBB permeability: preferably No (not a strict fail)
      4. GI absorption: High preferred
      5. PAINS alerts: 0
      6. Solubility: not "Insoluble" (ESOL Class)
    """
    passed = []
    failed = []
    warnings = []

    # 1. Lipinski Ro5
    try:
        lipinski_violations = int(row.get("Lipinski #violations", "0") or "0")
    except ValueError:
        lipinski_violations = 0
    if lipinski_violations == 0:
        passed.append("Lipinski Ro5 (0 violations)")
    else:
        failed.append(f"Lipinski Ro5 ({lipinski_violations} violations)")

    # 2. Bioavailability Score
    try:
        bioavail = float(row.get("Bioavailability Score", "0") or "0")
    except ValueError:
        bioavail = 0.0
    if bioavail >= 0.55:
        passed.append(f"Bioavailability Score ({bioavail:.2f} >= 0.55)")
    else:
        failed.append(f"Bioavailability Score ({bioavail:.2f} < 0.55)")

    # 3. PAINS alerts (hard fail)
    try:
        pains = int(row.get("PAINS #alerts", "0") or "0")
    except ValueError:
        pains = 0
    if pains == 0:
        passed.append("PAINS (no alerts)")
    else:
        failed.append(f"PAINS ({pains} alert(s))")

    # 4. GI absorption (soft criterion)
    gi = row.get("GI absorption", "").strip()
    if gi == "High":
        passed.append("GI absorption (High)")
    elif gi == "Low":
        warnings.append("GI absorption (Low)")
    else:
        warnings.append(f"GI absorption ({gi})")

    # 5. Solubility (ESOL Class)
    esol_class = row.get("ESOL Class", "").strip()
    if "Insoluble" in esol_class:
        failed.append(f"Solubility ({esol_class})")
    elif esol_class:
        passed.append(f"Solubility ({esol_class})")

    # 6. BBB permeability (informational)
    bbb = row.get("BBB permeant", "").strip()
    if bbb == "No":
        passed.append("BBB (No — preferable for antibiotics)")
    else:
        warnings.append(f"BBB permeant ({bbb})")

    # Classify
    hard_fails = [f for f in failed if "Lipinski" in f or "PAINS" in f or "Bioavailability" in f]
    soft_fails = [f for f in failed if f not in hard_fails]

    if len(hard_fails) == 0 and len(soft_fails) == 0:
        classification = "Green"
    elif len(hard_fails) == 0:
        classification = "Yellow"
    elif len(hard_fails) >= 2:
        classification = "Red"
    else:
        # 1 hard fail → Yellow (borderline)
        classification = "Yellow"

    return classification, passed, failed + warnings


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="admet_analysis.py — ADMET analysis of top docking candidates"
    )
    parser.add_argument(
        "--config", default=str(project_path("config.yaml")),
        help="Path to config.yaml"
    )
    parser.add_argument(
        "--top", type=int, default=None,
        help="Analyse top N compounds from ranking (overrides config top_n_final)"
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    config = load_config(config_path)

    output_dir = project_path(config.get("output_dir", "results/"))
    output_dir.mkdir(parents=True, exist_ok=True)

    ranking_csv = output_dir / "ranking.csv"
    compounds_sdf = project_path("data/compounds/compounds.sdf")
    admet_report_csv = output_dir / "admet_report.csv"
    final_candidates_csv = output_dir / "final_candidates.csv"

    top_n = args.top or int(config.get("top_n_final", 10))

    # ── Load ranking ──────────────────────────────────────────────────────────
    if not ranking_csv.exists():
        print(f"Error: {ranking_csv} not found. Run screen_and_rank.py first.")
        sys.exit(1)

    ranked_rows = []
    with open(ranking_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ranked_rows.append(row)

    top_rows = ranked_rows[:top_n]
    print(f"Analysing top {len(top_rows)} compounds from ranking.csv...")

    # ── Build SMILES map from SDF ─────────────────────────────────────────────
    if not compounds_sdf.exists():
        print(f"Error: {compounds_sdf} not found. Run fetch_compounds.py first.")
        sys.exit(1)

    print("Loading SMILES from compounds SDF...")
    smiles_map = build_smiles_map(compounds_sdf)
    print(f"  {len(smiles_map)} SMILES loaded.")

    # ── ADMET analysis loop ───────────────────────────────────────────────────
    admet_results = []

    for i, row in enumerate(top_rows, 1):
        compound_id = row["compound_id"]
        pdb_id = row["pdb_id"]
        dg = float(row["dg_kcal_mol"])

        smiles = smiles_map.get(compound_id)
        if not smiles:
            print(f"  [{i}/{len(top_rows)}] {compound_id}: SMILES not found, skipping.")
            continue

        print(f"  [{i}/{len(top_rows)}] {compound_id} (ΔG={dg:.3f})...")
        admet_row = query_swissadme(smiles)

        if admet_row is None:
            print(f"    SwissADME query failed.")
            admet_results.append({
                "rank": row["rank"],
                "pdb_id": pdb_id,
                "compound_id": compound_id,
                "dg_kcal_mol": dg,
                "smiles": smiles,
                "classification": "N/A",
                "mw": "",
                "logp": "",
                "hbd": "",
                "hba": "",
                "tpsa": "",
                "gi_absorption": "",
                "bbb_permeant": "",
                "pains_alerts": "",
                "lipinski_violations": "",
                "bioavailability_score": "",
                "esol_class": "",
                "passed": "",
                "failed": "",
                "swissadme_error": "API query failed",
            })
            continue

        classification, passed, failed = classify_admet(admet_row, smiles)
        color_map = {"Green": "✓", "Yellow": "~", "Red": "✗", "N/A": "?"}
        print(f"    {color_map.get(classification, '?')} {classification} — "
              f"MW={admet_row.get('MW','')} LogP={admet_row.get('Consensus Log P','')} "
              f"GI={admet_row.get('GI absorption','')} PAINS={admet_row.get('PAINS #alerts','')}")

        admet_results.append({
            "rank": row["rank"],
            "pdb_id": pdb_id,
            "compound_id": compound_id,
            "dg_kcal_mol": dg,
            "smiles": smiles,
            "classification": classification,
            "mw": admet_row.get("MW", ""),
            "logp": admet_row.get("Consensus Log P", ""),
            "hbd": admet_row.get("#H-bond donors", ""),
            "hba": admet_row.get("#H-bond acceptors", ""),
            "tpsa": admet_row.get("TPSA", ""),
            "gi_absorption": admet_row.get("GI absorption", ""),
            "bbb_permeant": admet_row.get("BBB permeant", ""),
            "pains_alerts": admet_row.get("PAINS #alerts", ""),
            "lipinski_violations": admet_row.get("Lipinski #violations", ""),
            "bioavailability_score": admet_row.get("Bioavailability Score", ""),
            "esol_class": admet_row.get("ESOL Class", ""),
            "passed": " | ".join(passed),
            "failed": " | ".join(failed),
            "swissadme_error": "",
        })

        # Rate-limit: 1 second between requests
        if i < len(top_rows):
            time.sleep(1.0)

    # ── Write reports ─────────────────────────────────────────────────────────
    fieldnames = [
        "rank", "pdb_id", "compound_id", "dg_kcal_mol", "smiles",
        "classification", "mw", "logp", "hbd", "hba", "tpsa",
        "gi_absorption", "bbb_permeant", "pains_alerts",
        "lipinski_violations", "bioavailability_score", "esol_class",
        "passed", "failed", "swissadme_error",
    ]

    # Full ADMET report
    with open(admet_report_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(admet_results)

    # Final candidates (Green + Yellow only)
    final = [r for r in admet_results if r["classification"] in ("Green", "Yellow")]
    with open(final_candidates_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final)

    # ── Summary ───────────────────────────────────────────────────────────────
    counts = {
        "Green": sum(1 for r in admet_results if r["classification"] == "Green"),
        "Yellow": sum(1 for r in admet_results if r["classification"] == "Yellow"),
        "Red": sum(1 for r in admet_results if r["classification"] == "Red"),
        "N/A": sum(1 for r in admet_results if r["classification"] == "N/A"),
    }
    print(f"\n{'='*60}")
    print(f"ADMET analysis complete.")
    print(f"  Analysed:          {len(admet_results)} compounds")
    print(f"  Green (priority):  {counts['Green']}")
    print(f"  Yellow (review):   {counts['Yellow']}")
    print(f"  Red (discard):     {counts['Red']}")
    print(f"  N/A (API failed):  {counts['N/A']}")
    print(f"\n  Full report:       {admet_report_csv}")
    print(f"  Final candidates:  {final_candidates_csv}")

    if final:
        print(f"\nFinal candidates ({len(final)}):")
        print(f"  {'Rank':<6} {'Compound':<35} {'ΔG':<10} {'Class':<8} {'MW':<8} {'LogP'}")
        print(f"  {'-'*6} {'-'*35} {'-'*10} {'-'*8} {'-'*8} {'-'*6}")
        for r in final:
            print(f"  {r['rank']:<6} {r['compound_id']:<35} {r['dg_kcal_mol']:<10} "
                  f"{r['classification']:<8} {r['mw']:<8} {r['logp']}")
    else:
        print("\nNo Green or Yellow candidates found. Consider relaxing criteria or fetching more compounds.")


if __name__ == "__main__":
    main()
