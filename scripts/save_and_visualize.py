#!/usr/bin/env python3
"""
save_and_visualize.py — Phase 6: Save complexes and generate visualizations

Usage:
    python save_and_visualize.py [--config CONFIG]

For each final candidate in results/final_candidates.csv:
  1. Combines protein + docking pose into results/complexes/PDBID_compound.pdb
  2. Generates interactive 3D HTML with py3Dmol
  3. Generates static PNG with PyMOL (if available)
  4. Creates results/final_report.html — full summary table

Skips PNG generation gracefully if PyMOL is not installed.
"""

import argparse
import csv
import html as html_lib
import json
import re
import sys
from pathlib import Path

import py3Dmol
import yaml

# ── Path helpers ──────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent


def project_path(p: str) -> Path:
    return PROJECT_DIR / p


# ── Config ────────────────────────────────────────────────────────────────────

def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# ── PDBQT → PDB (first model only) ───────────────────────────────────────────

def pdbqt_to_pdb_lines(pdbqt_path: Path) -> list[str]:
    """
    Extract the first MODEL from a PDBQT file and convert to basic PDB lines.
    Strips AutoDock-specific columns (partial charges, atom type at end).
    """
    in_model = False
    pdb_lines = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                break  # only first model
            if not in_model:
                continue
            if line.startswith(("ATOM", "HETATM")):
                # PDB ATOM record is 80 chars; PDBQT adds extra cols at end
                # Truncate to standard PDB columns (keep first 66 chars)
                pdb_line = line[:66].rstrip()
                # Change residue name to LIG for ligand atoms
                pdb_line = pdb_line[:17] + "LIG" + pdb_line[20:]
                pdb_lines.append(pdb_line)
            elif line.startswith(("REMARK", "BRANCH", "ENDBRANCH", "ROOT", "ENDROOT", "TORSDOF")):
                continue  # skip vina-specific records
    return pdb_lines


# ── Complex PDB generation ────────────────────────────────────────────────────

def save_complex(
    protein_pdb_path: Path,
    pose_pdbqt_path: Path,
    out_path: Path,
    compound_id: str,
) -> Path:
    """
    Combine protein PDB + ligand pose into a single complex PDB file.
    """
    # Read protein (skip HETATM records that aren't part of chain to avoid confusion)
    with open(protein_pdb_path) as f:
        protein_lines = [
            line.rstrip()
            for line in f
            if line.startswith(("ATOM", "TER", "SEQRES"))
        ]

    # Convert ligand PDBQT to PDB lines
    ligand_lines = pdbqt_to_pdb_lines(pose_pdbqt_path)

    with open(out_path, "w") as f:
        f.write(f"REMARK Complex: {protein_pdb_path.stem} + {compound_id}\n")
        for line in protein_lines:
            f.write(line + "\n")
        f.write("TER\n")
        f.write(f"REMARK Ligand: {compound_id}\n")
        for line in ligand_lines:
            # Mark ligand as HETATM in chain L
            if line.startswith("ATOM"):
                line = "HETATM" + line[6:21] + "L" + line[22:]
            f.write(line + "\n")
        f.write("END\n")

    return out_path


# ── py3Dmol HTML visualization ────────────────────────────────────────────────

def generate_html_visualization(
    complex_pdb_path: Path,
    compound_id: str,
    out_html_path: Path,
    dg: float,
    classification: str,
) -> Path:
    """
    Generate a standalone interactive 3D HTML visualization using py3Dmol.
    The protein is shown as cartoon (rainbow), the ligand as sticks.
    """
    with open(complex_pdb_path) as f:
        complex_pdb = f.read()

    view = py3Dmol.view(width=800, height=600)
    view.addModel(complex_pdb, "pdb")

    # Protein: cartoon, colored by chain (rainbow spectrum)
    view.setStyle(
        {"resn": ["LIG"]},
        {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}},
    )
    view.setStyle(
        {"not": {"resn": ["LIG"]}},
        {"cartoon": {"color": "spectrum", "opacity": 0.85}},
    )

    # Highlight binding pocket residues within 5 Å of ligand
    view.addSurface(
        py3Dmol.SES,
        {"opacity": 0.1, "color": "lightblue"},
        {"not": {"resn": ["LIG"]}},
    )

    view.zoomTo({"resn": ["LIG"]})
    view.setBackgroundColor("white")

    inner_html = view._make_html()

    # Wrap in a full page with title and metadata
    color_badge = {
        "Green": "#28a745",
        "Yellow": "#ffc107",
        "Red": "#dc3545",
    }.get(classification, "#6c757d")

    page_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Beta-Lactamase Inhibitor Candidate: {html_lib.escape(compound_id)}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; }}
    h1 {{ color: #333; }}
    .meta {{ background: white; padding: 15px; border-radius: 8px; margin-bottom: 20px;
              box-shadow: 0 1px 3px rgba(0,0,0,0.1); display: inline-block; }}
    .badge {{ display: inline-block; padding: 4px 12px; border-radius: 12px;
               color: white; font-weight: bold; background: {color_badge}; }}
    .viewer-container {{ background: white; padding: 10px; border-radius: 8px;
                          box-shadow: 0 1px 3px rgba(0,0,0,0.1); display: inline-block; }}
  </style>
</head>
<body>
  <h1>Beta-Lactamase Inhibitor Candidate</h1>
  <div class="meta">
    <p><strong>Compound:</strong> {html_lib.escape(compound_id)}</p>
    <p><strong>Binding affinity (ΔG):</strong> {dg:.3f} kcal/mol</p>
    <p><strong>ADMET classification:</strong> <span class="badge">{html_lib.escape(classification)}</span></p>
  </div>
  <br>
  <div class="viewer-container">
    {inner_html}
  </div>
  <p style="color:#666;font-size:12px;margin-top:10px;">
    Protein shown as rainbow cartoon. Ligand shown as green sticks.
    Rotate with left-click, zoom with scroll, translate with right-click.
  </p>
</body>
</html>"""

    out_html_path.write_text(page_html)
    return out_html_path


# ── PyMOL PNG (optional) ──────────────────────────────────────────────────────

def generate_pymol_png(
    complex_pdb_path: Path,
    compound_id: str,
    out_png_path: Path,
) -> bool:
    """
    Generate a static PNG using PyMOL headless. Returns True if successful.
    """
    try:
        import pymol
        from pymol import cmd

        cmd.reinitialize()
        cmd.load(str(complex_pdb_path), "complex")
        cmd.hide("everything")
        cmd.show("cartoon", "complex and not resn LIG")
        cmd.show("sticks", "complex and resn LIG")
        # color protein by residue index along chain (rainbow)
        cmd.spectrum("count", "rainbow", "complex and not resn LIG and name CA")
        cmd.color("green", "complex and resn LIG")
        cmd.zoom("resn LIG", buffer=8)
        cmd.bg_color("white")
        cmd.ray(800, 600)
        cmd.png(str(out_png_path), dpi=150, quiet=1)
        return True
    except ImportError:
        return False
    except Exception as e:
        print(f"    PyMOL error for {compound_id}: {e}", file=sys.stderr)
        return False


# ── Final HTML report ─────────────────────────────────────────────────────────

def generate_final_report(
    candidates: list[dict],
    visualization_dir: Path,
    out_path: Path,
):
    """
    Generate results/final_report.html — a summary table of all final candidates
    with ΔG, ADMET classification, and links to 3D visualizations.
    """
    rows_html = ""
    for r in candidates:
        compound_id = r["compound_id"]
        viz_html = visualization_dir / f"{r['pdb_id']}_{compound_id}.html"
        viz_link = (
            f'<a href="visualizations/{html_lib.escape(str(viz_html.name))}" target="_blank">View 3D</a>'
            if viz_html.exists()
            else "N/A"
        )
        color_badge = {
            "Green": "#28a745",
            "Yellow": "#ffc107",
            "Red": "#dc3545",
        }.get(r["classification"], "#6c757d")

        rows_html += f"""
        <tr>
          <td>{html_lib.escape(str(r['rank']))}</td>
          <td>{html_lib.escape(r['pdb_id'])}</td>
          <td style="font-family:monospace;font-size:0.9em">{html_lib.escape(compound_id)}</td>
          <td>{float(r['dg_kcal_mol']):.3f}</td>
          <td>{html_lib.escape(str(r.get('mw', '')))}</td>
          <td>{html_lib.escape(str(r.get('logp', '')))}</td>
          <td>{html_lib.escape(str(r.get('gi_absorption', '')))}</td>
          <td>{html_lib.escape(str(r.get('bbb_permeant', '')))}</td>
          <td>{html_lib.escape(str(r.get('pains_alerts', '')))}</td>
          <td>
            <span style="background:{color_badge};color:white;padding:3px 10px;
                         border-radius:10px;font-weight:bold">
              {html_lib.escape(r['classification'])}
            </span>
          </td>
          <td>{viz_link}</td>
        </tr>"""

    report_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Beta-Lactamase Virtual Screening — Final Report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 30px; background: #f8f9fa; }}
    h1 {{ color: #2c3e50; }}
    h2 {{ color: #555; font-size: 1.1em; font-weight: normal; margin-top: -10px; }}
    table {{ border-collapse: collapse; width: 100%; background: white;
              box-shadow: 0 1px 3px rgba(0,0,0,0.1); border-radius: 8px; overflow: hidden; }}
    th {{ background: #2c3e50; color: white; padding: 10px 14px; text-align: left; font-size: 0.9em; }}
    td {{ padding: 9px 14px; border-bottom: 1px solid #eee; font-size: 0.9em; }}
    tr:hover {{ background: #f1f3f5; }}
    .legend {{ margin-top: 20px; font-size: 0.85em; color: #666; }}
    .green {{ color: #28a745; font-weight: bold; }}
    .yellow {{ color: #856404; font-weight: bold; }}
    .red {{ color: #dc3545; font-weight: bold; }}
  </style>
</head>
<body>
  <h1>Beta-Lactamase Inhibitor Candidates</h1>
  <h2>Virtual Screening Pipeline — Final Report</h2>
  <table>
    <thead>
      <tr>
        <th>Rank</th>
        <th>Target (PDB)</th>
        <th>Compound ID</th>
        <th>ΔG (kcal/mol)</th>
        <th>MW (Da)</th>
        <th>LogP</th>
        <th>GI Absorption</th>
        <th>BBB</th>
        <th>PAINS</th>
        <th>ADMET</th>
        <th>3D View</th>
      </tr>
    </thead>
    <tbody>
      {rows_html}
    </tbody>
  </table>
  <div class="legend">
    <strong>Legend:</strong>
    <span class="green">■ Green</span> — passes all ADMET criteria, priority candidate for wet lab validation. &nbsp;
    <span class="yellow">■ Yellow</span> — passes most criteria, requires manual review. &nbsp;
    <span class="red">■ Red</span> — fails key criteria, discarded. <br>
    ΔG: binding free energy (more negative = stronger predicted binding affinity). <br>
    BBB: blood-brain barrier permeability (negative preferred for antibiotics targeting bacteria).
  </div>
</body>
</html>"""

    out_path.write_text(report_html)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="save_and_visualize.py — Save complexes and generate visualizations"
    )
    parser.add_argument(
        "--config", default=str(project_path("config.yaml")),
        help="Path to config.yaml"
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    config = load_config(config_path)

    output_dir = project_path(config.get("output_dir", "results/"))
    output_dir.mkdir(parents=True, exist_ok=True)

    final_candidates_csv = output_dir / "final_candidates.csv"
    structures_dir = project_path("data/structures")
    poses_dir = output_dir / "top50_poses"

    complexes_dir = output_dir / "complexes"
    viz_dir = output_dir / "visualizations"
    complexes_dir.mkdir(exist_ok=True)
    viz_dir.mkdir(exist_ok=True)

    final_report_html = output_dir / "final_report.html"

    # ── Load final candidates ─────────────────────────────────────────────────
    if not final_candidates_csv.exists():
        print(f"Error: {final_candidates_csv} not found. Run admet_analysis.py first.")
        sys.exit(1)

    candidates = []
    with open(final_candidates_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append(row)

    if not candidates:
        print("No final candidates found. Nothing to visualize.")
        sys.exit(0)

    print(f"Processing {len(candidates)} final candidates...")

    for i, cand in enumerate(candidates, 1):
        pdb_id = cand["pdb_id"]
        compound_id = cand["compound_id"]
        dg = float(cand["dg_kcal_mol"])
        classification = cand.get("classification", "N/A")

        print(f"\n[{i}/{len(candidates)}] {pdb_id} + {compound_id} ({classification}, ΔG={dg:.3f})")

        protein_pdb = structures_dir / f"{pdb_id}.pdb"
        pose_pdbqt = poses_dir / f"{pdb_id}__{compound_id}.pdbqt"

        if not protein_pdb.exists():
            print(f"  Warning: protein PDB {protein_pdb} not found, skipping.")
            continue
        if not pose_pdbqt.exists():
            print(f"  Warning: pose PDBQT {pose_pdbqt} not found, skipping.")
            continue

        # 1. Save complex PDB
        complex_path = complexes_dir / f"{pdb_id}_{compound_id}.pdb"
        save_complex(protein_pdb, pose_pdbqt, complex_path, compound_id)
        print(f"  Complex saved: {complex_path.name}")

        # 2. Generate interactive HTML
        html_out = viz_dir / f"{pdb_id}_{compound_id}.html"
        generate_html_visualization(complex_path, compound_id, html_out, dg, classification)
        print(f"  HTML viewer: {html_out.name}")

        # 3. Generate static PNG (optional)
        png_out = viz_dir / f"{pdb_id}_{compound_id}.png"
        png_ok = generate_pymol_png(complex_path, compound_id, png_out)
        if png_ok:
            print(f"  PNG image: {png_out.name}")
        else:
            print(f"  PNG skipped (PyMOL not available).")

    # ── Final report ──────────────────────────────────────────────────────────
    print(f"\nGenerating final report...")
    generate_final_report(candidates, viz_dir, final_report_html)
    print(f"  Report: {final_report_html}")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"Visualization complete.")
    print(f"  Complexes:     {complexes_dir}")
    print(f"  Visualizations: {viz_dir}")
    print(f"  Final report:  {final_report_html}")
    print(f"\nOpen {final_report_html} in a browser to explore results.")


if __name__ == "__main__":
    main()
