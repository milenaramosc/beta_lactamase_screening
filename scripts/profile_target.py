#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

# Add src to sys.path to allow imports
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_DIR / "src"))

# ── Cores ANSI ────────────────────────────────────────────────────────────────

GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"

try:
    from betalactamase_engine.target_profiler import TargetProfiler
except ImportError as e:
    print(f"{RED}Error: Could not import TargetProfiler. {e}{RESET}")
    sys.exit(1)

def _ensure_output_path(out_path: Path):
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        print(f"{RED}Error: could not create output directory {out_path.parent}: {exc}{RESET}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Generate a target profile for a beta-lactamase PDB.")
    parser.add_argument("--pdb", required=True, help="Path to the PDB file.")
    parser.add_argument("--out", required=True, help="Path to the output JSON file.")
    parser.add_argument("--radius", type=float, default=8.0, help="Radius for active site analysis (default: 8.0).")
    parser.add_argument("--chain", help="Chain ID to analyze (optional).")
    parser.add_argument("--active-site", help="User-provided active site residues (CHAIN:RESNAME:RESID,...)")
    parser.add_argument("--no-interactive", action="store_true", help="Disable interactive prompts.")
    parser.add_argument("--verbose", action="store_true", help="Display summary in terminal.")

    args = parser.parse_args()

    pdb_path = Path(args.pdb)
    out_path = Path(args.out)

    if not pdb_path.exists():
        print(f"{RED}Error: PDB file not found: {pdb_path}{RESET}")
        sys.exit(1)

    if pdb_path.stat().st_size == 0:
        print(f"{RED}Error: PDB file is empty: {pdb_path}{RESET}")
        sys.exit(1)

    _ensure_output_path(out_path)

    profiler = TargetProfiler(pdb_path, radius=args.radius, chain=args.chain)

    if args.active_site:
        result = profiler.set_user_active_site(args.active_site)
        if result.get("format_errors") and "empty input" not in result["format_errors"]:
            print(f"{YELLOW}Invalid active site format: {', '.join(result['format_errors'])}{RESET}")
        if result.get("missing"):
            profiler.record_missing_active_site_entries(result["missing"])
        if not result.get("valid"):
            print(f"{YELLOW}No valid active site residues found; falling back to automatic detection.{RESET}")
    else:
        site_records = TargetProfiler.parse_site_records_from_file(pdb_path)
        if not site_records and not args.no_interactive:
            print("O arquivo PDB nao possui registros SITE anotados.")
            print("Voce sabe informar os residuos do sitio ativo?")
            print("Digite no formato CHAIN:RESNAME:RESID separados por virgula.")
            print("Exemplo: A:SER:70,A:LYS:73,A:GLU:166")
            user_input = input("Ou pressione ENTER para tentar deteccao automatica: ").strip()
            if user_input:
                result = profiler.set_user_active_site(user_input)
                invalid_items = [item for item in result.get("format_errors", []) if item != "empty input"]
                if invalid_items or not result.get("valid"):
                    if result.get("missing"):
                        profiler.record_missing_active_site_entries(result["missing"])
                    if invalid_items:
                        print(f"{YELLOW}Residuos invalidos: {', '.join(invalid_items)}{RESET}")
                    retry = input("Deseja tentar novamente? (s/N): ").strip().lower()
                    if retry == "s":
                        retry_input = input("Informe os residuos no formato CHAIN:RESNAME:RESID: ").strip()
                        if retry_input:
                            result = profiler.set_user_active_site(retry_input)
                            if result.get("missing"):
                                profiler.record_missing_active_site_entries(result["missing"])
                        else:
                            print("Usando deteccao automatica.")
                else:
                    if result.get("missing"):
                        profiler.record_missing_active_site_entries(result["missing"])

    profile_data = profiler.run()
    profiler.save_json(out_path)

    if args.verbose:
        target = profile_data.target
        classification = profile_data.classification
        active_site = profile_data.active_site

        print(f"{GREEN}\nTarget profile generated successfully.{RESET}")
        print(f"PDB: {target.get('pdb_name')}.pdb")
        print(f"Predicted class: {classification.get('predicted_class')}")
        print(f"Active site method: {active_site.get('detection_method')}")

        center = active_site.get('center', [0.0, 0.0, 0.0])
        print(f"Active site center: {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")

        metals = profile_data.metals
        print(f"Metals detected: {len(metals)}")

        catalytic = active_site.get('candidate_catalytic_residues', [])
        print(f"Candidate catalytic residues: {len(catalytic)}")

        print(f"Output: {out_path}")

        if profile_data.warnings:
            print(f"{YELLOW}\nWarnings:{RESET}")
            for warning in profile_data.warnings:
                print(f" - {warning}")

if __name__ == "__main__":
    main()
