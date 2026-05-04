#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_DIR / "src"))

GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"

try:
    from betalactamase_engine.scoring.interaction_scoring import (
        score_interactions,
        write_csv,
        write_json,
    )
except ImportError as exc:
    print(f"{RED}Error: Could not import scoring module. {exc}{RESET}")
    sys.exit(1)


def _ensure_output_path(out_path: Path):
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        print(f"{RED}Error: could not create output directory {out_path.parent}: {exc}{RESET}")
        sys.exit(1)


def _print_summary(result, out_path: Path):
    print(f"{GREEN}Interaction scoring completed.{RESET}")
    print(f"Compound: {result.compound_id}")
    print(f"Evaluated site: {result.evaluated_site}")
    if result.distance_to_site_center is None:
        print("Distance to active site center: N/A")
    else:
        print(f"Distance to active site center: {result.distance_to_site_center:.2f} Å")
    print(f"Nearby residue contacts: {result.contact_summary.get('nearby_residue_contact_count', 0)}")
    print(f"Critical residue contacts: {result.contact_summary.get('critical_residue_contact_count', 0)}")
    print(f"Metal contacts: {result.contact_summary.get('metal_contact_count', 0)}")
    print(f"Interaction score: {result.interaction_score:.2f}")
    print(f"Pose classification: {result.pose_classification}")
    print(f"Output: {out_path}")
    if result.warnings:
        print(f"{YELLOW}Warnings:{RESET}")
        for warning in result.warnings:
            print(f" - {warning}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze docking pose interactions against the active site."
    )
    parser.add_argument("--target-profile", required=True, help="Path to target_profile.json")
    parser.add_argument("--protein", required=True, help="Path to protein PDB file")
    parser.add_argument("--ligand-pose", required=True, help="Path to ligand pose file")
    parser.add_argument("--out", required=True, help="Path to write analysis output")
    parser.add_argument(
        "--compound-id",
        help="Compound identifier (defaults to ligand file name)",
    )
    parser.add_argument(
        "--site",
        default="active_site",
        choices=["active_site"],
        help="Site identifier (default: active_site)",
    )
    parser.add_argument(
        "--distance-threshold",
        type=float,
        default=4.0,
        help="Distance threshold for residue contacts in Å (default: 4.0)",
    )
    parser.add_argument(
        "--center-threshold",
        type=float,
        default=8.0,
        help="Distance threshold for site center in Å (default: 8.0)",
    )
    parser.add_argument(
        "--metal-threshold",
        type=float,
        default=3.0,
        help="Distance threshold for metal contacts in Å (default: 3.0)",
    )
    parser.add_argument(
        "--format",
        choices=["json", "csv"],
        default="json",
        help="Output format (default: json)",
    )
    parser.add_argument("--verbose", action="store_true", help="Show summary output")

    args = parser.parse_args()

    ligand_path = Path(args.ligand_pose)
    compound_id = args.compound_id or ligand_path.stem

    out_path = Path(args.out)
    _ensure_output_path(out_path)

    result = score_interactions(
        target_profile_path=Path(args.target_profile),
        protein_path=Path(args.protein),
        ligand_pose_path=ligand_path,
        compound_id=compound_id,
        site=args.site,
        distance_threshold=args.distance_threshold,
        center_threshold=args.center_threshold,
        metal_threshold=args.metal_threshold,
    )

    if args.format == "json":
        write_json(result, out_path)
    else:
        write_csv(result, out_path)

    if args.verbose:
        _print_summary(result, out_path)


if __name__ == "__main__":
    main()
