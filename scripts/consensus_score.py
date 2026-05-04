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
    from betalactamase_engine.scoring.consensus_score import (
        compute_consensus,
        load_admet,
        load_interactions,
        load_ranking,
        sort_consensus_rows,
        summarize_rows,
        write_consensus_csv,
        write_summary_json,
    )
except ImportError as exc:
    print(f"{RED}Error: Could not import consensus scoring module. {exc}{RESET}")
    sys.exit(1)


def _print_summary(
    ranking_path: Path,
    interactions_path: Path,
    admet_path: Path,
    out_path: Path,
    counts,
):
    interactions_label = str(interactions_path) if interactions_path else "not provided"
    admet_label = str(admet_path) if admet_path else "not provided"
    print(f"{GREEN}Consensus scoring completed.{RESET}")
    print(f"Docking ranking: {ranking_path}")
    print(f"Interactions: {interactions_label}")
    print(f"ADMET: {admet_label}")
    print(f"Total candidates: {counts.get('total', 0)}")
    print(f"With interaction scores: {counts.get('with_interactions', 0)}")
    print(f"With ADMET: {counts.get('with_admet', 0)}")
    print(f"High priority candidates: {counts.get('high', 0)}")
    print(f"Medium priority candidates: {counts.get('medium', 0)}")
    print(f"Low priority candidates: {counts.get('low', 0)}")
    print(f"Deprioritized candidates: {counts.get('deprioritized', 0)}")
    print(f"Output: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a multiobjective consensus ranking."
    )
    parser.add_argument("--ranking", required=True, help="Path to docking ranking CSV")
    parser.add_argument("--out", required=True, help="Path to final CSV output")
    parser.add_argument("--interactions", help="Path to interaction scoring CSV/JSON/dir")
    parser.add_argument("--admet", help="Path to ADMET CSV")
    parser.add_argument("--compound-column", help="Compound id column in ranking")
    parser.add_argument("--vina-column", help="Docking score column in ranking")
    parser.add_argument("--top-n", type=int, help="Limit number of results")
    parser.add_argument(
        "--min-interaction-score",
        type=float,
        help="Filter below this interaction score (penalizes if missing)",
    )
    parser.add_argument("--verbose", action="store_true", help="Show summary")

    args = parser.parse_args()

    ranking_path = Path(args.ranking)
    out_path = Path(args.out)
    interactions_path = Path(args.interactions) if args.interactions else None
    admet_path = Path(args.admet) if args.admet else None

    try:
        ranking_rows, compound_col, vina_col, columns = load_ranking(
            ranking_path,
            compound_col_override=args.compound_column,
            vina_col_override=args.vina_column,
        )

        interactions, interaction_warnings, interactions_available, _ = load_interactions(
            interactions_path
        )
        interactions_provided = interactions_path is not None
        admet, admet_warnings, admet_available = load_admet(admet_path)
        admet_provided = admet_path is not None

        output_rows, counts, weights, warnings = compute_consensus(
            ranking_rows,
            compound_col,
            vina_col,
            interactions,
            interactions_available,
            interactions_provided,
            admet,
            admet_available,
            admet_provided,
            min_interaction_score=args.min_interaction_score,
        )

        warnings.extend(interaction_warnings)
        warnings.extend(admet_warnings)

        output_rows = sort_consensus_rows(output_rows)

        if args.top_n:
            output_rows = output_rows[: args.top_n]

        counts = summarize_rows(
            output_rows,
            interactions_available,
            admet_available,
            interactions_provided,
            admet_provided,
        )

        write_consensus_csv(output_rows, out_path, columns)
        write_summary_json(out_path, counts, weights, warnings)

        if args.verbose:
            _print_summary(
                ranking_path,
                interactions_path,
                admet_path,
                out_path,
                counts,
            )
            if warnings:
                print(f"{YELLOW}Warnings:{RESET}")
                for warning in warnings:
                    print(f" - {warning}")
    except (FileNotFoundError, ValueError) as exc:
        print(f"{RED}Error: {exc}{RESET}")
        sys.exit(1)


if __name__ == "__main__":
    main()
