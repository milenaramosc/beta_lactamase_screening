#!/usr/bin/env python3
"""
Generate de novo candidates from docking-ranked seed compounds.

Compatibility wrapper for the unified genetic optimizer. It reads
`results/ranking.csv` plus `data/compounds/compounds.sdf` and writes
`results/generated_candidates.sdf` and `results/generated_candidates.csv`.
"""

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
    from betalactamase_engine.generation import (
        load_ranking_seeds,
        run_genetic_optimization_from_seeds,
    )
except ImportError as exc:
    print(f"{RED}Error: Could not import genetic optimizer module. {exc}{RESET}")
    sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate candidates with the unified AG from docking ranking"
    )
    parser.add_argument("--ranking", default=str(PROJECT_DIR / "results" / "ranking.csv"))
    parser.add_argument("--compounds", default=str(PROJECT_DIR / "data" / "compounds" / "compounds.sdf"))
    parser.add_argument("--out", default=str(PROJECT_DIR / "results" / "generated_candidates.sdf"))
    parser.add_argument("--out-csv", default=str(PROJECT_DIR / "results" / "generated_candidates.csv"))
    parser.add_argument("--generations", type=int, default=30)
    parser.add_argument("--population", type=int, default=50)
    parser.add_argument("--top-seed", type=int, default=20, dest="top_seed")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--mutation-rate", type=float, default=0.25)
    parser.add_argument("--crossover-rate", type=float, default=0.50)
    parser.add_argument("--elite-size", type=int, default=5)
    parser.add_argument("--max-molecular-weight", type=float, default=500.0)
    parser.add_argument("--max-logp", type=float, default=5.0)
    parser.add_argument("--max-tpsa", type=float, default=250.0)
    parser.add_argument("--max-hbd", type=int, default=5)
    parser.add_argument("--max-hba", type=int, default=10)
    parser.add_argument("--min-qed", type=float, default=0.05)
    parser.add_argument("--top-out", type=int, default=5, dest="top_out")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    try:
        seeds, warnings = load_ranking_seeds(
            Path(args.ranking),
            Path(args.compounds),
            top_n=args.top_seed,
        )
        summary = run_genetic_optimization_from_seeds(
            seeds=seeds,
            out_sdf=Path(args.out),
            out_csv=Path(args.out_csv),
            generations=args.generations,
            population_size=args.population,
            mutation_rate=args.mutation_rate,
            crossover_rate=args.crossover_rate,
            elite_size=args.elite_size,
            seed=args.seed,
            max_molecular_weight=args.max_molecular_weight,
            max_logp=args.max_logp,
            max_tpsa=args.max_tpsa,
            max_hbd=args.max_hbd,
            max_hba=args.max_hba,
            min_qed=args.min_qed,
            initial_warnings=warnings,
            output_limit=args.top_out,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"{RED}Erro ao gerar candidatos:{RESET} {exc}")
        sys.exit(1)

    print(f"{GREEN}AG concluido.{RESET}")
    print(f"Sementes selecionadas: {summary.seed_selected}")
    print(f"Moléculas válidas geradas: {summary.generated_valid}")
    print(f"Moléculas filtradas: {summary.filtered}")
    print(f"SDF: {Path(args.out)}")
    print(f"CSV: {Path(args.out_csv)}")
    if summary.filtered_csv:
        print(f"CSV filtradas: {summary.filtered_csv}")
    if args.verbose and summary.warnings:
        print(f"{YELLOW}Warnings:{RESET}")
        for warning in summary.warnings:
            print(f" - {warning}")
    print("Próxima etapa: dockar o SDF gerado com scripts/screen_and_rank.py --ligands.")


if __name__ == "__main__":
    main()
