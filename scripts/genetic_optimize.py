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
    from betalactamase_engine.generation.genetic_optimizer import (
        run_genetic_optimization,
    )
except ImportError as exc:
    print(f"{RED}Error: Could not import genetic optimizer module. {exc}{RESET}")
    sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Genetic optimization guided by consensus final_score"
    )
    parser.add_argument("--final-candidates", required=True)
    parser.add_argument("--compounds", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--out-csv")
    parser.add_argument("--top-n-seeds", type=int, default=10)
    parser.add_argument("--min-final-score", type=float, default=0.0)
    parser.add_argument("--generations", type=int, default=10)
    parser.add_argument("--population-size", type=int, default=30)
    parser.add_argument("--mutation-rate", type=float, default=0.25)
    parser.add_argument("--crossover-rate", type=float, default=0.50)
    parser.add_argument("--elite-size", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--run-id")
    parser.add_argument("--max-molecular-weight", type=float, default=650.0)
    parser.add_argument("--max-tpsa", type=float, default=250.0)
    parser.add_argument("--max-hbd", type=int, default=8)
    parser.add_argument("--max-hba", type=int, default=15)
    parser.add_argument("--min-qed", type=float, default=0.05)
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    final_candidates_path = Path(args.final_candidates)
    compounds_path = Path(args.compounds)
    out_path = Path(args.out)
    out_csv_path = Path(args.out_csv) if args.out_csv else None

    if not final_candidates_path.exists():
        raise SystemExit("final_candidates.csv not found")
    if not compounds_path.exists():
        raise SystemExit("compounds.sdf not found")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_csv_path:
        out_csv_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        print(
            f"{YELLOW}Warning: --out-csv not provided; skipping CSV summary output.{RESET}"
        )

    try:
        summary = run_genetic_optimization(
            final_candidates_path=final_candidates_path,
            compounds_path=compounds_path,
            out_sdf=out_path,
            out_csv=out_csv_path,
            top_n_seeds=args.top_n_seeds,
            min_final_score=args.min_final_score,
            generations=args.generations,
            population_size=args.population_size,
            mutation_rate=args.mutation_rate,
            crossover_rate=args.crossover_rate,
            elite_size=args.elite_size,
            seed=args.seed,
            run_id=args.run_id,
            max_molecular_weight=args.max_molecular_weight,
            max_tpsa=args.max_tpsa,
            max_hbd=args.max_hbd,
            max_hba=args.max_hba,
            min_qed=args.min_qed,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(
            f"{RED}Arquivo de candidatos inválido.{RESET} Esta etapa requer o CSV "
            "gerado por scripts/consensus_score.py contendo compound_id, "
            "final_score e final_classification."
        )
        print(f"{YELLOW}Detalhe:{RESET} {exc}")
        sys.exit(1)

    if args.verbose:
        print(f"{GREEN}Genetic optimization completed.{RESET}")
        print(f"Final candidates input: {final_candidates_path}")
        print(f"Seed molecules selected: {summary.seed_selected}")
        print(f"Seed molecules found in SDF: {summary.seed_mapped}")
        print(f"Generations: {args.generations}")
        print(f"Population size: {args.population_size}")
        print(f"Generated valid molecules: {summary.generated_valid}")
        print(f"Unique generated molecules: {summary.generated_unique}")
        print(f"Filtered molecules: {summary.filtered}")
        print(f"Output SDF: {out_path}")
        if out_csv_path:
            print(f"Output CSV: {out_csv_path}")
            print(f"Filtered CSV: {summary.filtered_csv}")
        else:
            print("Output CSV: not provided")
        if summary.warnings:
            print("Warnings:")
            for warning in summary.warnings:
                print(f" - {warning}")

    print(
        "Os candidatos gerados devem retornar ao ciclo de docking e consensus scoring."
    )


if __name__ == "__main__":
    main()
