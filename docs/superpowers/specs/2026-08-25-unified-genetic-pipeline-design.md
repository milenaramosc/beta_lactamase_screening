# Unified Genetic Pipeline Design

## Goal

Unify the three genetic-algorithm entry points into one shared architecture that supports candidate generation from docking rankings, consensus-ranked molecules, or a protein structure, while preserving the existing script names.

## Requirements

- Keep `scripts/generate_candidate.py`, `scripts/generate_inhibitors_from_protein.py`, and `scripts/genetic_optimize.py` available as user-facing commands.
- Move shared GA behavior into `src/betalactamase_engine/generation/`.
- Support three seed modes:
  - `ranking`: seed from `results/ranking.csv` and `compounds.sdf`.
  - `consensus`: seed from `final_candidates.csv` and `compounds.sdf`.
  - `protein`: seed from known beta-lactamase inhibitors and score with protein-site information.
- Apply selection, mutation, crossover, elite preservation, uniqueness handling, and final chemical filtering through the same engine.
- Include Lipinski RO5 and physicochemical/structural properties in output.
- Make generated SDF files usable as docking input.
- Let `scripts/screen_and_rank.py` accept a ligand SDF path with `--ligands`.
- Preserve traceability from generated candidates to seed molecules and generation operations.

## Architecture

The new architecture keeps the existing `genetic_optimizer.py` module as the shared GA engine and extends it to accept a typed `SeedRecord` list and a configurable fitness strategy. Loader functions convert ranking, consensus, and protein inputs into this shared seed format. Thin CLIs parse legacy arguments and call the same engine.

`screen_and_rank.py` will gain a `--ligands` option that defaults to `data/compounds/compounds.sdf`, allowing generated SDFs to be docked without replacing the main compound library.

## Compatibility

Existing commands continue to work. The old scripts become compatibility wrappers over the unified engine. Outputs remain SDF and CSV, with a richer common schema.

## Validation

- `python -m compileall scripts src`
- `python scripts/generate_candidate.py --help`
- `python scripts/generate_inhibitors_from_protein.py --help`
- `python scripts/genetic_optimize.py --help`
- Short `/tmp` runs for ranking/consensus/protein modes where local inputs exist.
