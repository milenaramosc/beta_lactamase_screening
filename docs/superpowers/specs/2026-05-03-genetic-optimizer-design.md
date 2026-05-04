# Genetic Optimizer (Etapa 4) — Design Spec

## Goal

Implement a new genetic optimization step that uses `final_score` from `final_candidates.csv` to drive seed selection, parent choice, and elite preservation. Generated molecules are **preliminary** and must re-enter the pipeline (preparation → docking → interaction scoring → consensus scoring). Existing scripts must remain intact.

## Scope

In scope:
- New module `src/betalactamase_engine/generation/genetic_optimizer.py`.
- New CLI `scripts/genetic_optimize.py`.
- Documentation `docs/genetic_optimizer.md` (user-facing).
- Minimal, robust genetic operators with RDKit.

Out of scope:
- Docking, ADMET, alosteric sites, molecular dynamics, pipeline refactor, removal of existing scripts.

## Architecture

### Module: `genetic_optimizer.py`

Responsibilities:
- Parse and validate `final_candidates.csv`.
- Rank and select seeds using `final_score` and `final_classification`.
- Map seed compound IDs to molecules in `compounds.sdf`.
- Generate initial population (seeds + simple variants).
- Run genetic optimization with reproducibility.
- Assign fitness:
  - Seeds inherit `final_score` (real fitness).
  - New molecules receive `inherited_or_estimated_fitness` (preliminary).
- Deduplicate by canonical SMILES.
- Write SDF + CSV with required properties.
- Emit warnings and summary.

Key functions (proposed):
- `load_final_candidates(path)`
- `select_seeds(rows, top_n, min_score)`
- `load_compound_library(sdf_path)`
- `map_seeds_to_molecules(seed_rows, library)`
- `build_initial_population(seeds, population_size, rng)`
- `run_genetic_optimization(population, generations, rates, rng)`
- `evaluate_fitness(mol, seed_refs)`
- `write_sdf(...)` and `write_summary_csv(...)`

### CLI: `scripts/genetic_optimize.py`

Responsibilities:
- Parse CLI args and validate paths.
- Call module functions.
- Print verbose summary and next-step guidance.

## Data Flow

1. Read `final_candidates.csv` and validate columns: `compound_id`, `final_score`, `final_classification`.
2. Sort by `final_score` desc.
3. Select seeds:
   - Validate `final_score` is numeric; drop rows with invalid/missing score and warn.
   - Normalize `final_score` to 0–1 for filtering and selection (same rules as fitness normalization).
   - If all normalized scores equal (0.5), bypass `--min-final-score` filtering and keep all valid rows.
   - Otherwise, filter by normalized `final_score >= --min-final-score` (default 0.0 keeps all valid rows).
   - If no rows remain after filtering, fail with a clear error.
   - Classification priority order:
     - `high_priority_candidate`
     - `medium_priority_candidate`
     - `low_priority_candidate`
     - `deprioritized_candidate`
   - Any missing/unknown classification is treated as `low_priority_candidate` (do not drop).
   - Within each class, order by normalized `final_score` desc; tie-break by `compound_id` ascending.
   - Take up to `--top-n-seeds`; include deprioritized only if needed to reach `--top-n-seeds`.
4. Load `compounds.sdf` and map seeds by ID using these property keys (in order):
   - `compound_id`
   - `chembl_id`
   - `molecule_id`
   - `ChEMBL ID`
   - `ID`
   - `name`
   - `_Name`
   If multiple SDF entries match the same `compound_id`, pick the first and warn about duplicates.
   Seeds that cannot be mapped are dropped with a warning; continue if at least one seed maps.
5. Build initial population from mapped seeds; if short, create conservative variants.
   - If seeds > population_size: keep top seeds by the same seed-selection order (class priority then normalized final_score).
   - If seeds < population_size: generate variants per seed in round-robin until population_size is reached.
   - Variants per seed max: 3; if still short, allow duplicates with a warning.
   - Conservative variants use the same mutation operator list with a single successful mutation per variant; if mutation fails, clone the seed.
6. Run GA for `--generations` with selection weighted by normalized fitness, mutation/crossover rates, elite preservation.
   - Deduplicate by canonical SMILES at the end of each generation; keep the candidate with highest fitness (tie-break: seed > generated, then shorter SMILES). Discard invalid molecules.
   - After dedup/validation, if population shrinks below `--population-size`, refill with best remaining candidates or clones of elite (warn on refill).
     - "Best remaining candidates" refers to the current generation's non-elite pool sorted by fitness.
7. Final dedup by canonical SMILES before output; discard invalid molecules.
8. Write SDF and CSV with required properties and notes.
9. Print reminder that generated molecules must re-enter docking and scoring.

## Genetic Operators

- Mutation: small substitutions using this SMARTS list (max 8 attempts per mutation) with RDKit AllChem.ReactionFromSmarts:
  - "[c:1][F]>>[c:1][Cl]"
  - "[c:1][Cl]>>[c:1][F]"
  - "[c:1][Br]>>[c:1][Cl]"
  - "[c:1][CH3]>>[c:1][OH]"
  - "[c:1][OH]>>[c:1][NH2]"
  - "[c:1][NH2]>>[c:1][CH3]"
  - "[C:1](=O)[OH]>>[C:1](=O)[NH2]"
  Must sanitize and validate.
- Crossover: BRICS when possible (RDKit Chem.BRICS). Fallback: return the better parent unchanged (no fragment recombination) when BRICS fails (max 8 attempts). Must sanitize and validate.
- Always discard invalid molecules; never save invalid outputs.

## Fitness Strategy

- Seeds: use `final_score` as real fitness, normalized to 0–1 if needed:
  - If any `final_score` < 0 or > 1.0, normalize all by min/max to 0–1.
  - If all scores equal, use 0.5 for all.
  - Store normalized value in `inherited_or_estimated_fitness` for seeds.
- New molecules: `inherited_or_estimated_fitness` computed as:
  - `0.70 * max_tanimoto` to mapped seed molecules only (Morgan FP, radius=2, nBits=2048)
  - `0.20 * qed` (RDKit QED.qed)
  - `-0.10 * size_penalty` where `size_penalty = max(0, (MolWt - 500) / 200)`
  Clamp to 0.0–1.0. Mark as preliminary.

## Output Requirements

- Canonical SMILES must use RDKit MolToSmiles with `isomericSmiles=True`, `canonical=True` and `kekuleSmiles=False`.

### SDF properties (per molecule)

- `generated_id` (format: GEN_<run_id>_<counter>, where run_id defaults to the first 8 hex chars of SHA256(str(--seed)); optional `--run-id` overrides)
- `parent_1`
- `parent_2`
- `generation`
- `operation` (seed|mutation|crossover|elite)
- `inherited_or_estimated_fitness` (seed: normalized final_score; generated: estimated fitness)
- `source_final_score` (seed: real final_score; generated: propagated from best parent seed)
- `source_final_classification` (seed: real classification; generated: propagated from best parent seed)
- `canonical_smiles`
- `valid_molecule` (always true for written entries; invalid molecules are discarded)
- `notes` (semicolon-separated reasons like: missing_seed_match; mutation_failed; crossover_failed; duplicate_removed; estimated_fitness)

### CSV columns

- `generated_id`
- `canonical_smiles`
- `generation`
- `operation`
- `parent_1`
- `parent_2`
- `inherited_or_estimated_fitness`
- `source_final_score`
- `source_final_classification`
- `molecular_weight` (RDKit Descriptors.MolWt, in Da)
- `logp` (RDKit Descriptors.MolLogP)
- `hbd` (RDKit Descriptors.NumHDonors)
- `hba` (RDKit Descriptors.NumHAcceptors)
- `tpsa` (RDKit rdMolDescriptors.CalcTPSA)
- `qed` (RDKit QED.qed)
- `valid_molecule`
- `notes`

## Error Handling

Friendly errors for:
- Missing/empty input files.
- Missing required columns.
- Missing/empty SDF.
- No seeds mapped: fail with a clear error (cannot build population).
- If all seeds become invalid after sanitization/deduplication and population is empty: fail with a clear error.
- RDKit unavailable.
- Mutation/crossover failures (warn and continue).

## Reproducibility

- Use `--seed` to control Python `random` (and `numpy` if used).
- RDKit seeding requirements:
  - AllChem.EmbedMolecule(..., randomSeed=seed)
  - Any random fragment/parent selection uses Python `random` seeded by --seed.
  - BRICS fragment selection and mutation choice are driven by Python `random` only.
- RDKit steps should avoid random seeds unless explicitly set: use `AllChem.EmbedMolecule(..., randomSeed=seed)` for any 3D generation and avoid stochastic operations otherwise.
- `generated_id` is deterministic across runs with the same `--seed` unless `--run-id` is provided.

## Testing (Manual)

- Run with real `final_candidates.csv` and `compounds.sdf` as specified in the plan.
- Verify SDF and CSV outputs and required columns.
- Run twice with same `--seed` to check consistency.

## Constraints

- Do not modify existing scripts.
- No docking, ADMET local, alosteric sites, molecular dynamics.
- No git worktree usage.
## CLI Interface

Required:
- `--final-candidates` (CSV)
- `--compounds` (SDF)
- `--out` (SDF output)

Optional:
- `--out-csv` (CSV summary; if omitted, still write SDF, and skip CSV output with a warning)
- `--top-n-seeds` (default 10)
- `--min-final-score` (default 0.0)
- `--generations` (default 10)
- `--population-size` (default 30)
- `--mutation-rate` (default 0.25)
- `--crossover-rate` (default 0.50)
- `--elite-size` (default 5)
- `--seed` (default 42)
- `--run-id` (optional deterministic override for generated_id)
- `--verbose`

Outputs overwrite existing files by default; the CLI creates missing output directories. If output path is read-only, fail with a clear error (both --out and --out-csv). If --out-csv is omitted, skip CSV output and still print verbose summary.

## Selection & Elite Preservation

- Parent selection uses roulette-wheel sampling weighted by normalized fitness (`inherited_or_estimated_fitness`) for both seeds and generated molecules.
- If all fitness values are zero, fall back to uniform random selection.
- Elite handling: elites are carried forward unchanged each generation and are excluded from mutation/crossover for that generation.
- Elite size is fixed by `--elite-size`; if elite_size > population_size, clamp to population_size.
- Tie-breaks: higher fitness, then seed > generated, then shorter SMILES length.
- Mutation/crossover rates are applied per offspring by proportional weights:
  - Draw a random number and choose operation with probabilities:
    - crossover = crossover_rate
    - mutation = mutation_rate
    - clone = max(0.0, 1 - crossover_rate - mutation_rate)
  - If rates sum to >1.0, normalize to sum 1.0.
- For ancestry tracking, each molecule carries `source_final_score` and `source_final_classification` from the best seed in its lineage (propagate from parents; if both parents have values, take the higher `source_final_score`).
  - If parent values are missing, fall back to the nearest available seed in lineage; if none, leave blank and add a `notes` entry.
  - Lineage tracking uses `source_seed_id` stored for each molecule; for seeds, it equals their compound_id. For generated molecules, propagate the parent with higher `source_final_score`.
- Missing/invalid `final_score`: warn and skip those rows. Missing/invalid `final_classification` is treated as low priority.
   - If duplicate `compound_id` rows exist, keep the entry with highest normalized `final_score` (tie-break: higher priority class).
