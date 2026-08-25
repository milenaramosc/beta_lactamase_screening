# Unified Genetic Pipeline Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor the molecular generation scripts so all AG modes use one shared engine and generated candidates can return directly to docking.

**Architecture:** Extend `src/betalactamase_engine/generation/genetic_optimizer.py` into the common AG engine, with shared seed records, configurable filters, and loaders for ranking, consensus, and protein modes. Keep the three existing scripts as compatibility CLIs that delegate to this engine. Add `--ligands` to `scripts/screen_and_rank.py`.

**Tech Stack:** Python 3, RDKit, Bio.PDB, AutoDock Vina, Meeko/OpenBabel, CSV/SDF outputs.

---

## Chunk 1: Shared AG Engine

### Task 1: Add common seed and mode APIs

**Files:**
- Modify: `src/betalactamase_engine/generation/genetic_optimizer.py`

- [x] Add a `SeedRecord` dataclass with `compound_id`, `mol`, `source_score`, `source_classification`, and `source_mode`.
- [x] Add RO5-aware `ChemicalFilters` fields for `max_logp`.
- [x] Add `load_ranking_seeds`, `load_consensus_seeds`, and `build_known_inhibitor_seeds`.
- [x] Add a shared `run_genetic_optimization_from_seeds` function used by all modes.
- [x] Preserve the existing `run_genetic_optimization` wrapper for consensus mode.

## Chunk 2: Script Wrappers

### Task 2: Convert old scripts into wrappers

**Files:**
- Modify: `scripts/generate_candidate.py`
- Modify: `scripts/generate_inhibitors_from_protein.py`
- Modify: `scripts/genetic_optimize.py`

- [x] Make `generate_candidate.py` call the shared engine with ranking seeds.
- [x] Make `genetic_optimize.py` call the shared engine with consensus seeds.
- [x] Make `generate_inhibitors_from_protein.py` call the shared engine with protein-based fitness.
- [x] Keep legacy argument names where practical.

## Chunk 3: Docking Reintegration

### Task 3: Allow generated SDF docking

**Files:**
- Modify: `scripts/screen_and_rank.py`
- Modify: `README.md`
- Modify: `docs/genetic_optimizer.md`

- [x] Add `--ligands` to `screen_and_rank.py`.
- [x] Use the provided SDF instead of `data/compounds/compounds.sdf` when supplied.
- [x] Derive ligand names from `generated_id`, `_Name`, `compound_id`, or fallback index.
- [x] Document the AG -> docking -> interaction -> consensus cycle.

## Chunk 4: Validation

### Task 4: Verify behavior

**Files:**
- No persistent project outputs.

- [x] Run `python -m compileall scripts src`.
- [x] Run help for all modified CLIs.
- [x] Run a short consensus AG execution to `/tmp`.
- [x] Run a short ranking AG execution to `/tmp`.
- [x] Confirm Git status and avoid modifying unrelated user changes.
