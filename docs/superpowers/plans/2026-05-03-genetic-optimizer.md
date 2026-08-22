# Genetic Optimizer Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the new genetic optimizer that uses `final_score` seeds to generate preliminary candidates, producing SDF/CSV outputs with required metadata.

**Architecture:** Add a new `genetic_optimizer.py` module for core logic and a thin CLI in `scripts/genetic_optimize.py`. Preserve existing scripts and pipeline behavior. Outputs are deterministic with `--seed` and explicitly marked preliminary.

**Tech Stack:** Python, RDKit (Chem, AllChem, BRICS, QED, Descriptors, rdMolDescriptors), CSV, pathlib.

---

## Chunk 1: Core Module Implementation

### Task 1: Create module skeleton and data loading

**Files:**
- Create: `src/betalactamase_engine/generation/__init__.py`
- Create: `src/betalactamase_engine/generation/genetic_optimizer.py`

- [ ] **Step 1: Add package init**

```python
"""Genetic optimization utilities."""

from .genetic_optimizer import run_genetic_optimization

__all__ = ["run_genetic_optimization"]
```

- [ ] **Step 2: Add core imports and constants**

```python
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import csv
import random
import hashlib

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, QED
from rdkit.Chem import BRICS
```

- [ ] **Step 3: Implement final_candidates loader with required column validation**

```python
REQUIRED_FINAL_COLUMNS = {"compound_id", "final_score", "final_classification"}

def load_final_candidates(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    if not path.exists():
        raise FileNotFoundError(f"final_candidates.csv not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError("final_candidates.csv is empty")
    with open(path, "r") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("final_candidates.csv has no header")
        missing = REQUIRED_FINAL_COLUMNS.difference(set(reader.fieldnames))
        if missing:
            raise ValueError(f"final_candidates.csv missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    return rows, []
```

- [ ] **Step 4: Implement seed selection (normalized scores, class priority)**

```python
CLASS_PRIORITY = {
    "high_priority_candidate": 0,
    "medium_priority_candidate": 1,
    "low_priority_candidate": 2,
    "deprioritized_candidate": 3,
}

def normalize_scores(values: List[float]) -> List[float]:
    if not values:
        return []
    min_v, max_v = min(values), max(values)
    if min_v == max_v:
        return [0.5 for _ in values]
    return [(v - min_v) / (max_v - min_v) for v in values]

def select_seeds(rows: List[Dict[str, str]], top_n: int, min_score: float) -> Tuple[List[Dict[str, str]], Dict[str, float], List[str]]:
    warnings = []
    parsed = []
    for row in rows:
        raw = row.get("final_score", "")
        try:
            score = float(raw)
        except ValueError:
            warnings.append(f"Invalid final_score skipped: {raw}")
            continue
        parsed.append((row, score))
    if not parsed:
        raise ValueError("No valid final_score values found")
    scores = [s for _, s in parsed]
    use_norm = any(s < 0 or s > 1 for s in scores) or len(set(scores)) == 1
    norm_scores = normalize_scores(scores) if use_norm else scores
    if use_norm and all(v == 0.5 for v in norm_scores):
        filtered = list(zip(parsed, norm_scores))
    else:
        filtered = [(p, n) for p, n in zip(parsed, norm_scores) if n >= min_score]
    if not filtered:
        raise ValueError("No candidates left after min-final-score filtering")
    dedup = {}
    for (row, _), norm in filtered:
        cid = row.get("compound_id", "").strip()
        if not cid:
            continue
        cls = row.get("final_classification") or "low_priority_candidate"
        cls = cls if cls in CLASS_PRIORITY else "low_priority_candidate"
        prev = dedup.get(cid)
        if not prev or norm > prev[1] or (norm == prev[1] and CLASS_PRIORITY[cls] < CLASS_PRIORITY[prev[2]]):
            dedup[cid] = (row, norm, cls)
    sorted_rows = sorted(
        dedup.values(),
        key=lambda x: (CLASS_PRIORITY[x[2]], -x[1], x[0].get("compound_id", ""))
    )
    seeds = []
    norms = {}
    for row, norm, cls in sorted_rows:
        if cls == "deprioritized_candidate" and len(seeds) < top_n:
            seeds.append(row)
            norms[row.get("compound_id", "").strip()] = norm
        elif cls != "deprioritized_candidate":
            seeds.append(row)
            norms[row.get("compound_id", "").strip()] = norm
        if len(seeds) >= top_n:
            break
    return seeds, norms, warnings
```

- [ ] **Step 5: Implement SDF mapping and warnings**

```python
SDF_KEYS = ["compound_id", "chembl_id", "molecule_id", "ChEMBL ID", "ID", "name", "_Name"]

def load_sdf_library(path: Path) -> List[Chem.Mol]:
    if not path.exists():
        raise FileNotFoundError(f"compounds.sdf not found: {path}")
    suppl = Chem.SDMolSupplier(str(path), removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        raise ValueError("compounds.sdf has no valid molecules")
    return mols

def map_seeds_to_mols(seed_rows: List[Dict[str, str]], mols: List[Chem.Mol]) -> Tuple[List[Dict], List[str]]:
    warnings = []
    index = {key: {} for key in SDF_KEYS}
    for mol in mols:
        for key in SDF_KEYS:
            if mol.HasProp(key):
                value = mol.GetProp(key).strip()
                index[key].setdefault(value, []).append(mol)
    mapped = []
    for row in seed_rows:
        cid = row.get("compound_id", "").strip()
        if not cid:
            continue
        match = None
        match_count = 0
        for key in SDF_KEYS:
            hits = index.get(key, {}).get(cid, [])
            if hits:
                match = hits[0]
                match_count = len(hits)
                break
        if match is None:
            warnings.append(f"Seed not found in SDF: {cid}")
            continue
        if match_count > 1:
            warnings.append(f"Duplicate SDF entries for seed: {cid}")
        mapped.append({"row": row, "mol": match})
    if not mapped:
        raise ValueError("No seeds mapped to SDF entries")
    return mapped, warnings
```

### Task 2: Genetic operators and fitness

**Files:**
- Modify: `src/betalactamase_engine/generation/genetic_optimizer.py`

- [ ] **Step 1: Implement canonical SMILES & dedup**

```python
def canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True, kekuleSmiles=False)

def append_note(entry: Dict, note: str) -> None:
    existing = entry.get("notes", "")
    parts = [p for p in existing.split(";") if p]
    if note not in parts:
        parts.append(note)
    entry["notes"] = ";".join(parts)
```

- [ ] **Step 2: Implement mutation and crossover**

```python
TRANSFORMS = [
    "[c:1][F]>>[c:1][Cl]",
    "[c:1][Cl]>>[c:1][F]",
    "[c:1][Br]>>[c:1][Cl]",
    "[c:1][CH3]>>[c:1][OH]",
    "[c:1][OH]>>[c:1][NH2]",
    "[c:1][NH2]>>[c:1][CH3]",
    "[C:1](=O)[OH]>>[C:1](=O)[NH2]",
]

def mutate(mol: Chem.Mol, rng: random.Random) -> Optional[Chem.Mol]:
    for _ in range(8):
        smarts = rng.choice(TRANSFORMS)
        rxn = AllChem.ReactionFromSmarts(smarts)
        products = rxn.RunReactants((mol,))
        if not products:
            continue
        candidate = products[0][0]
        try:
            Chem.SanitizeMol(candidate)
            return candidate
        except Exception:
            continue
    return None

def crossover(mol_a: Chem.Mol, mol_b: Chem.Mol, rng: random.Random, fit_a: float, fit_b: float) -> Optional[Chem.Mol]:
    for _ in range(8):
        try:
            frags_a = list(BRICS.BRICSDecompose(mol_a))
            frags_b = list(BRICS.BRICSDecompose(mol_b))
            if not frags_a or not frags_b:
                return mol_a if fit_a >= fit_b else mol_b
            frag = rng.choice(frags_b)
            combined = frags_a[: max(1, len(frags_a) // 2)] + [frag]
            frag_mols = [Chem.MolFromSmiles(f) for f in combined if Chem.MolFromSmiles(f)]
            if not frag_mols:
                return mol_a if fit_a >= fit_b else mol_b
            built = list(BRICS.BRICSBuild(frag_mols))
            if not built:
                return mol_a if fit_a >= fit_b else mol_b
            candidate = rng.choice(built)
            Chem.SanitizeMol(candidate)
            return candidate
        except Exception:
            return mol_a if fit_a >= fit_b else mol_b
    return None
```

- [ ] **Step 3: Implement fitness**

```python
def seed_fitness(scores: List[float]) -> List[float]:
    use_norm = any(s < 0 or s > 1 for s in scores) or len(set(scores)) == 1
    return normalize_scores(scores) if use_norm else scores

def estimated_fitness(mol: Chem.Mol, seed_fps: List) -> float:
    if not seed_fps:
        return 0.0
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    from rdkit import DataStructs
    max_sim = max(DataStructs.BulkTanimotoSimilarity(fp, seed_fps))
    qed = QED.qed(mol)
    size_penalty = max(0.0, (Descriptors.MolWt(mol) - 500.0) / 200.0)
    score = 0.70 * max_sim + 0.20 * qed - 0.10 * size_penalty
    return max(0.0, min(score, 1.0))
```

### Task 3: GA loop, outputs, and API

**Files:**
- Modify: `src/betalactamase_engine/generation/genetic_optimizer.py`

- [ ] **Step 1: Normalize mapped seeds into internal records**

```python
def normalize_mapped_seeds(mapped: List[Dict], seed_norms: Dict[str, float], seed_order: List[str]) -> List[Dict]:
    normalized = []
    for entry in mapped:
        row = entry["row"]
        cid = row.get("compound_id", "").strip()
        norm = seed_norms.get(cid, 0.0)
        normalized.append({
            "compound_id": cid,
            "mol": entry["mol"],
            "operation": "seed",
            "generation": 0,
            "parent_1": "",
            "parent_2": "",
            "source_seed_id": cid,
            "source_final_score": row.get("final_score", ""),
            "source_final_classification": row.get("final_classification", ""),
            "inherited_or_estimated_fitness": norm,
            "valid_molecule": True,
            "notes": "",
        })
    # preserve original seed selection order
    normalized.sort(key=lambda s: seed_order.index(s["compound_id"]) if s["compound_id"] in seed_order else len(seed_order))
    return normalized
```

- [ ] **Step 2: Implement initial population builder**

```python
def build_initial_population(seed_records: List[Dict], population_size: int, rng: random.Random, seed_fps: List) -> Tuple[List[Dict], List[str]]:
    # seeds sorted by selection order
    population = []
    warnings = []
    # add seeds first
    for seed in seed_records:
        population.append(seed)
        if len(population) >= population_size:
            return population, warnings
    # round-robin variants (max 3 per seed)
    variant_counts = {s["compound_id"]: 0 for s in seed_records}
    while len(population) < population_size:
        progressed = False
        for seed in seed_records:
            if len(population) >= population_size:
                break
            cid = seed["compound_id"]
            if variant_counts[cid] >= 3:
                continue
            variant = mutate(seed["mol"], rng)
            notes = "estimated_fitness"
            if variant is None:
                variant = seed["mol"]
                notes = "estimated_fitness;mutation_failed"
            est_fit = estimated_fitness(variant, seed_fps)
            population.append({
                "mol": variant,
                "operation": "mutation",
                "generation": 0,
                "parent_1": cid,
                "parent_2": "",
                "source_seed_id": cid,
                "source_final_score": seed["source_final_score"],
                "source_final_classification": seed["source_final_classification"],
                "inherited_or_estimated_fitness": est_fit,
                "valid_molecule": True,
                "notes": notes,
            })
            variant_counts[cid] += 1
            progressed = True
        if not progressed:
            # allow duplicates with warning at caller
            warnings.append("Population filled with duplicate seeds")
            population.append(seed_records[0])
    return population, warnings
```

- [ ] **Step 3: Implement GA loop and lineage propagation**

```python
def run_genetic_optimization(...):
    # build initial population with source_seed_id/source_final_score/source_final_classification
    # per-generation:
    # - roulette parent selection on normalized fitness (uniform fallback if all zero)
    # - normalize rates if mutation+crossover > 1.0
    # - clamp elite_size to population_size and mark elite operation
    # - elite carry-forward (excluded from mutation/crossover)
    # - generate offspring by rates; mutate/crossover/clone
    #   - for any generated molecule, notes += "estimated_fitness"
    #   - if mutation fails, notes += "mutation_failed"
    #   - if crossover falls back, notes += "crossover_failed"
    # - propagate source_seed_id/source_final_score/source_final_classification from best parent (higher source_final_score)
    # - if parent values missing, fallback to source_seed_id lineage; if none, notes += "missing_seed_match"
    # - discard invalid molecules (notes += "invalid_molecule")
    # - per-generation dedup by canonical SMILES; keep highest fitness, tie-break seed>generated then shorter smiles
    # - refill with best remaining non-elite pool or clone elite if needed (warn on refill, notes += "duplicate_removed")
    # final dedup before output
    return population
```

- [ ] **Step 4: Implement SDF/CSV writers with required properties**

```python
def write_sdf(mols, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_path))
    for entry in mols:
        mol = entry["mol"]
        for key, value in entry["props"].items():
            mol.SetProp(key, str(value))
        writer.write(mol)
    writer.close()

def write_csv(rows, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "generated_id",
        "canonical_smiles",
        "generation",
        "operation",
        "parent_1",
        "parent_2",
        "inherited_or_estimated_fitness",
        "source_final_score",
        "source_final_classification",
        "molecular_weight",
        "logp",
        "hbd",
        "hba",
        "tpsa",
        "qed",
        "valid_molecule",
        "notes",
    ]
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

- [ ] **Step 5: Add generated_id/run_id and notes semantics**

```python
def build_run_id(seed: int, run_id: Optional[str]) -> str:
    if run_id:
        return run_id
    digest = hashlib.sha256(str(seed).encode("utf-8")).hexdigest()
    return digest[:8]

def build_generated_id(run_id: str, counter: int) -> str:
    return f"GEN_{run_id}_{counter:05d}"

def build_props(entry, counter, run_id):
    return {
        "generated_id": build_generated_id(run_id, counter),
        "canonical_smiles": canonical_smiles(entry["mol"]),
        "generation": entry.get("generation", 0),
        "operation": entry.get("operation", ""),
        "parent_1": entry.get("parent_1", ""),
        "parent_2": entry.get("parent_2", ""),
        "inherited_or_estimated_fitness": entry.get("inherited_or_estimated_fitness", ""),
        "source_final_score": entry.get("source_final_score", ""),
        "source_final_classification": entry.get("source_final_classification", ""),
        "valid_molecule": True,
        "notes": entry.get("notes", ""),
    }

def build_csv_row(entry):
    mol = entry["mol"]
    return {
        "generated_id": entry["props"]["generated_id"],
        "canonical_smiles": entry["props"]["canonical_smiles"],
        "generation": entry["props"]["generation"],
        "operation": entry["props"]["operation"],
        "parent_1": entry["props"]["parent_1"],
        "parent_2": entry["props"]["parent_2"],
        "inherited_or_estimated_fitness": entry["props"]["inherited_or_estimated_fitness"],
        "source_final_score": entry["props"]["source_final_score"],
        "source_final_classification": entry["props"]["source_final_classification"],
        "molecular_weight": f"{Descriptors.MolWt(mol):.2f}",
        "logp": f"{Descriptors.MolLogP(mol):.2f}",
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "qed": f"{QED.qed(mol):.3f}",
        "valid_molecule": True,
        "notes": entry["props"].get("notes", ""),
    }
```

- [ ] **Step 5: Add error handling for RDKit availability and empty populations**

```python
try:
    from rdkit import Chem
except ImportError as exc:
    raise RuntimeError("RDKit not available") from exc

if not population:
    raise ValueError("Population empty after sanitization/dedup")
```
```

- [ ] **Step 4: Expose top-level API**

```python
def run_genetic_optimization(...):
    # orchestrate load → select → map → evolve → write
    return summary
```

---

## Chunk 2: CLI + Documentation + Manual Tests

### Task 4: CLI script

**Files:**
- Create: `scripts/genetic_optimize.py`

- [ ] **Step 1: Create argparse CLI matching spec**

```python
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
parser.add_argument("--verbose", action="store_true")
```

- [ ] **Step 1b: Validate input paths before running**

```python
if not Path(args.final_candidates).exists():
    raise SystemExit("final_candidates.csv not found")
if not Path(args.compounds).exists():
    raise SystemExit("compounds.sdf not found")
```

- [ ] **Step 2: Wire CLI to module and print summary**

```python
if args.verbose:
    print("Genetic optimization completed.")
    print(f"Final candidates input: {args.final_candidates}")
    print(f"Seed molecules selected: {summary['seed_selected']}")
    print(f"Seed molecules found in SDF: {summary['seed_mapped']}")
    print(f"Generations: {args.generations}")
    print(f"Population size: {args.population_size}")
    print(f"Generated valid molecules: {summary['generated_valid']}")
    print(f"Unique generated molecules: {summary['generated_unique']}")
    print(f"Output SDF: {args.out}")
    if args.out_csv:
        print(f"Output CSV: {args.out_csv}")
    else:
        print("Output CSV: not provided")
    if summary.get('warnings'):
        print("Warnings:")
        for warning in summary['warnings']:
            print(f" - {warning}")
print("Generated candidates should be submitted to docking and consensus scoring in the next validation cycle.")
```

- [ ] **Step 3: Add CLI output handling rules**

```python
# create output directories (overwrite by default)
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
if args.out_csv:
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
# fail clearly if output paths are read-only
# if --out-csv omitted, warn and skip CSV (always warn, even without --verbose)
```

### Task 5: Documentation

**Files:**
- Create: `docs/genetic_optimizer.md`
- Modify: `README.md`

- [ ] **Step 1: Write docs with usage, limitations, and output interpretation**
- [ ] **Step 2: Update README pipeline step**

### Task 6: Manual tests

**Files:**
- Test: `scripts/genetic_optimize.py`

- [ ] **Step 1: Run minimal test (per plan)**

Run:
```bash
python scripts/genetic_optimize.py \
  --final-candidates results/experiments/teste/final_candidates.csv \
  --compounds data/compounds/compounds.sdf \
  --out results/experiments/teste/generated_candidates.sdf \
  --out-csv results/experiments/teste/generated_candidates_summary.csv \
  --generations 3 \
  --population-size 10 \
  --verbose
```

Expected:
- `generated_candidates.sdf` created
- `generated_candidates_summary.csv` created
- Summary printed with counts

- [ ] **Step 1b: Verify output columns/properties**

Check CSV has: `generated_id`, `canonical_smiles`, `generation`, `operation`, `parent_1`, `parent_2`, `inherited_or_estimated_fitness`, `source_final_score`, `source_final_classification`, `molecular_weight`, `logp`, `hbd`, `hba`, `tpsa`, `qed`, `valid_molecule`, `notes`.

Check SDF properties include: `generated_id`, `parent_1`, `parent_2`, `generation`, `operation`, `inherited_or_estimated_fitness`, `source_final_score`, `source_final_classification`, `canonical_smiles`, `valid_molecule`, `notes`.

- [ ] **Step 1c: Run without --out-csv and confirm warning**

Run the same command without `--out-csv` and confirm a warning is printed and SDF still written.

- [ ] **Step 2: Reproducibility check**

Run the same command twice with `--seed 42` and confirm the output counts are identical.

---

Plan complete and saved to `docs/superpowers/plans/2026-05-03-genetic-optimizer.md`. Ready to execute?
