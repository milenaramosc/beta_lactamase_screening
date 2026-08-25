from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple
import csv
import hashlib
import random

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, BRICS, Descriptors, QED, rdMolDescriptors
    from rdkit import DataStructs
except ImportError as exc:
    raise RuntimeError("RDKit not available") from exc


REQUIRED_FINAL_COLUMNS = {"compound_id", "final_score", "final_classification"}

CLASS_PRIORITY = {
    "high_priority_candidate": 0,
    "medium_priority_candidate": 1,
    "low_priority_candidate": 2,
    "deprioritized_candidate": 3,
}

SDF_KEYS = [
    "compound_id",
    "chembl_id",
    "molecule_id",
    "ChEMBL ID",
    "ID",
    "name",
    "_Name",
]

TRANSFORMS = [
    "[c:1][F]>>[c:1][Cl]",
    "[c:1][Cl]>>[c:1][F]",
    "[c:1][Br]>>[c:1][Cl]",
    "[c:1][CH3]>>[c:1][OH]",
    "[c:1][OH]>>[c:1][NH2]",
    "[c:1][NH2]>>[c:1][CH3]",
    "[C:1](=O)[OH]>>[C:1](=O)[NH2]",
]

FitnessFunction = Callable[[Chem.Mol, List], float]


@dataclass
class OptimizationSummary:
    seed_selected: int
    seed_mapped: int
    generated_valid: int
    generated_unique: int
    filtered: int
    filtered_csv: Optional[Path]
    warnings: List[str]


@dataclass(frozen=True)
class SeedRecord:
    compound_id: str
    mol: Chem.Mol
    source_score: float
    source_classification: str = ""
    source_mode: str = ""


@dataclass(frozen=True)
class ChemicalFilters:
    max_molecular_weight: float = 500.0
    max_logp: float = 5.0
    max_tpsa: float = 250.0
    max_hbd: int = 5
    max_hba: int = 10
    min_qed: float = 0.05


FILTERED_FIELDNAMES = [
    "generated_id",
    "canonical_smiles",
    "generation",
    "operation",
    "parent_1",
    "parent_2",
    "filter_reason",
    "molecular_weight",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "qed",
    "notes",
    "valid_molecule",
    "filtered",
    "source_final_score",
    "source_final_classification",
    "source_mode",
]


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
            raise ValueError(
                "final_candidates.csv missing columns: " + ", ".join(sorted(missing))
            )
        rows = list(reader)
    if not rows:
        raise ValueError("final_candidates.csv is empty")
    return rows, []


def normalize_scores(values: List[float]) -> List[float]:
    if not values:
        return []
    min_v, max_v = min(values), max(values)
    if min_v == max_v:
        return [0.5 for _ in values]
    return [(v - min_v) / (max_v - min_v) for v in values]


def select_seeds(
    rows: List[Dict[str, str]],
    top_n: int,
    min_score: float,
) -> Tuple[List[Dict[str, str]], Dict[str, float], List[str]]:
    warnings: List[str] = []
    parsed: List[Tuple[Dict[str, str], float]] = []
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
    dedup: Dict[str, Tuple[Dict[str, str], float, str]] = {}
    for (row, _), norm in filtered:
        cid = row.get("compound_id", "").strip()
        if not cid:
            continue
        cls = row.get("final_classification") or "low_priority_candidate"
        cls = cls if cls in CLASS_PRIORITY else "low_priority_candidate"
        prev = dedup.get(cid)
        if not prev or norm > prev[1] or (
            norm == prev[1] and CLASS_PRIORITY[cls] < CLASS_PRIORITY[prev[2]]
        ):
            dedup[cid] = (row, norm, cls)
    sorted_rows = sorted(
        dedup.values(),
        key=lambda x: (CLASS_PRIORITY[x[2]], -x[1], x[0].get("compound_id", "")),
    )
    seeds: List[Dict[str, str]] = []
    norms: Dict[str, float] = {}
    for row, norm, cls in sorted_rows:
        if cls == "deprioritized_candidate" and len(seeds) < top_n:
            seeds.append(row)
            norms[row.get("compound_id", "").strip()] = norm
        elif cls != "deprioritized_candidate":
            seeds.append(row)
            norms[row.get("compound_id", "").strip()] = norm
        if len(seeds) >= top_n:
            break
    if not seeds:
        raise ValueError("No candidates available for seeding")
    return seeds, norms, warnings


def load_sdf_library(path: Path) -> List[Chem.Mol]:
    if not path.exists():
        raise FileNotFoundError(f"compounds.sdf not found: {path}")
    suppl = Chem.SDMolSupplier(str(path), removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        raise ValueError("compounds.sdf has no valid molecules")

    def sort_key(mol: Chem.Mol) -> Tuple[str, ...]:
        values = []
        for key in SDF_KEYS:
            values.append(mol.GetProp(key).strip() if mol.HasProp(key) else "")
        try:
            values.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
        except Exception:
            values.append("")
        return tuple(values)

    return sorted(mols, key=sort_key)


def _mol_name_index(mols: List[Chem.Mol]) -> Dict[str, List[Chem.Mol]]:
    index: Dict[str, List[Chem.Mol]] = {}
    for mol in mols:
        for key in SDF_KEYS:
            if mol.HasProp(key):
                value = mol.GetProp(key).strip()
                if value:
                    index.setdefault(value, []).append(mol)
    return index


def _lookup_mol(index: Dict[str, List[Chem.Mol]], compound_id: str) -> Optional[Chem.Mol]:
    hits = index.get(compound_id.strip())
    return hits[0] if hits else None


def load_ranking_seeds(
    ranking_path: Path,
    compounds_path: Path,
    top_n: int = 20,
) -> Tuple[List[SeedRecord], List[str]]:
    if not ranking_path.exists():
        raise FileNotFoundError(f"ranking.csv not found: {ranking_path}")
    if ranking_path.stat().st_size == 0:
        raise ValueError("ranking.csv is empty")

    with open(ranking_path, "r") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("ranking.csv has no header")
        required = {"compound_id", "dg_kcal_mol"}
        missing = required.difference(set(reader.fieldnames))
        if missing:
            raise ValueError("ranking.csv missing columns: " + ", ".join(sorted(missing)))
        parsed = []
        for row in reader:
            compound_id = (row.get("compound_id") or "").strip()
            if not compound_id:
                continue
            try:
                dg = float(row.get("dg_kcal_mol", ""))
            except ValueError:
                continue
            parsed.append((compound_id, dg))

    if not parsed:
        raise ValueError("No valid docking scores found in ranking.csv")

    parsed = sorted(parsed, key=lambda item: item[1])[:top_n]
    values = [dg for _, dg in parsed]
    best, worst = min(values), max(values)
    if best == worst:
        scores = [0.5 for _ in values]
    else:
        scores = [(worst - value) / (worst - best) for value in values]

    mols = load_sdf_library(compounds_path)
    index = _mol_name_index(mols)
    seeds: List[SeedRecord] = []
    warnings: List[str] = []
    for (compound_id, _), score in zip(parsed, scores):
        mol = _lookup_mol(index, compound_id)
        if mol is None:
            warnings.append(f"Seed not found in SDF: {compound_id}")
            continue
        seeds.append(
            SeedRecord(
                compound_id=compound_id,
                mol=mol,
                source_score=score,
                source_classification="docking_seed",
                source_mode="ranking",
            )
        )

    if not seeds:
        raise ValueError("No ranking seeds mapped to SDF entries")
    return seeds, warnings


def load_consensus_seeds(
    final_candidates_path: Path,
    compounds_path: Path,
    top_n: int = 10,
    min_final_score: float = 0.0,
) -> Tuple[List[SeedRecord], List[str]]:
    rows, warnings = load_final_candidates(final_candidates_path)
    seed_rows, seed_norms, select_warnings = select_seeds(
        rows, top_n=top_n, min_score=min_final_score
    )
    warnings.extend(select_warnings)
    mols = load_sdf_library(compounds_path)
    mapped, map_warnings = map_seeds_to_mols(seed_rows, mols)
    warnings.extend(map_warnings)

    seeds: List[SeedRecord] = []
    for entry in mapped:
        row = entry["row"]
        compound_id = str(row.get("compound_id", "")).strip()
        seeds.append(
            SeedRecord(
                compound_id=compound_id,
                mol=entry["mol"],
                source_score=seed_norms.get(compound_id, 0.0),
                source_classification=row.get("final_classification", ""),
                source_mode="consensus",
            )
        )
    if not seeds:
        raise ValueError("No consensus seeds mapped to SDF entries")
    return seeds, warnings


def build_known_inhibitor_seeds(inhibitors: Dict[str, str]) -> Tuple[List[SeedRecord], List[str]]:
    seeds: List[SeedRecord] = []
    warnings: List[str] = []
    for name, smiles in inhibitors.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            warnings.append(f"Known inhibitor SMILES invalid: {name}")
            continue
        seeds.append(
            SeedRecord(
                compound_id=name,
                mol=mol,
                source_score=1.0,
                source_classification="known_inhibitor",
                source_mode="protein",
            )
        )
    if not seeds:
        raise ValueError("No valid known inhibitor seeds")
    return seeds, warnings


def map_seeds_to_mols(
    seed_rows: List[Dict[str, str]],
    mols: List[Chem.Mol],
) -> Tuple[List[Dict[str, object]], List[str]]:
    warnings: List[str] = []
    index: Dict[str, Dict[str, List[Chem.Mol]]] = {key: {} for key in SDF_KEYS}
    for mol in mols:
        for key in SDF_KEYS:
            if mol.HasProp(key):
                value = mol.GetProp(key).strip()
                index[key].setdefault(value, []).append(mol)
    mapped: List[Dict[str, object]] = []
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


def canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(
        mol,
        isomericSmiles=True,
        canonical=True,
        kekuleSmiles=False,
    )


def append_note(entry: Dict[str, object], note: str) -> None:
    existing = str(entry.get("notes", ""))
    parts = [p for p in existing.split(";") if p]
    if note not in parts:
        parts.append(note)
    entry["notes"] = ";".join(parts)


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


def crossover(
    mol_a: Chem.Mol,
    mol_b: Chem.Mol,
    rng: random.Random,
    fit_a: float,
    fit_b: float,
) -> Optional[Chem.Mol]:
    for _ in range(8):
        try:
            frags_a = sorted(BRICS.BRICSDecompose(mol_a))
            frags_b = sorted(BRICS.BRICSDecompose(mol_b))
            if not frags_a or not frags_b:
                return mol_a if fit_a >= fit_b else mol_b
            frag = rng.choice(frags_b)
            combined = frags_a[: max(1, len(frags_a) // 2)] + [frag]
            frag_mols = [Chem.MolFromSmiles(f) for f in combined if Chem.MolFromSmiles(f)]
            if not frag_mols:
                return mol_a if fit_a >= fit_b else mol_b
            built = sorted(
                BRICS.BRICSBuild(frag_mols),
                key=lambda candidate: Chem.MolToSmiles(
                    candidate, canonical=True, isomericSmiles=True
                ),
            )
            if not built:
                return mol_a if fit_a >= fit_b else mol_b
            candidate = rng.choice(built)
            Chem.SanitizeMol(candidate)
            return candidate
        except Exception:
            return mol_a if fit_a >= fit_b else mol_b
    return None


def seed_fitness(scores: List[float]) -> List[float]:
    use_norm = any(s < 0 or s > 1 for s in scores) or len(set(scores)) == 1
    return normalize_scores(scores) if use_norm else scores


def estimated_fitness(mol: Chem.Mol, seed_fps: List) -> float:
    if not seed_fps:
        return 0.0
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    max_sim = max(DataStructs.BulkTanimotoSimilarity(fp, seed_fps))
    qed = QED.qed(mol)
    mw = Descriptors.MolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    penalty_mw = max(0.0, (mw - 500.0) / 200.0)
    penalty_tpsa = max(0.0, (tpsa - 180.0) / 80.0)
    smiles = canonical_smiles(mol)
    fragments = len(Chem.GetMolFrags(mol, asMols=False))
    score = 0.60 * max_sim + 0.20 * qed - 0.12 * penalty_mw - 0.06 * penalty_tpsa
    if "." in smiles:
        score -= 0.10
    if fragments > 3:
        score -= 0.05
    return max(0.0, min(score, 1.0))


def _score_molecule(
    mol: Chem.Mol,
    seed_fps: List,
    fitness_fn: Optional[FitnessFunction] = None,
) -> float:
    if fitness_fn is not None:
        return max(0.0, min(float(fitness_fn(mol, seed_fps)), 1.0))
    return estimated_fitness(mol, seed_fps)


def normalize_mapped_seeds(
    mapped: List[Dict[str, object]],
    seed_norms: Dict[str, float],
    seed_order: List[str],
) -> List[Dict[str, object]]:
    normalized = []
    for entry in mapped:
        row = entry["row"]
        cid = str(row.get("compound_id", "")).strip()
        norm = seed_norms.get(cid, 0.0)
        normalized.append(
            {
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
            }
        )
    normalized.sort(
        key=lambda s: seed_order.index(s["compound_id"])
        if s["compound_id"] in seed_order
        else len(seed_order)
    )
    return normalized


def build_initial_population(
    seed_records: List[Dict[str, object]],
    population_size: int,
    rng: random.Random,
    seed_fps: List,
    fitness_fn: Optional[FitnessFunction] = None,
) -> Tuple[List[Dict[str, object]], List[str]]:
    population: List[Dict[str, object]] = []
    warnings: List[str] = []
    for seed in seed_records:
        population.append(seed)
        if len(population) >= population_size:
            return population, warnings
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
            est_fit = _score_molecule(variant, seed_fps, fitness_fn)
            population.append(
                {
                    "compound_id": cid,
                    "mol": variant,
                    "operation": "mutation",
                    "generation": 0,
                    "parent_1": cid,
                    "parent_2": "",
                    "source_seed_id": cid,
                    "source_final_score": seed.get("source_final_score", ""),
                    "source_final_classification": seed.get(
                        "source_final_classification", ""
                    ),
                    "inherited_or_estimated_fitness": est_fit,
                    "valid_molecule": True,
                    "notes": notes,
                }
            )
            variant_counts[cid] += 1
            progressed = True
        if not progressed:
            warnings.append("Population filled with duplicate seeds")
            population.append(seed_records[0])
    return population, warnings


def _weighted_choice(entries: List[Dict[str, object]], rng: random.Random) -> Dict[str, object]:
    weights = [float(entry.get("inherited_or_estimated_fitness", 0.0)) for entry in entries]
    if not weights or max(weights) <= 0:
        return rng.choice(entries)
    return rng.choices(entries, weights=weights, k=1)[0]


def _best_parent(a: Dict[str, object], b: Dict[str, object]) -> Dict[str, object]:
    def _score(entry: Dict[str, object]) -> float:
        raw = entry.get("source_final_score")
        try:
            return float(raw) if raw not in (None, "") else -1.0
        except ValueError:
            return -1.0

    score_a = _score(a)
    score_b = _score(b)
    if score_a == score_b:
        if a.get("source_seed_id") and not b.get("source_seed_id"):
            return a
        if b.get("source_seed_id") and not a.get("source_seed_id"):
            return b
        return a
    return a if score_a > score_b else b


def _dedup_population(population: List[Dict[str, object]]) -> Tuple[List[Dict[str, object]], int]:
    seen: Dict[str, Dict[str, object]] = {}
    removed = 0
    for entry in population:
        mol = entry.get("mol")
        if mol is None:
            removed += 1
            continue
        smiles = canonical_smiles(mol)
        existing = seen.get(smiles)
        if existing is None:
            seen[smiles] = entry
            continue
        best = entry
        if float(existing.get("inherited_or_estimated_fitness", 0.0)) > float(
            entry.get("inherited_or_estimated_fitness", 0.0)
        ):
            best = existing
        if existing.get("operation") == "seed" and entry.get("operation") != "seed":
            best = existing
        elif entry.get("operation") == "seed" and existing.get("operation") != "seed":
            best = entry
        else:
            if len(canonical_smiles(existing.get("mol"))) <= len(smiles):
                best = existing
        if best is not existing:
            removed += 1
        seen[smiles] = best
    return list(seen.values()), removed


def build_run_id(seed: int, run_id: Optional[str]) -> str:
    if run_id:
        return run_id
    digest = hashlib.sha256(str(seed).encode("utf-8")).hexdigest()
    return digest[:8]


def build_generated_id(run_id: str, counter: int) -> str:
    return f"GEN_{run_id}_{counter:05d}"


def build_props(entry: Dict[str, object], counter: int, run_id: str) -> Dict[str, object]:
    mol = entry["mol"]
    return {
        "generated_id": build_generated_id(run_id, counter),
        "canonical_smiles": canonical_smiles(mol),
        "generation": entry.get("generation", 0),
        "operation": entry.get("operation", ""),
        "parent_1": entry.get("parent_1", ""),
        "parent_2": entry.get("parent_2", ""),
        "inherited_or_estimated_fitness": entry.get(
            "inherited_or_estimated_fitness", ""
        ),
        "source_final_score": entry.get("source_final_score", ""),
        "source_final_classification": entry.get("source_final_classification", ""),
        "source_mode": entry.get("source_mode", ""),
        "valid_molecule": True,
        "filter_reason": "",
        "notes": entry.get("notes", ""),
    }


def build_csv_row(entry: Dict[str, object]) -> Dict[str, object]:
    mol = entry["mol"]
    return {
        "generated_id": entry["props"]["generated_id"],
        "canonical_smiles": entry["props"]["canonical_smiles"],
        "generation": entry["props"]["generation"],
        "operation": entry["props"]["operation"],
        "parent_1": entry["props"]["parent_1"],
        "parent_2": entry["props"]["parent_2"],
        "inherited_or_estimated_fitness": entry["props"][
            "inherited_or_estimated_fitness"
        ],
        "source_final_score": entry["props"]["source_final_score"],
        "source_final_classification": entry["props"]["source_final_classification"],
        "source_mode": entry["props"].get("source_mode", ""),
        "molecular_weight": f"{Descriptors.MolWt(mol):.2f}",
        "logp": f"{Descriptors.MolLogP(mol):.2f}",
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "qed": f"{QED.qed(mol):.3f}",
        "valid_molecule": True,
        "filter_reason": "",
        "notes": entry["props"].get("notes", ""),
    }


def molecule_metrics(mol: Chem.Mol) -> Dict[str, object]:
    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "qed": QED.qed(mol),
    }


def filter_reasons(mol: Chem.Mol, filters: ChemicalFilters) -> List[str]:
    metrics = molecule_metrics(mol)
    reasons = []
    if metrics["molecular_weight"] > filters.max_molecular_weight:
        reasons.append(f"mw>{filters.max_molecular_weight:g}")
    if metrics["logp"] > filters.max_logp:
        reasons.append(f"logp>{filters.max_logp:g}")
    if metrics["tpsa"] > filters.max_tpsa:
        reasons.append(f"tpsa>{filters.max_tpsa:g}")
    if metrics["hbd"] > filters.max_hbd:
        reasons.append(f"hbd>{filters.max_hbd:g}")
    if metrics["hba"] > filters.max_hba:
        reasons.append(f"hba>{filters.max_hba:g}")
    if metrics["qed"] < filters.min_qed:
        reasons.append(f"qed<{filters.min_qed:g}")
    return reasons


def build_filtered_row(
    entry: Dict[str, object],
    reason: str,
    generated_id: str = "",
) -> Dict[str, object]:
    mol = entry.get("mol")
    metrics = {}
    smiles = ""
    if mol is not None and reason != "invalid_molecule":
        try:
            metrics = molecule_metrics(mol)
            smiles = canonical_smiles(mol)
        except Exception:
            reason = "invalid_molecule"
            metrics = {}
            smiles = ""

    return {
        "generated_id": generated_id,
        "canonical_smiles": smiles,
        "generation": entry.get("generation", ""),
        "operation": entry.get("operation", ""),
        "parent_1": entry.get("parent_1", ""),
        "parent_2": entry.get("parent_2", ""),
        "filter_reason": reason,
        "molecular_weight": f"{metrics['molecular_weight']:.2f}" if metrics else "",
        "logp": f"{metrics['logp']:.2f}" if metrics else "",
        "tpsa": f"{metrics['tpsa']:.2f}" if metrics else "",
        "hbd": metrics["hbd"] if metrics else "",
        "hba": metrics["hba"] if metrics else "",
        "qed": f"{metrics['qed']:.3f}" if metrics else "",
        "notes": entry.get("notes", ""),
        "valid_molecule": "false",
        "filtered": "true",
        "source_final_score": entry.get("source_final_score", ""),
        "source_final_classification": entry.get("source_final_classification", ""),
        "source_mode": entry.get("source_mode", ""),
    }


def write_filtered_csv(rows: List[Dict[str, object]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FILTERED_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def write_sdf(entries: List[Dict[str, object]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_path))
    for entry in entries:
        mol = entry["mol"]
        for key, value in entry["props"].items():
            mol.SetProp(key, str(value))
        mol.SetProp("filter_reason", "")
        writer.write(mol)
    writer.close()


def write_csv(rows: List[Dict[str, object]], out_path: Path) -> None:
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
        "source_mode",
        "molecular_weight",
        "logp",
        "hbd",
        "hba",
        "tpsa",
        "qed",
        "valid_molecule",
        "filter_reason",
        "notes",
    ]
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _sort_population_for_output(population: List[Dict[str, object]]) -> List[Dict[str, object]]:
    return sorted(
        population,
        key=lambda entry: (
            -float(entry.get("inherited_or_estimated_fitness", 0.0)),
            canonical_smiles(entry["mol"]),
        ),
    )


def _try_replenish_unique_population(
    population: List[Dict[str, object]],
    population_size: int,
    rng: random.Random,
    seed_fps: List,
    filters: ChemicalFilters,
    warnings: List[str],
    filtered_rows: List[Dict[str, object]],
    fitness_fn: Optional[FitnessFunction] = None,
) -> List[Dict[str, object]]:
    if population_size <= 1 or len(population) >= population_size or not population:
        return population

    seen = {canonical_smiles(entry["mol"]) for entry in population}
    attempts = 0
    max_attempts = population_size * 8

    while len(population) < population_size and attempts < max_attempts:
        attempts += 1
        parent_a = _weighted_choice(population, rng)
        parent_b = _weighted_choice(population, rng)
        operation = "mutation" if rng.random() < 0.5 else "crossover"

        if operation == "mutation":
            child_mol = mutate(parent_a["mol"], rng)
        else:
            child_mol = crossover(
                parent_a["mol"],
                parent_b["mol"],
                rng,
                float(parent_a.get("inherited_or_estimated_fitness", 0.0)),
                float(parent_b.get("inherited_or_estimated_fitness", 0.0)),
            )

        best_parent = _best_parent(parent_a, parent_b)
        child_entry: Dict[str, object] = {
            "compound_id": best_parent.get("compound_id", ""),
            "mol": child_mol,
            "operation": operation,
            "generation": "replenish",
            "parent_1": parent_a.get("source_seed_id", ""),
            "parent_2": parent_b.get("source_seed_id", "") if operation == "crossover" else "",
            "source_seed_id": best_parent.get("source_seed_id", ""),
            "source_final_score": best_parent.get("source_final_score", ""),
            "source_final_classification": best_parent.get("source_final_classification", ""),
            "inherited_or_estimated_fitness": 0.0,
            "valid_molecule": True,
            "notes": "estimated_fitness;diversity_replenishment",
        }

        if child_mol is None:
            child_entry["valid_molecule"] = False
            filtered_rows.append(build_filtered_row(child_entry, "invalid_molecule"))
            continue

        try:
            Chem.SanitizeMol(child_mol)
            smiles = canonical_smiles(child_mol)
        except Exception:
            child_entry["valid_molecule"] = False
            filtered_rows.append(build_filtered_row(child_entry, "invalid_molecule"))
            continue

        if smiles in seen:
            continue

        reasons = filter_reasons(child_mol, filters)
        if reasons:
            filtered_rows.append(build_filtered_row(child_entry, ";".join(reasons)))
            continue

        child_entry["inherited_or_estimated_fitness"] = _score_molecule(
            child_mol, seed_fps, fitness_fn
        )
        population.append(child_entry)
        seen.add(smiles)

    if len(population) < population_size:
        warnings.append(
            f"Only {len(population)} unique molecules could be generated after deduplication."
        )
    return population


def _seed_records_to_population(seeds: List[SeedRecord]) -> List[Dict[str, object]]:
    normalized = []
    for seed in seeds:
        normalized.append(
            {
                "compound_id": seed.compound_id,
                "mol": seed.mol,
                "operation": "seed",
                "generation": 0,
                "parent_1": "",
                "parent_2": "",
                "source_seed_id": seed.compound_id,
                "source_final_score": f"{seed.source_score:.6g}",
                "source_final_classification": seed.source_classification,
                "source_mode": seed.source_mode,
                "inherited_or_estimated_fitness": seed.source_score,
                "valid_molecule": True,
                "notes": "",
            }
        )
    return normalized


def run_genetic_optimization_from_seeds(
    seeds: List[SeedRecord],
    out_sdf: Path,
    out_csv: Optional[Path],
    generations: int = 10,
    population_size: int = 30,
    mutation_rate: float = 0.25,
    crossover_rate: float = 0.50,
    elite_size: int = 5,
    seed: int = 42,
    run_id: Optional[str] = None,
    max_molecular_weight: float = 500.0,
    max_logp: float = 5.0,
    max_tpsa: float = 250.0,
    max_hbd: int = 5,
    max_hba: int = 10,
    min_qed: float = 0.05,
    fitness_fn: Optional[FitnessFunction] = None,
    initial_warnings: Optional[List[str]] = None,
    output_limit: Optional[int] = None,
) -> OptimizationSummary:
    if population_size < 1:
        raise ValueError("population_size must be at least 1")
    if not seeds:
        raise ValueError("No seed molecules provided")

    rng = random.Random(seed)
    filters = ChemicalFilters(
        max_molecular_weight=max_molecular_weight,
        max_logp=max_logp,
        max_tpsa=max_tpsa,
        max_hbd=max_hbd,
        max_hba=max_hba,
        min_qed=min_qed,
    )
    warnings: List[str] = list(initial_warnings or [])
    filtered_rows: List[Dict[str, object]] = []

    seed_records = _seed_records_to_population(seeds)
    seed_fps = []
    for seed_entry in seed_records:
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(
                seed_entry["mol"], radius=2, nBits=2048
            )
            seed_fps.append(fp)
        except Exception:
            append_note(seed_entry, "seed_fingerprint_failed")

    for seed_entry in seed_records:
        seed_entry["inherited_or_estimated_fitness"] = _score_molecule(
            seed_entry["mol"], seed_fps, fitness_fn
        )

    population, init_warnings = build_initial_population(
        seed_records, population_size, rng, seed_fps, fitness_fn
    )
    warnings.extend(init_warnings)

    if mutation_rate + crossover_rate > 1.0:
        total = mutation_rate + crossover_rate
        mutation_rate = mutation_rate / total
        crossover_rate = crossover_rate / total
        warnings.append("Mutation and crossover rates normalized to sum to 1.0")

    elite_size = max(1, min(elite_size, population_size))

    for gen in range(1, generations + 1):
        population.sort(
            key=lambda x: float(x.get("inherited_or_estimated_fitness", 0.0)),
            reverse=True,
        )
        elite = population[:elite_size]
        for entry in elite:
            entry["operation"] = "elite"
            entry["generation"] = gen
        next_population: List[Dict[str, object]] = [entry.copy() for entry in elite]

        attempts = 0
        max_attempts = population_size * 6
        while len(next_population) < population_size and attempts < max_attempts:
            attempts += 1
            parent_a = _weighted_choice(population, rng)
            if len(population) > 1:
                parent_b = _weighted_choice(
                    [entry for entry in population if entry is not parent_a], rng
                )
            else:
                parent_b = parent_a

            op_roll = rng.random()
            child_mol: Optional[Chem.Mol] = None
            operation = "clone"
            mutation_failed = False
            crossover_failed = False
            if op_roll < mutation_rate:
                child_mol = mutate(parent_a["mol"], rng)
                operation = "mutation"
                if child_mol is None:
                    mutation_failed = True
            elif op_roll < mutation_rate + crossover_rate:
                child_mol = crossover(
                    parent_a["mol"],
                    parent_b["mol"],
                    rng,
                    float(parent_a.get("inherited_or_estimated_fitness", 0.0)),
                    float(parent_b.get("inherited_or_estimated_fitness", 0.0)),
                )
                operation = "crossover"
                if child_mol is None:
                    crossover_failed = True
            if child_mol is None:
                child_mol = parent_a["mol"]
                operation = "clone"

            best_parent = _best_parent(parent_a, parent_b)
            child_entry: Dict[str, object] = {
                "compound_id": best_parent.get("compound_id", ""),
                "mol": child_mol,
                "operation": operation,
                "generation": gen,
                "parent_1": parent_a.get("source_seed_id", ""),
                "parent_2": parent_b.get("source_seed_id", "") if operation == "crossover" else "",
                "source_seed_id": best_parent.get("source_seed_id", ""),
                "source_final_score": best_parent.get("source_final_score", ""),
                "source_final_classification": best_parent.get(
                    "source_final_classification", ""
                ),
                "source_mode": best_parent.get("source_mode", ""),
                "inherited_or_estimated_fitness": 0.0,
                "valid_molecule": True,
                "notes": "estimated_fitness",
            }
            if mutation_failed:
                append_note(child_entry, "mutation_failed")
            if crossover_failed:
                append_note(child_entry, "crossover_failed")
            if operation == "crossover" and parent_a is parent_b:
                append_note(child_entry, "crossover_same_parent")
            if not child_entry.get("source_seed_id"):
                append_note(child_entry, "missing_seed_match")

            try:
                Chem.SanitizeMol(child_entry["mol"])
            except Exception:
                append_note(child_entry, "invalid_molecule")
                child_entry["valid_molecule"] = False
                filtered_rows.append(build_filtered_row(child_entry, "invalid_molecule"))
                continue

            child_entry["inherited_or_estimated_fitness"] = _score_molecule(
                child_entry["mol"], seed_fps, fitness_fn
            )
            next_population.append(child_entry)

        deduped, removed = _dedup_population(next_population)
        if removed:
            warnings.append(f"Removed {removed} duplicate molecules at generation {gen}")
        deduped.sort(
            key=lambda x: float(x.get("inherited_or_estimated_fitness", 0.0)),
            reverse=True,
        )
        if len(deduped) < population_size:
            warnings.append(f"Population refilled with clones at generation {gen}")
            while len(deduped) < population_size:
                clone = deduped[0].copy()
                append_note(clone, "duplicate_removed")
                deduped.append(clone)
        population = deduped[:population_size]

    population, removed_final = _dedup_population(population)
    if removed_final:
        warnings.append(f"Removed {removed_final} duplicates before output")
    if not population:
        raise ValueError("Population empty after sanitization/dedup")

    population = _try_replenish_unique_population(
        _sort_population_for_output(population),
        population_size,
        rng,
        seed_fps,
        filters,
        warnings,
        filtered_rows,
        fitness_fn,
    )

    run_tag = build_run_id(seed, run_id)
    entries_with_props: List[Dict[str, object]] = []
    csv_rows: List[Dict[str, object]] = []
    for entry in _sort_population_for_output(population):
        if output_limit is not None and len(entries_with_props) >= output_limit:
            break
        try:
            reasons = filter_reasons(entry["mol"], filters)
        except Exception:
            reasons = ["invalid_molecule"]
        if reasons:
            append_note(entry, "filtered_from_output")
            filtered_rows.append(build_filtered_row(entry, ";".join(reasons)))
            continue

        props = build_props(entry, len(entries_with_props) + 1, run_tag)
        entry["props"] = props
        entries_with_props.append(entry)
        csv_rows.append(build_csv_row(entry))

    if not entries_with_props:
        raise ValueError("Population empty after chemical filtering")

    for counter, row in enumerate(filtered_rows, 1):
        if not row.get("generated_id"):
            row["generated_id"] = build_generated_id(run_tag, len(entries_with_props) + counter)

    filtered_csv_path = None
    if out_csv is not None:
        filtered_csv_path = out_csv.parent / "generated_candidates_filtered.csv"
        write_filtered_csv(filtered_rows, filtered_csv_path)
        if not filtered_rows:
            warnings.append("No filtered molecules were generated.")

    unique_smiles = {entry["props"]["canonical_smiles"] for entry in entries_with_props}

    write_sdf(entries_with_props, out_sdf)
    if out_csv is not None:
        write_csv(csv_rows, out_csv)

    return OptimizationSummary(
        seed_selected=len(seeds),
        seed_mapped=len(seeds),
        generated_valid=len(entries_with_props),
        generated_unique=len(unique_smiles),
        filtered=len(filtered_rows),
        filtered_csv=filtered_csv_path,
        warnings=warnings,
    )


def run_genetic_optimization(
    final_candidates_path: Path,
    compounds_path: Path,
    out_sdf: Path,
    out_csv: Optional[Path],
    top_n_seeds: int = 10,
    min_final_score: float = 0.0,
    generations: int = 10,
    population_size: int = 30,
    mutation_rate: float = 0.25,
    crossover_rate: float = 0.50,
    elite_size: int = 5,
    seed: int = 42,
    run_id: Optional[str] = None,
    max_molecular_weight: float = 500.0,
    max_logp: float = 5.0,
    max_tpsa: float = 250.0,
    max_hbd: int = 5,
    max_hba: int = 10,
    min_qed: float = 0.05,
) -> OptimizationSummary:
    seeds, warnings = load_consensus_seeds(
        final_candidates_path,
        compounds_path,
        top_n=top_n_seeds,
        min_final_score=min_final_score,
    )
    return run_genetic_optimization_from_seeds(
        seeds=seeds,
        out_sdf=out_sdf,
        out_csv=out_csv,
        generations=generations,
        population_size=population_size,
        mutation_rate=mutation_rate,
        crossover_rate=crossover_rate,
        elite_size=elite_size,
        seed=seed,
        run_id=run_id,
        max_molecular_weight=max_molecular_weight,
        max_logp=max_logp,
        max_tpsa=max_tpsa,
        max_hbd=max_hbd,
        max_hba=max_hba,
        min_qed=min_qed,
        initial_warnings=warnings,
    )
