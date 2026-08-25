import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple


COMPOUND_COLUMNS = ["compound_id", "ligand_id", "molecule_id", "chembl_id", "name", "title"]
VINA_COLUMNS = [
    "vina_score",
    "binding_affinity",
    "affinity",
    "docking_score",
    "score",
    "best_score",
    "dg_kcal_mol",
]


POSE_BONUS = {
    "mechanistically_plausible": 1.0,
    "weak_active_site_pose": 0.7,
    "outside_active_site": 0.2,
    "unknown": 0.4,
}


def _lower_map(columns: List[str]) -> Dict[str, str]:
    return {col.lower(): col for col in columns}


def _detect_column(candidates: List[str], columns: List[str]) -> Optional[str]:
    lower_map = _lower_map(columns)
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    return None


def _parse_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    raw = str(value).strip()
    if not raw:
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _parse_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    raw = str(value).strip()
    if not raw:
        return None
    try:
        return int(raw)
    except ValueError:
        digits = "".join(ch for ch in raw if ch.isdigit())
        if digits:
            try:
                return int(digits)
            except ValueError:
                return None
    return None


def _read_csv_rows(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV file not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"CSV file is empty: {path}")
    with open(path, "r") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV file has no header: {path}")
        rows = list(reader)
        return rows, reader.fieldnames


def _normalize_docking(values: List[float]) -> Tuple[List[float], List[str]]:
    warnings = []
    if not values:
        return [], ["No docking scores found to normalize"]
    best = min(values)
    worst = max(values)
    if best == worst:
        warnings.append("All docking scores are identical; using docking_norm = 0.5")
        return [0.5 for _ in values], warnings
    normalized = [(worst - value) / (worst - best) for value in values]
    return normalized, warnings


def load_ranking(
    path: Path,
    compound_col_override: Optional[str] = None,
    vina_col_override: Optional[str] = None,
) -> Tuple[List[Dict[str, str]], str, str, List[str]]:
    rows, columns = _read_csv_rows(path)
    if compound_col_override:
        if compound_col_override not in columns:
            raise ValueError(
                f"Compound column '{compound_col_override}' not found in ranking."
            )
        compound_col = compound_col_override
    else:
        compound_col = _detect_column(COMPOUND_COLUMNS, columns)

    if vina_col_override:
        if vina_col_override not in columns:
            raise ValueError(f"Vina column '{vina_col_override}' not found in ranking.")
        vina_col = vina_col_override
    else:
        vina_col = _detect_column(VINA_COLUMNS, columns)

    if not compound_col or not vina_col:
        missing = []
        if not compound_col:
            missing.append("compound-column")
        if not vina_col:
            missing.append("vina-column")
        raise ValueError(
            "Could not detect required columns. "
            f"Missing: {', '.join(missing)}. "
            f"Available columns: {', '.join(columns)}."
        )
    return rows, compound_col, vina_col, columns


def _load_interaction_json(path: Path) -> Optional[Dict[str, Optional[float]]]:
    with open(path, "r") as handle:
        data = json.load(handle)
    compound_id = data.get("compound_id")
    if not compound_id:
        return None
    return {
        "compound_id": compound_id,
        "interaction_score": _parse_float(data.get("interaction_score")),
        "pose_classification": data.get("pose_classification") or "unknown",
        "distance_to_site_center": _parse_float(data.get("distance_to_site_center")),
        "critical_residue_contact_count": _parse_int(
            (data.get("contact_summary") or {}).get("critical_residue_contact_count")
        ),
        "nearby_residue_contact_count": _parse_int(
            (data.get("contact_summary") or {}).get("nearby_residue_contact_count")
        ),
        "metal_contact_count": _parse_int(
            (data.get("contact_summary") or {}).get("metal_contact_count")
        ),
    }


def load_interactions(path: Optional[Path]) -> Tuple[Dict[str, Dict[str, Optional[float]]], List[str], bool, List[str]]:
    if not path:
        return {}, ["Interactions not provided; using docking-based ranking."], False, []
    if not path.exists():
        return {}, [f"Interactions path not found: {path}"], False, []

    warnings = []
    duplicates = []
    interactions: Dict[str, Dict[str, Optional[float]]] = {}

    def _record(entry: Dict[str, Optional[float]]):
        compound_id = entry.get("compound_id")
        if not compound_id:
            return
        aliases = [str(compound_id)]
        if "__" in str(compound_id):
            aliases.append(str(compound_id).split("__", 1)[1])
        for alias in aliases:
            alias_entry = dict(entry)
            alias_entry["compound_id"] = alias
            score = alias_entry.get("interaction_score")
            current = interactions.get(alias)
            if current:
                current_score = current.get("interaction_score")
                if score is not None and (current_score is None or score > current_score):
                    interactions[alias] = alias_entry
                duplicates.append(alias)
            else:
                interactions[alias] = alias_entry

    if path.is_dir():
        json_files = sorted(path.glob("*.json"))
        if not json_files:
            return {}, [f"Interactions directory is empty: {path}"], False, []
        for json_path in json_files:
            try:
                entry = _load_interaction_json(json_path)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid interaction JSON: {json_path}: {exc}")
            if entry:
                _record(entry)
    elif path.suffix.lower() == ".json":
        try:
            entry = _load_interaction_json(path)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid interaction JSON: {path}: {exc}")
        if entry:
            _record(entry)
    elif path.suffix.lower() == ".csv":
        rows, columns = _read_csv_rows(path)
        compound_col = _detect_column(["compound_id"], columns)
        if not compound_col:
            raise ValueError(
                f"Interactions CSV missing compound_id column. Available: {', '.join(columns)}"
            )
        for row in rows:
            compound_id = row.get(compound_col, "").strip()
            if not compound_id:
                continue
            entry = {
                "compound_id": compound_id,
                "interaction_score": _parse_float(row.get("interaction_score")),
                "pose_classification": row.get("pose_classification") or "unknown",
                "distance_to_site_center": _parse_float(row.get("distance_to_site_center")),
                "critical_residue_contact_count": _parse_int(
                    row.get("critical_residue_contact_count")
                ),
                "nearby_residue_contact_count": _parse_int(
                    row.get("nearby_residue_contact_count")
                ),
                "metal_contact_count": _parse_int(row.get("metal_contact_count")),
            }
            _record(entry)
    else:
        raise ValueError("Unsupported interactions format. Use CSV, JSON, or a directory of JSONs.")

    if duplicates:
        unique_dupes = sorted(set(duplicates))
        warnings.append(
            "Duplicate interaction entries found for: " + ", ".join(unique_dupes)
        )

    return interactions, warnings, True, duplicates


def _compute_admet_component(row: Dict[str, str]) -> Optional[float]:
    admet_score = _parse_float(row.get("admet_score"))
    if admet_score is not None:
        if admet_score > 1.0:
            admet_score = admet_score / 10.0
        return max(0.0, min(admet_score, 1.0))

    score = 0.5
    qed = _parse_float(row.get("qed"))
    if qed is not None:
        score = 0.5 * score + 0.5 * max(0.0, min(qed, 1.0))

    lipinski = _parse_int(row.get("lipinski_violations"))
    if lipinski:
        score -= min(0.4, 0.1 * lipinski)

    pains = _parse_int(row.get("pains"))
    if pains is None:
        pains = _parse_int(row.get("pains_alerts"))
    if pains and pains > 0:
        score -= 0.2

    score = max(0.0, min(score, 1.0))
    return score


def load_admet(path: Optional[Path]) -> Tuple[Dict[str, Dict[str, Optional[float]]], List[str], bool]:
    if not path:
        return {}, ["ADMET not provided; redistributing weights."], False
    if not path.exists():
        return {}, [f"ADMET file not found: {path}"], False

    rows, columns = _read_csv_rows(path)
    compound_col = _detect_column(COMPOUND_COLUMNS, columns)
    if not compound_col:
        raise ValueError(
            f"ADMET CSV missing compound identifier column. Available: {', '.join(columns)}"
        )

    admet_data = {}
    for row in rows:
        compound_id = row.get(compound_col, "").strip()
        if not compound_id:
            continue
        component = _compute_admet_component(row)
        admet_data[compound_id] = {
            "admet_component": component,
            "raw": row,
        }
    return admet_data, [], True


def compute_consensus(
    ranking_rows: List[Dict[str, str]],
    compound_col: str,
    vina_col: str,
    interactions: Dict[str, Dict[str, Optional[float]]],
    interactions_available: bool,
    interactions_provided: bool,
    admet: Dict[str, Dict[str, Optional[float]]],
    admet_available: bool,
    admet_provided: bool,
    min_interaction_score: Optional[float] = None,
) -> Tuple[List[Dict[str, str]], Dict[str, int], Dict[str, float], List[str]]:
    warnings: List[str] = []
    ranking_ids = {row.get(compound_col, "").strip() for row in ranking_rows}
    ranking_ids.discard("")

    if interactions_available and interactions:
        match_count = sum(1 for cid in interactions.keys() if cid in ranking_ids)
        if match_count == 0:
            warnings.append(
                "Interaction scores were provided, but no compounds matched the docking ranking."
            )
            interactions_available = False
            interactions = {}
    elif interactions_provided and not interactions_available:
        warnings.append(
            "Interaction scores were provided, but no compounds matched the docking ranking."
        )

    if admet_available and admet:
        match_count = sum(1 for cid in admet.keys() if cid in ranking_ids)
        if match_count == 0:
            warnings.append(
                "ADMET file was provided, but no compounds matched the docking ranking."
            )
            admet_available = False
            admet = {}
    elif admet_provided and not admet_available:
        warnings.append(
            "ADMET file was provided, but no compounds matched the docking ranking."
        )
    docking_values = []
    for row in ranking_rows:
        value = _parse_float(row.get(vina_col))
        if value is None:
            raise ValueError(
                f"Non-numeric docking score found in column {vina_col}: {row.get(vina_col)}"
            )
        docking_values.append(value)

    docking_norms, norm_warnings = _normalize_docking(docking_values)
    warnings.extend(norm_warnings)

    weights = {"docking": 0.40, "interaction": 0.35, "admet": 0.20, "pose": 0.05}
    if interactions_available and not admet_available:
        weights = {"docking": 0.50, "interaction": 0.45, "admet": 0.0, "pose": 0.05}
    elif admet_available and not interactions_available:
        weights = {"docking": 0.75, "interaction": 0.0, "admet": 0.20, "pose": 0.05}
    elif not interactions_available and not admet_available:
        weights = {"docking": 0.95, "interaction": 0.0, "admet": 0.0, "pose": 0.05}
        warnings.append("No interaction or ADMET data provided; ranking based mostly on docking.")

    output_rows: List[Dict[str, str]] = []
    counts = {
        "total": 0,
        "with_interactions": 0,
        "with_admet": 0,
        "high": 0,
        "medium": 0,
        "low": 0,
        "deprioritized": 0,
    }

    for idx, row in enumerate(ranking_rows):
        compound_id = row.get(compound_col, "").strip()
        if not compound_id:
            warnings.append("Missing compound_id in ranking row; skipping entry")
            continue

        docking_norm = docking_norms[idx]
        interaction_entry = interactions.get(compound_id)
        admet_entry = admet.get(compound_id)

        interaction_score = None
        pose_classification = "unknown"
        missing_interaction_score = False
        if interactions_available:
            if interaction_entry:
                interaction_score = interaction_entry.get("interaction_score")
                pose_classification = interaction_entry.get("pose_classification") or "unknown"
            else:
                missing_interaction_score = True
                interaction_score = 0.0
                pose_classification = "unknown"
        else:
            interaction_score = 0.0
            if interactions_provided:
                missing_interaction_score = True

        if interaction_score is None:
            interaction_score = 0.0
            missing_interaction_score = True

        if min_interaction_score is not None and interaction_score < min_interaction_score:
            missing_interaction_score = True

        admet_component = None
        missing_admet = False
        if admet_available:
            if admet_entry:
                admet_component = admet_entry.get("admet_component")
            else:
                missing_admet = True
        if admet_component is None:
            admet_component = 0.0
            if admet_available or admet_provided:
                missing_admet = True

        pose_bonus = POSE_BONUS.get(pose_classification, 0.4)

        interaction_component = interaction_score

        final_score = (
            weights["docking"] * docking_norm
            + weights["interaction"] * interaction_component
            + weights["admet"] * admet_component
            + weights["pose"] * pose_bonus
        )

        penalty_notes = []
        if missing_interaction_score and interactions_available:
            penalty_notes.append("missing_interaction_score")
        if missing_admet and admet_available:
            penalty_notes.append("missing_admet")
        if pose_classification == "outside_active_site":
            penalty_notes.append("pose_outside_active_site")
        if pose_classification == "unknown":
            penalty_notes.append("pose_unknown")

        final_score = max(0.0, min(final_score, 1.0))

        if final_score >= 0.75 and (
            pose_classification == "mechanistically_plausible" or interaction_score >= 0.65
        ):
            final_classification = "high_priority_candidate"
        elif final_score >= 0.55:
            final_classification = "medium_priority_candidate"
        elif final_score >= 0.35:
            final_classification = "low_priority_candidate"
        else:
            final_classification = "deprioritized_candidate"

        if pose_classification == "outside_active_site" and interaction_score < 0.35:
            final_classification = "deprioritized_candidate"

        if interactions_available and not missing_interaction_score:
            counts["with_interactions"] += 1
        if admet_available and not missing_admet:
            counts["with_admet"] += 1

        if final_classification == "high_priority_candidate":
            counts["high"] += 1
        elif final_classification == "medium_priority_candidate":
            counts["medium"] += 1
        elif final_classification == "low_priority_candidate":
            counts["low"] += 1
        else:
            counts["deprioritized"] += 1

        counts["total"] += 1

        output = {
            "compound_id": compound_id,
            "vina_score": row.get(vina_col, ""),
            "docking_norm": f"{docking_norm:.3f}",
            "interaction_score": f"{interaction_score:.3f}",
            "pose_classification": pose_classification,
            "interaction_component": f"{interaction_component:.3f}",
            "admet_component": f"{admet_component:.3f}",
            "pose_bonus": f"{pose_bonus:.3f}",
            "final_score": f"{final_score:.3f}",
            "final_classification": final_classification,
            "missing_interaction_score": str(missing_interaction_score).lower(),
            "missing_admet": str(missing_admet).lower(),
            "penalty_notes": ";".join(penalty_notes),
        }

        for key, value in row.items():
            if key not in output:
                output[key] = value
        output_rows.append(output)

    return output_rows, counts, weights, warnings


def sort_consensus_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    class_rank = {
        "high_priority_candidate": 0,
        "medium_priority_candidate": 1,
        "low_priority_candidate": 2,
        "deprioritized_candidate": 3,
    }

    def _sort_key(row: Dict[str, str]):
        final_score = _parse_float(row.get("final_score")) or 0.0
        docking_norm = _parse_float(row.get("docking_norm")) or 0.0
        classification = row.get("final_classification") or "deprioritized_candidate"
        rank = class_rank.get(classification, 3)
        return (-final_score, rank, -docking_norm)

    return sorted(rows, key=_sort_key)


def summarize_rows(
    rows: List[Dict[str, str]],
    interactions_available: bool,
    admet_available: bool,
    interactions_provided: bool,
    admet_provided: bool,
) -> Dict[str, int]:
    counts = {
        "total": 0,
        "with_interactions": 0,
        "with_admet": 0,
        "high": 0,
        "medium": 0,
        "low": 0,
        "deprioritized": 0,
    }
    for row in rows:
        counts["total"] += 1
        classification = row.get("final_classification")
        if classification == "high_priority_candidate":
            counts["high"] += 1
        elif classification == "medium_priority_candidate":
            counts["medium"] += 1
        elif classification == "low_priority_candidate":
            counts["low"] += 1
        else:
            counts["deprioritized"] += 1

        if interactions_available and row.get("missing_interaction_score") == "false":
            counts["with_interactions"] += 1
        if admet_available and row.get("missing_admet") == "false":
            counts["with_admet"] += 1
    if interactions_provided and not interactions_available:
        counts["with_interactions"] = 0
    if admet_provided and not admet_available:
        counts["with_admet"] = 0
    return counts


def write_consensus_csv(rows: List[Dict[str, str]], out_path: Path, extra_fields: List[str]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    required_fields = [
        "compound_id",
        "vina_score",
        "docking_norm",
        "interaction_score",
        "pose_classification",
        "interaction_component",
        "admet_component",
        "pose_bonus",
        "final_score",
        "final_classification",
        "missing_interaction_score",
        "missing_admet",
        "penalty_notes",
    ]
    fieldnames = required_fields + [f for f in extra_fields if f not in required_fields]
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_summary_json(
    out_path: Path,
    counts: Dict[str, int],
    weights: Dict[str, float],
    warnings: List[str],
) -> None:
    summary_path = out_path.parent / "final_candidates_summary.json"
    payload = {
        "total_candidates": counts.get("total", 0),
        "candidates_with_interaction_score": counts.get("with_interactions", 0),
        "candidates_with_admet": counts.get("with_admet", 0),
        "high_priority_count": counts.get("high", 0),
        "medium_priority_count": counts.get("medium", 0),
        "low_priority_count": counts.get("low", 0),
        "deprioritized_count": counts.get("deprioritized", 0),
        "weights_used": weights,
        "warnings": warnings,
    }
    with open(summary_path, "w") as handle:
        json.dump(payload, handle, indent=2)
