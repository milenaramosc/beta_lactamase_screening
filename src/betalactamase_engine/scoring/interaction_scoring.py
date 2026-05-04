import csv
import json
import math
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from betalactamase_engine.utils.pdb_utils import load_structure


@dataclass
class AtomRecord:
    x: float
    y: float
    z: float
    element: str = ""


@dataclass
class ResidueContact:
    chain: str
    residue_name: str
    residue_id: int
    min_distance: float


@dataclass
class ScoringResult:
    compound_id: str
    evaluated_site: str
    site_id: str
    ligand_pose_file: str
    protein_file: str
    target_profile_file: str
    ligand_center: List[float]
    site_center: List[float]
    distance_to_site_center: Optional[float]
    distance_threshold: float
    center_threshold: float
    metal_threshold: float
    nearby_residue_contacts: List[ResidueContact]
    critical_residue_contacts: List[ResidueContact]
    metal_contacts: List[Dict[str, float]]
    contact_summary: Dict[str, int]
    interaction_score: float
    pose_classification: str
    warnings: List[str] = field(default_factory=list)


SUPPORTED_LIGAND_FORMATS = (".pdbqt", ".pdb", ".sdf", ".mol2")


def load_target_profile(path: Path) -> Tuple[Optional[Dict], List[str]]:
    warnings = []
    if not path.exists():
        return None, [f"Target profile not found: {path}"]
    if path.stat().st_size == 0:
        return None, ["Target profile JSON is empty"]
    try:
        with open(path, "r") as handle:
            return json.load(handle), warnings
    except json.JSONDecodeError as exc:
        return None, [f"Target profile JSON is invalid: {exc}"]


def _distance(coord_a: List[float], coord_b: List[float]) -> float:
    return math.sqrt(
        (coord_a[0] - coord_b[0]) ** 2
        + (coord_a[1] - coord_b[1]) ** 2
        + (coord_a[2] - coord_b[2]) ** 2
    )


def _centroid(atoms: List[AtomRecord]) -> Optional[List[float]]:
    if not atoms:
        return None
    return [
        sum(atom.x for atom in atoms) / len(atoms),
        sum(atom.y for atom in atoms) / len(atoms),
        sum(atom.z for atom in atoms) / len(atoms),
    ]


def _parse_pdb_atoms(lines: List[str]) -> List[AtomRecord]:
    atoms = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if len(line) < 54:
            continue
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue
        element = line[76:78].strip() if len(line) >= 78 else ""
        atoms.append(AtomRecord(x=x, y=y, z=z, element=element))
    return atoms


def _parse_pdbqt_atoms(lines: List[str]) -> List[AtomRecord]:
    atoms = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if len(line) < 54:
            continue
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue
        element = ""
        if len(line) >= 78:
            element = line[76:78].strip()
        atoms.append(AtomRecord(x=x, y=y, z=z, element=element))
    return atoms


def _parse_sdf_atoms(lines: List[str]) -> List[AtomRecord]:
    if len(lines) < 4:
        return []
    counts_line = lines[3]
    if len(counts_line) < 6:
        return []
    try:
        atom_count = int(counts_line[0:3].strip())
    except ValueError:
        return []
    atoms = []
    start = 4
    for idx in range(atom_count):
        if start + idx >= len(lines):
            break
        line = lines[start + idx]
        if len(line) < 34:
            continue
        try:
            x = float(line[0:10].strip())
            y = float(line[10:20].strip())
            z = float(line[20:30].strip())
        except ValueError:
            continue
        element = line[31:34].strip()
        atoms.append(AtomRecord(x=x, y=y, z=z, element=element))
    return atoms


def _parse_mol2_atoms(lines: List[str]) -> List[AtomRecord]:
    atoms = []
    in_atoms = False
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>") and in_atoms:
            break
        if not in_atoms:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
        except ValueError:
            continue
        element = parts[5]
        atoms.append(AtomRecord(x=x, y=y, z=z, element=element))
    return atoms


def load_ligand_atoms(path: Path) -> Tuple[List[AtomRecord], List[str]]:
    warnings = []
    if not path.exists():
        return [], [f"Ligand pose not found: {path}"]
    if path.stat().st_size == 0:
        return [], ["Ligand pose file is empty"]
    suffix = path.suffix.lower()
    if suffix not in SUPPORTED_LIGAND_FORMATS:
        return [], [f"Ligand format not supported: {suffix}"]
    try:
        with open(path, "r") as handle:
            lines = handle.readlines()
    except OSError as exc:
        return [], [f"Could not read ligand file: {exc}"]
    parser_map = {
        ".pdb": _parse_pdb_atoms,
        ".pdbqt": _parse_pdbqt_atoms,
        ".sdf": _parse_sdf_atoms,
        ".mol2": _parse_mol2_atoms,
    }
    atoms = parser_map[suffix](lines)
    if not atoms:
        warnings.append("No valid atoms found in ligand pose")
    return atoms, warnings


def load_protein_atoms(pdb_path: Path) -> Tuple[Optional[object], List[str]]:
    warnings = []
    if not pdb_path.exists():
        return None, [f"Protein PDB not found: {pdb_path}"]
    if pdb_path.stat().st_size == 0:
        return None, ["Protein PDB is empty"]
    try:
        structure = load_structure(pdb_path)
    except Exception as exc:
        return None, [f"Error loading protein PDB: {exc}"]
    return structure, warnings


def _extract_candidate_residues(active_site: Dict) -> Dict[Tuple[str, int, str], Dict]:
    residues = {}
    for entry in active_site.get("candidate_catalytic_residues", []) or []:
        try:
            resname = entry["resname"].strip()
            chain = entry.get("chain") or ""
            resid = int(entry["id"])
        except (KeyError, ValueError, TypeError):
            continue
        residues[(chain, resid, resname)] = entry
    return residues


def _collect_metal_coords(profile: Dict) -> List[List[float]]:
    metals = []
    for entry in profile.get("metals", []) or []:
        coord = entry.get("coord")
        if coord and len(coord) == 3:
            metals.append(coord)
    hetero_metals = profile.get("heteroatoms", {}).get("metals", []) or []
    for entry in hetero_metals:
        coord = entry.get("coord")
        if coord and len(coord) == 3:
            metals.append(coord)
    return metals


def _iter_residues(structure, selected_chain: Optional[str]):
    if not selected_chain:
        for residue in structure.get_residues():
            yield residue
        return
    try:
        for residue in structure[0][selected_chain].get_residues():
            yield residue
    except Exception:
        for residue in structure.get_residues():
            yield residue


def _compute_residue_contacts(structure, ligand_atoms: List[AtomRecord], distance_threshold: float, selected_chain: Optional[str]) -> List[ResidueContact]:
    contacts = {}
    for residue in _iter_residues(structure, selected_chain):
        if residue.get_id()[0] != " ":
            continue
        chain_id = residue.get_parent().id
        residue_id = residue.get_id()[1]
        residue_name = residue.get_resname().strip()
        min_distance = None
        for atom in residue.get_atoms():
            coord = atom.get_coord()
            for ligand in ligand_atoms:
                distance = _distance(coord, [ligand.x, ligand.y, ligand.z])
                if distance <= distance_threshold:
                    if min_distance is None or distance < min_distance:
                        min_distance = distance
        if min_distance is not None:
            key = (chain_id, residue_id, residue_name)
            contacts[key] = ResidueContact(
                chain=chain_id,
                residue_name=residue_name,
                residue_id=int(residue_id),
                min_distance=float(min_distance),
            )
    return list(contacts.values())


def _compute_metal_contacts(ligand_atoms: List[AtomRecord], metal_coords: List[List[float]], threshold: float) -> List[Dict[str, float]]:
    contacts = []
    for metal in metal_coords:
        min_distance = None
        for atom in ligand_atoms:
            distance = _distance([atom.x, atom.y, atom.z], metal)
            if distance <= threshold:
                if min_distance is None or distance < min_distance:
                    min_distance = distance
        if min_distance is not None:
            contacts.append({"metal_coord": metal, "min_distance": float(min_distance)})
    return contacts


def _score_interactions(
    distance_to_center: Optional[float],
    center_threshold: float,
    critical_contact_count: int,
    metal_contacts: int,
    metal_present: bool,
    nearby_contact_count: int,
) -> float:
    score = 0.0

    center_component = 0.0
    if distance_to_center is not None and distance_to_center <= center_threshold:
        center_component = max(
            0.0, (center_threshold - distance_to_center) / center_threshold
        )

    catalytic_component = 0.0
    if critical_contact_count >= 2:
        catalytic_component = 1.0
    elif critical_contact_count == 1:
        catalytic_component = 0.625

    contact_component = 0.0
    if nearby_contact_count >= 8:
        contact_component = 1.0
    elif nearby_contact_count >= 4:
        contact_component = 0.6
    elif nearby_contact_count >= 1:
        contact_component = 0.3

    metal_component = 1.0 if metal_contacts > 0 else 0.0

    weights = {
        "center": 0.35,
        "catalytic": 0.40,
        "contacts": 0.10,
        "metal": 0.15,
    }

    if not metal_present:
        redistributed = weights["metal"] / 3.0
        weights["center"] += redistributed
        weights["catalytic"] += redistributed
        weights["contacts"] += redistributed
        weights["metal"] = 0.0

    score += weights["center"] * center_component
    score += weights["catalytic"] * catalytic_component
    score += weights["contacts"] * contact_component
    score += weights["metal"] * metal_component

    return min(score, 1.0)


def _classify_pose(
    distance_to_center: Optional[float],
    center_threshold: float,
    critical_contact_count: int,
    interaction_score: float,
    critical_issue: bool,
) -> str:
    if critical_issue:
        return "unknown"
    if distance_to_center is None:
        return "unknown"
    if distance_to_center > center_threshold:
        return "outside_active_site"
    if interaction_score >= 0.65 and critical_contact_count >= 1:
        return "mechanistically_plausible"
    if 0.35 <= interaction_score < 0.65:
        return "weak_active_site_pose"
    return "weak_active_site_pose"


def score_interactions(
    target_profile_path: Path,
    protein_path: Path,
    ligand_pose_path: Path,
    compound_id: str,
    site: str = "active_site",
    distance_threshold: float = 4.0,
    center_threshold: float = 8.0,
    metal_threshold: float = 3.0,
) -> ScoringResult:
    warnings: List[str] = []
    profile, profile_warnings = load_target_profile(target_profile_path)
    warnings.extend(profile_warnings)

    structure, structure_warnings = load_protein_atoms(protein_path)
    warnings.extend(structure_warnings)

    ligand_atoms, ligand_warnings = load_ligand_atoms(ligand_pose_path)
    warnings.extend(ligand_warnings)

    site_center = None
    selected_chain = None
    candidate_residues = {}
    metal_coords = []

    critical_issue = False

    if profile:
        active_site = profile.get(site, {})
        if not active_site:
            warnings.append(f"Active site entry '{site}' not found in target profile")
            critical_issue = True
        else:
            site_center = active_site.get("center")
            if not site_center:
                warnings.append("Active site center is missing in target profile")
                critical_issue = True
            candidate_residues = _extract_candidate_residues(active_site)
        selected_chain = profile.get("target", {}).get("selected_chain")
        metal_coords = _collect_metal_coords(profile)
    else:
        critical_issue = True

    ligand_center = _centroid(ligand_atoms) if ligand_atoms else None

    distance_to_center = None
    if ligand_center and site_center and len(site_center) == 3:
        distance_to_center = _distance(ligand_center, site_center)
    elif ligand_center and site_center:
        warnings.append("Active site center has invalid coordinates")
        critical_issue = True

    nearby_contacts: List[ResidueContact] = []
    critical_contacts: List[ResidueContact] = []
    metal_contacts: List[Dict[str, float]] = []

    if structure and ligand_atoms:
        nearby_contacts = _compute_residue_contacts(
            structure, ligand_atoms, distance_threshold, selected_chain
        )
        if candidate_residues:
            for contact in nearby_contacts:
                key = (contact.chain, contact.residue_id, contact.residue_name)
                if key in candidate_residues:
                    critical_contacts.append(contact)
        metal_contacts = _compute_metal_contacts(ligand_atoms, metal_coords, metal_threshold)
    else:
        if not structure:
            warnings.append("Protein structure unavailable for contact analysis")
            critical_issue = True
        if not ligand_atoms:
            warnings.append("Ligand atoms unavailable for contact analysis")
            critical_issue = True

    interaction_score = _score_interactions(
        distance_to_center,
        center_threshold,
        len(critical_contacts),
        len(metal_contacts),
        bool(metal_coords),
        len(nearby_contacts),
    )

    if distance_to_center is None:
        critical_issue = True

    pose_classification = _classify_pose(
        distance_to_center,
        center_threshold,
        len(critical_contacts),
        interaction_score,
        critical_issue,
    )

    contact_summary = {
        "nearby_residue_contact_count": len(nearby_contacts),
        "critical_residue_contact_count": len(critical_contacts),
        "metal_contact_count": len(metal_contacts),
    }

    return ScoringResult(
        compound_id=compound_id,
        evaluated_site=site,
        site_id=site,
        ligand_pose_file=str(ligand_pose_path),
        protein_file=str(protein_path),
        target_profile_file=str(target_profile_path),
        ligand_center=ligand_center or [0.0, 0.0, 0.0],
        site_center=site_center or [0.0, 0.0, 0.0],
        distance_to_site_center=distance_to_center,
        distance_threshold=distance_threshold,
        center_threshold=center_threshold,
        metal_threshold=metal_threshold,
        nearby_residue_contacts=nearby_contacts,
        critical_residue_contacts=critical_contacts,
        metal_contacts=metal_contacts,
        contact_summary=contact_summary,
        interaction_score=interaction_score,
        pose_classification=pose_classification,
        warnings=warnings,
    )


def write_json(result: ScoringResult, out_path: Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(result)
    payload["nearby_residue_contacts"] = [asdict(item) for item in result.nearby_residue_contacts]
    payload["critical_residue_contacts"] = [asdict(item) for item in result.critical_residue_contacts]
    with open(out_path, "w") as handle:
        json.dump(payload, handle, indent=2)


def write_csv(result: ScoringResult, out_path: Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "compound_id",
        "evaluated_site",
        "pose_classification",
        "interaction_score",
        "distance_to_site_center",
        "nearby_residue_contact_count",
        "critical_residue_contact_count",
        "metal_contact_count",
        "ligand_pose_file",
        "protein_file",
        "target_profile_file",
        "warnings",
    ]
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(
            {
                "compound_id": result.compound_id,
                "evaluated_site": result.evaluated_site,
                "pose_classification": result.pose_classification,
                "interaction_score": f"{result.interaction_score:.3f}",
                "distance_to_site_center": (
                    f"{result.distance_to_site_center:.3f}"
                    if result.distance_to_site_center is not None
                    else ""
                ),
                "nearby_residue_contact_count": result.contact_summary.get(
                    "nearby_residue_contact_count", 0
                ),
                "critical_residue_contact_count": result.contact_summary.get(
                    "critical_residue_contact_count", 0
                ),
                "metal_contact_count": result.contact_summary.get("metal_contact_count", 0),
                "ligand_pose_file": result.ligand_pose_file,
                "protein_file": result.protein_file,
                "target_profile_file": result.target_profile_file,
                "warnings": "; ".join(result.warnings),
            }
        )
