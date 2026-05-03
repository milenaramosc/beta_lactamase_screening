import json
import math
from pathlib import Path
from typing import List, Dict, Any
from dataclasses import dataclass, asdict, field

from betalactamase_engine.utils.pdb_utils import load_structure, calculate_centroid

@dataclass
class TargetProfile:
    target: Dict[str, Any] = field(default_factory=dict)
    classification: Dict[str, Any] = field(default_factory=dict)
    active_site: Dict[str, Any] = field(default_factory=dict)
    metals: List[Dict[str, Any]] = field(default_factory=list)
    heteroatoms: Dict[str, List[Any]] = field(default_factory=lambda: {
        "waters": [],
        "metals": [],
        "organic_ligands": [],
        "others": []
    })
    pocket_properties: Dict[str, Any] = field(default_factory=dict)
    quality_control: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

class TargetProfiler:
    def __init__(self, pdb_path: Path, radius: float = 8.0, chain: str = None):
        self.pdb_path = Path(pdb_path)
        self.radius = radius
        self.chain = chain.strip() if isinstance(chain, str) and chain.strip() else None
        self.structure = None
        self.profile = TargetProfile()
        self.profile.quality_control = {
            "has_site_records": False,
            "has_organic_ligand": False,
            "has_metals": False,
            "used_fallback": False,
            "notes": []
        }
        self.selected_chain = None
        self._water_entries = []
        
        # Typical residues for serine-beta-lactamases
        self.CATALYTIC_RES_TYPES = ["SER", "LYS", "GLU", "ASN", "ARG", "ASP"]
        # Typical metals
        self.METAL_TYPES = ["ZN", "MG", "MN", "FE", "FE2", "FE3", "CA", "CO", "NI", "CU", "NA", "K"]

    def run(self) -> TargetProfile:
        if not self.pdb_path.exists():
            self.profile.warnings.append(f"PDB file not found: {self.pdb_path}")
            return self.profile

        if self.pdb_path.stat().st_size == 0:
            self.profile.warnings.append("PDB file is empty")
            return self.profile

        try:
            self.structure = load_structure(self.pdb_path)
        except Exception as e:
            self.profile.warnings.append(f"Error loading PDB structure: {str(e)}")
            return self.profile

        self._extract_basic_info()
        if self.profile.target.get("atom_count", 0) == 0:
            self.profile.warnings.append("No atoms found in PDB")
            self._set_default_active_site()
            return self.profile

        self._detect_heteroatoms()
        self._detect_metals()
        self._estimate_active_site()
        self._summarize_waters()
        self._characterize_pocket()
        self._infer_class()
        
        return self.profile

    def _extract_basic_info(self):
        chains = [chain.id for chain in self.structure.get_chains()]
        atoms = list(self.structure.get_atoms())
        residues = list(self.structure.get_residues())
        self._select_chain(chains)

        self.profile.target = {
            "pdb_file": str(self.pdb_path),
            "pdb_name": self.pdb_path.stem,
            "chains": chains,
            "chain_count": len(chains),
            "atom_count": len(atoms),
            "residue_count": len(residues),
            "selected_chain": self.selected_chain
        }

    def _select_chain(self, chains: List[str]):
        if not chains:
            self.selected_chain = None
            self.profile.warnings.append("No chains detected in PDB")
            self.profile.quality_control = {
                "has_site_records": False,
                "has_organic_ligand": False,
                "has_metals": False,
                "used_fallback": True,
                "notes": ["No chains available for selection"]
            }
            return

        sorted_chains = sorted(chains)
        if self.chain:
            if self.chain in chains:
                self.selected_chain = self.chain
            else:
                self.selected_chain = sorted_chains[0]
                self.profile.warnings.append(f"Requested chain {self.chain} not found; using {self.selected_chain}")
                self.profile.quality_control["used_fallback"] = True
                self.profile.quality_control["notes"].append(
                    f"Requested chain {self.chain} not found; selected {self.selected_chain}"
                )
        else:
            if len(sorted_chains) > 1:
                self.selected_chain = sorted_chains[0]
                self.profile.quality_control["notes"].append(
                    f"Multiple chains detected; selected {self.selected_chain} automatically"
                )
            else:
                self.selected_chain = sorted_chains[0]

    def _detect_heteroatoms(self):
        waters = []
        metals = []
        organic_ligands = []
        others = []
        
        for residue in self._iter_residues():
            res_name = residue.get_resname().strip()
            res_id = residue.get_id()
            het_flag = res_id[0]
            chain_id = residue.get_parent().id

            if het_flag == " ":
                continue  # Normal amino acid

            entry = {
                "name": res_name,
                "chain": chain_id,
                "id": res_id[1],
                "icode": res_id[2].strip() if isinstance(res_id[2], str) else ""
            }

            if het_flag == "W" or res_name in {"HOH", "WAT"}:
                entry["coord"] = calculate_centroid(list(residue.get_atoms()))
                waters.append(entry)
                continue

            if het_flag.startswith("H_"):
                is_metal = res_name in self.METAL_TYPES
                if is_metal:
                    entry["coord"] = calculate_centroid(list(residue.get_atoms()))
                    entry["hetflag"] = het_flag
                    metals.append(entry)
                else:
                    entry["atoms"] = len(list(residue.get_atoms()))
                    entry["hetflag"] = het_flag
                    organic_ligands.append(entry)
            else:
                others.append(entry)

        self._water_entries = waters
        self.profile.heteroatoms = {
            "waters": [],
            "metals": metals,
            "organic_ligands": organic_ligands,
            "others": others
        }

        if not organic_ligands:
            self.profile.warnings.append("No organic ligands detected in PDB")
        else:
            self.profile.quality_control["has_organic_ligand"] = True

    def _detect_metals(self):
        # Already done in _detect_heteroatoms, just linking for clarity
        self.profile.metals = self.profile.heteroatoms["metals"]

        if not self.profile.metals:
            self.profile.warnings.append("No metals detected in PDB")
        else:
            self.profile.quality_control["has_metals"] = True

    def _set_default_active_site(self):
        self.profile.active_site = {
            "center": [0.0, 0.0, 0.0],
            "detection_method": "GEOMETRIC_FALLBACK",
            "radius_angstrom": self.radius,
            "candidate_catalytic_residues": []
        }

    def _read_site_records(self):
        site_residues = []
        try:
            with open(self.pdb_path, "r") as handle:
                for line in handle:
                    if not line.startswith("SITE "):
                        continue
                    site_residues.extend(self._parse_site_line(line))
        except OSError as exc:
            self.profile.warnings.append(f"Could not read PDB file for SITE records: {exc}")
            return []

        if not site_residues:
            self.profile.warnings.append("No SITE records found in PDB")
        else:
            self.profile.quality_control["has_site_records"] = True

        return site_residues

    @staticmethod
    def _parse_site_line(line: str):
        residues = []
        slots = [
            (18, 21, 22, 23, 27, 27),
            (29, 32, 33, 34, 38, 38),
            (40, 43, 44, 45, 49, 49),
            (51, 54, 55, 56, 60, 60)
        ]

        for res_start, res_end, chain_idx, seq_start, seq_end, icode_idx in slots:
            if len(line) < seq_end:
                continue
            resname = line[res_start:res_end].strip()
            chain_id = line[chain_idx].strip()
            resseq = line[seq_start:seq_end].strip()
            icode = line[icode_idx].strip() if len(line) > icode_idx else ""
            if not resname or not resseq:
                continue
            try:
                resseq_int = int(resseq)
            except ValueError:
                continue
            residues.append({
                "resname": resname,
                "chain": chain_id,
                "resseq": resseq_int,
                "icode": icode
            })

        return residues

    def _iter_residues(self):
        if not self.selected_chain:
            return self.structure.get_residues()
        try:
            return self.structure[0][self.selected_chain].get_residues()
        except Exception:
            return self.structure.get_residues()

    def _iter_atoms(self):
        if not self.selected_chain:
            return self.structure.get_atoms()
        try:
            return self.structure[0][self.selected_chain].get_atoms()
        except Exception:
            return self.structure.get_atoms()

    @staticmethod
    def _distance(coord_a, coord_b):
        return math.sqrt(
            (coord_a[0] - coord_b[0]) ** 2
            + (coord_a[1] - coord_b[1]) ** 2
            + (coord_a[2] - coord_b[2]) ** 2
        )

    def _estimate_active_site(self):
        center = [0.0, 0.0, 0.0]
        method = "GEOMETRIC_FALLBACK"

        # 1. Try SITE records
        site_residues = self._read_site_records()
        if site_residues:
            site_atoms = []
            for residue_data in site_residues:
                chain_id = residue_data["chain"]
                resseq = residue_data["resseq"]
                icode = residue_data["icode"] or " "
                if self.selected_chain and chain_id and chain_id != self.selected_chain:
                    continue
                try:
                    residue = self.structure[0][chain_id][(" ", resseq, icode)]
                except Exception:
                    continue
                site_atoms.extend(list(residue.get_atoms()))

            if site_atoms:
                centroid = calculate_centroid(site_atoms)
                if centroid:
                    center = centroid
                    method = "SITE_RECORD"
            else:
                self.profile.warnings.append("SITE records found but residues could not be resolved")
                self.profile.quality_control["used_fallback"] = True

        # 2. Try organic ligand
        if method == "GEOMETRIC_FALLBACK" and self.profile.heteroatoms["organic_ligands"]:
            largest_ligand = max(self.profile.heteroatoms["organic_ligands"], key=lambda x: x.get("atoms", 0))
            hetflag = largest_ligand.get("hetflag") or f"H_{largest_ligand['name']}"
            try:
                res = self.structure[0][largest_ligand["chain"]][(hetflag, largest_ligand["id"], largest_ligand.get("icode", " ") or " ")]
                centroid = calculate_centroid(list(res.get_atoms()))
            except Exception:
                centroid = None

            if centroid:
                center = centroid
                method = "COCRYSTAL_LIGAND"
        
        # 3. Try metal center (if ZN is present, it's likely a Metallo-beta-lactamase site)
        elif method == "GEOMETRIC_FALLBACK" and self.profile.metals:
            coords = [m["coord"] for m in self.profile.metals if m.get("coord")]
            if coords:
                center = [sum(c) / len(coords) for c in zip(*coords)]
                method = "METAL_CENTER"
            
        # 4. Try catalytic residues
        if method == "GEOMETRIC_FALLBACK":
            catalytic_candidates = []
            for residue in self._iter_residues():
                if residue.get_resname() in self.CATALYTIC_RES_TYPES and residue.get_id()[0] == " ":
                    catalytic_candidates.append(residue)

            if catalytic_candidates:
                ser_residues = [r for r in catalytic_candidates if r.get_resname() == "SER"]
                best_pair = None
                best_distance = None

                for ser in ser_residues:
                    ser_centroid = calculate_centroid(list(ser.get_atoms()))
                    if not ser_centroid:
                        continue
                    for partner in catalytic_candidates:
                        if partner is ser:
                            continue
                        if partner.get_resname() not in {"LYS", "GLU", "ASN", "ARG", "ASP"}:
                            continue
                        partner_centroid = calculate_centroid(list(partner.get_atoms()))
                        if not partner_centroid:
                            continue
                        distance = self._distance(ser_centroid, partner_centroid)
                        if distance <= 6.0 and (best_distance is None or distance < best_distance):
                            best_distance = distance
                            best_pair = (ser, partner)

                if best_pair:
                    atoms = list(best_pair[0].get_atoms()) + list(best_pair[1].get_atoms())
                    centroid = calculate_centroid(atoms)
                    if centroid:
                        center = centroid
                        method = "CATALYTIC_RESIDUES"

            if method == "GEOMETRIC_FALLBACK":
                atoms = list(self._iter_atoms())
                centroid = calculate_centroid(atoms)
                if centroid:
                    center = centroid
        
        self.profile.active_site = {
            "center": center,
            "detection_method": method,
            "radius_angstrom": self.radius,
            "candidate_catalytic_residues": [] # Will populate in pocket characterization
        }

        if method == "GEOMETRIC_FALLBACK":
            self.profile.quality_control["used_fallback"] = True
            self.profile.quality_control["notes"].append("Active site estimated with geometric fallback")

    def _summarize_waters(self):
        center = self.profile.active_site.get("center", [0.0, 0.0, 0.0])
        near = []
        for water in self._water_entries:
            coord = water.get("coord")
            if not coord:
                continue
            distance = self._distance(coord, center)
            if distance <= self.radius:
                near.append({
                    "name": water["name"],
                    "chain": water["chain"],
                    "id": water["id"],
                    "icode": water.get("icode", ""),
                    "distance": float(distance)
                })

        self.profile.heteroatoms["waters"] = {
            "count": len(self._water_entries),
            "near_active_site_count": len(near),
            "near_active_site": near
        }

    def _characterize_pocket(self):
        center = self.profile.active_site["center"]
        radius = self.radius
        
        nearby_residues = []
        
        # Residue classifications
        hydrophobic = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"]
        polar = ["SER", "THR", "CYS", "TYR", "ASN", "GLN"]
        positive = ["LYS", "ARG", "HIS"]
        negative = ["ASP", "GLU"]
        aromatic = ["PHE", "TYR", "TRP", "HIS"]
        
        counts = {
            "hydrophobic": 0,
            "polar": 0,
            "positive": 0,
            "negative": 0,
            "aromatic": 0
        }
        
        catalytic_candidates = []

        for residue in self._iter_residues():
            if residue.get_id()[0] != " ":
                continue # Skip heteroatoms
            
            res_centroid = calculate_centroid(list(residue.get_atoms()))
            if not res_centroid:
                continue
            distance = self._distance(res_centroid, center)
            
            if distance <= radius:
                res_name = residue.get_resname()
                nearby_residues.append(residue)
                
                if res_name in hydrophobic: counts["hydrophobic"] += 1
                if res_name in polar: counts["polar"] += 1
                if res_name in positive: counts["positive"] += 1
                if res_name in negative: counts["negative"] += 1
                if res_name in aromatic: counts["aromatic"] += 1
                
                if res_name in self.CATALYTIC_RES_TYPES:
                    catalytic_candidates.append({
                        "resname": res_name,
                        "chain": residue.get_parent().id,
                        "id": residue.get_id()[1],
                        "distance": float(distance)
                    })

        total = len(nearby_residues)
        if total > 0:
            props = {
                "nearby_residue_count": total,
                "hydrophobic_ratio": counts["hydrophobic"] / total,
                "polar_ratio": counts["polar"] / total,
                "charged_ratio": (counts["positive"] + counts["negative"]) / total,
                "aromatic_ratio": counts["aromatic"] / total,
                "positive_residue_count": counts["positive"],
                "negative_residue_count": counts["negative"],
                "polar_residue_count": counts["polar"],
                "hydrophobic_residue_count": counts["hydrophobic"],
                "aromatic_residue_count": counts["aromatic"]
            }
        else:
            props = {
                "nearby_residue_count": 0,
                "hydrophobic_ratio": 0.0,
                "polar_ratio": 0.0,
                "charged_ratio": 0.0,
                "aromatic_ratio": 0.0,
                "positive_residue_count": 0,
                "negative_residue_count": 0,
                "polar_residue_count": 0,
                "hydrophobic_residue_count": 0,
                "aromatic_residue_count": 0
            }
            
        self.profile.pocket_properties = props
        self.profile.active_site["candidate_catalytic_residues"] = catalytic_candidates

    def _infer_class(self):
        predicted_class = "unknown_beta_lactamase_class"
        confidence = "low"
        evidence = []
        
        # Check for metals (Class B)
        metal_present = any(m["name"] in self.METAL_TYPES for m in self.profile.metals)
        if metal_present:
            predicted_class = "possible_class_B_metallo_beta_lactamase"
            confidence = "medium"
            evidence.append("Possible metal ions detected in structure")

            # Check if ZN is near the active site center
            center = self.profile.active_site["center"]
            metal_near_site = False
            for m in self.profile.metals:
                if m.get("coord") and self._distance(m["coord"], center) < 5.0:
                    metal_near_site = True
                    break

            if metal_near_site:
                confidence = "high"
                evidence.append("Metal ion possibly located near active site center")

        # Check for serine-beta-lactamase motifs (Classes A, C, D)
        # Typically: SER-X-X-LYS (not easy with distance only, but we can check cluster)
        else:
            catalytic = self.profile.active_site["candidate_catalytic_residues"]
            has_ser = any(c["resname"] == "SER" for c in catalytic)
            has_lys = any(c["resname"] == "LYS" for c in catalytic)
            
            if has_ser and has_lys:
                predicted_class = "possible_serine_beta_lactamase"
                confidence = "medium"
                evidence.append("Active site shows candidate Serine and Lysine residues")
                
                # Check for Class A specific (Glu166 proxy)
                has_glu = any(c["resname"] == "GLU" for c in catalytic)
                if has_glu:
                    confidence = "high"
                    evidence.append("Glutamic acid possibly near Ser/Lys cluster")
            elif has_ser:
                predicted_class = "possible_serine_beta_lactamase"
                confidence = "low"
                evidence.append("Candidate serine residue detected in active site")

        self.profile.classification = {
            "predicted_class": predicted_class,
            "confidence": confidence,
            "evidence": evidence
        }

    def save_json(self, out_path: Path):
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        with open(out_path, "w") as f:
            json.dump(asdict(self.profile), f, indent=2)
