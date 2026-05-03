#!/usr/bin/env python3
"""
generate_inhibitors_from_protein.py — GA Guiado por Estrutura da Beta-Lactamase

Analisa a estrutura 3D da beta-lactamase (PDB) e usa características do sítio ativo
para guiar o algoritmo genético na geração de inibidores específicos.

Workflow:
1. Ler PDB da beta-lactamase
2. Identificar sítio ativo (SITE records ou resíduos catalíticos)
3. Analisar características do sítio (carga, hidrofobicidade, tamanho)
4. Criar farmacóforo baseado em inibidores co-cristalizados (se houver)
5. Executar GA com fitness baseada em complementaridade ao sítio ativo

Uso:
    python scripts/generate_inhibitors_from_protein.py --pdb data/structures/1ZG4.pdb
    python scripts/generate_inhibitors_from_protein.py --pdb data/prepared/1ZG4_prepared.pdb --generations 100
"""

import argparse
import csv
import random
import sys
import warnings
from pathlib import Path
from collections import Counter

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, QED
from rdkit.Chem import SDWriter
from rdkit import RDLogger
from rdkit import DataStructs

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR  = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

# ── ANSI colours ──────────────────────────────────────────────────────────────

GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"


# ── Resíduos padrão ───────────────────────────────────────────────────────────

STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
}

# Seeds: inibidores conhecidos
KNOWN_INHIBITORS = {
    "Clavulanato": "O=C(O)[C@H]1N2C(=C(CO)CO[C@@H]12)C(=O)O",
    "Sulbactam": "CC1(C)S[C@@H]2[C@H](NC(=O)C(=O)O)C(=O)N2[C@H]1C(=O)O",
    "Tazobactam": "CN1C(=O)N2[C@H]([C@H](C)S(=O)(=O)N(C)C2(C)C)C1C(=O)O",
    "Avibactam": "C[C@H]1CN(C(=O)N(O1)[C@@H]2CNC(=O)N2)S(=O)(=O)N",
    "Relebactam": "O=C(O)C1NC(=O)N[C@H]1c1cnccn1",
    "Vaborbactam": "CC(C)(C)C(=O)N[C@@H]1B(O)OCC1(C)C",
}


# ── Análise do sítio ativo da proteína ────────────────────────────────────────

def parse_site_records(pdb_path: Path) -> list:
    """Lê registros SITE do PDB."""
    site_residues = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("SITE"):
                continue
            i = 18
            while i + 10 <= len(line.rstrip()):
                group = line[i: i + 11]
                resname = group[0:3].strip()
                chain = group[4].strip() if len(group) > 4 else ""
                resnum_str = group[5:9].strip() if len(group) > 5 else ""
                if resname and chain and resnum_str:
                    try:
                        site_residues.append((resname, chain, int(resnum_str)))
                    except ValueError:
                        pass
                i += 11
    return site_residues


def prompt_catalytic_selection() -> tuple[str, list[str]]:
    if not sys.stdin.isatty():
        print(f"{YELLOW}Aviso: stdin nao interativo. Usando fallback all.{RESET}")
        return "all", []

    def parse_yes_no(value: str) -> str | None:
        value = value.strip().lower()
        if value in {"s", "sim", "y", "yes"}:
            return "yes"
        if value in {"n", "nao", "no"}:
            return "no"
        return None

    for _ in range(3):
        try:
            raw = input("Usar residuos cataliticos conhecidos? (s/n): ")
        except (EOFError, KeyboardInterrupt):
            print(f"{YELLOW}Aviso: entrada interrompida. Usando fallback all.{RESET}")
            return "all", []
        decision = parse_yes_no(raw)
        if decision == "yes":
            for _ in range(3):
                try:
                    res_raw = input("Informe residuos (ex.: SER,LYS,GLU): ")
                except (EOFError, KeyboardInterrupt):
                    print(f"{YELLOW}Aviso: entrada interrompida. Usando fallback all.{RESET}")
                    return "all", []
                tokens = [t.strip().upper() for t in res_raw.split(",") if t.strip()]
                valid = [t for t in tokens if len(t) == 3 and t.isalpha()]
                if valid:
                    print(f"{CYAN}Fallback conhecido: {', '.join(valid)}{RESET}")
                    return "known", valid
                print(f"{YELLOW}Entrada invalida. Tente novamente.{RESET}")
            print(f"{YELLOW}Limite de tentativas. Usando fallback all.{RESET}")
            return "all", []
        if decision == "no":
            print(f"{CYAN}Fallback: usar todos os residuos do PDB.{RESET}")
            return "all", []
        print(f"{YELLOW}Entrada invalida. Tente novamente.{RESET}")

    print(f"{YELLOW}Limite de tentativas. Usando fallback all.{RESET}")
    return "all", []


def identify_catalytic_residues(
    structure,
    site_residues: list,
    fallback_mode: str = "all",
    fallback_resnames: list[str] | None = None,
) -> list:
    """
    Identifica resíduos catalíticos baseado em:
    1. SITE records do PDB
    2. Fallback baseado na escolha do usuário
    """
    model = structure[0]
    catalytic_atoms = []
    if fallback_resnames is None:
        fallback_resnames = []
    
    # 1. Usar SITE records se disponível
    if site_residues:
        for resname, chain, resnum in site_residues:
            try:
                res = model[chain][resnum]
                for atom in res:
                    catalytic_atoms.append(atom)
            except KeyError:
                pass
    
    # 2. Fallback: usar escolha do usuário
    if not catalytic_atoms:
        print(f"{YELLOW}Aviso: SITE records ausentes/invalidos. Usando fallback.{RESET}")
        if fallback_mode == "known":
            fallback_set = {name.upper() for name in fallback_resnames if name}
            for chain in model:
                for res in chain:
                    if res.id[0] != " ":
                        continue
                    if res.get_resname() in fallback_set:
                        for atom in res:
                            catalytic_atoms.append(atom)
        else:
            for chain in model:
                for res in chain:
                    if res.id[0] != " ":
                        continue
                    if res.get_resname() in STANDARD_RESIDUES:
                        for atom in res:
                            catalytic_atoms.append(atom)
    
    return catalytic_atoms


def analyze_binding_site(
    pdb_path: Path,
    radius: float = 10.0,
    fallback_mode: str = "all",
    fallback_resnames: list[str] | None = None,
) -> dict:
    """
    Analisa o sítio de ligação da beta-lactamase.
    
    Retorna dict com:
    - center: coordenadas do centro do sítio
    - volume: volume aproximado (Å³)
    - hydrophobic_ratio: fração de resíduos hidrofóbicos
    - charged_residues: tipos de resíduos carregados
    - h_bond_donors/acceptors: número de doadores/aceptores de H
    - key_residues: resíduos importantes no sítio
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    model = structure[0]
    
    # 1. Identificar sítio ativo
    site_residues = parse_site_records(pdb_path)
    if fallback_resnames is None:
        fallback_resnames = []
    catalytic_atoms = identify_catalytic_residues(
        structure,
        site_residues,
        fallback_mode=fallback_mode,
        fallback_resnames=fallback_resnames,
    )
    
    if not catalytic_atoms:
        print(f"{YELLOW}Aviso: Sítio ativo não identificado. Usando geometria central.{RESET}")
        all_atoms = [atom for chain in model for res in chain for atom in res if atom.element != "H"]
        catalytic_atoms = all_atoms[:50]  # Primeiros 50 átomos como aproximação
    
    # 2. Calcular centro do sítio
    coords = np.array([atom.get_coord() for atom in catalytic_atoms])
    center = coords.mean(axis=0)
    
    # 3. Encontrar resíduos dentro do raio
    all_atoms = [atom for chain in model for res in chain for atom in res]
    ns = NeighborSearch(all_atoms)
    nearby_residues = set()
    
    for cat_atom in catalytic_atoms:
        neighbors = ns.search(cat_atom.get_coord(), radius, level="R")
        nearby_residues.update(neighbors)
    
    # 4. Analisar propriedades dos resíduos
    hydrophobic = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"]
    charged_positive = ["LYS", "ARG", "HIS"]
    charged_negative = ["ASP", "GLU"]
    h_bond_donors_res = ["SER", "THR", "TYR", "ASN", "GLN", "TRP", "HIS", "LYS", "ARG"]
    h_bond_acceptors_res = ["SER", "THR", "TYR", "ASN", "GLN", "ASP", "GLU", "HIS"]
    
    residue_types = [res.get_resname() for res in nearby_residues]
    residue_counts = Counter(residue_types)
    
    hydrophobic_count = sum(residue_counts.get(r, 0) for r in hydrophobic)
    charged_pos_count = sum(residue_counts.get(r, 0) for r in charged_positive)
    charged_neg_count = sum(residue_counts.get(r, 0) for r in charged_negative)
    hbd_count = sum(residue_counts.get(r, 0) for r in h_bond_donors_res)
    hba_count = sum(residue_counts.get(r, 0) for r in h_bond_acceptors_res)
    
    total = len(residue_types)
    
    # 5. Volume aproximado (bounding box)
    res_coords = []
    for res in nearby_residues:
        for atom in res:
            res_coords.append(atom.get_coord())
    res_coords = np.array(res_coords)
    
    if len(res_coords) > 0:
        extents = res_coords.max(axis=0) - res_coords.min(axis=0)
        volume = np.prod(extents)
    else:
        volume = 1000.0  # Default
    
    return {
        "center": center.tolist(),
        "radius": radius,
        "volume": float(volume),
        "num_residues": total,
        "hydrophobic_ratio": hydrophobic_count / total if total > 0 else 0.0,
        "charged_positive": charged_pos_count,
        "charged_negative": charged_neg_count,
        "h_bond_donors": hbd_count,
        "h_bond_acceptors": hba_count,
        "key_residues": dict(residue_counts.most_common(10)),
        "catalytic_residues": [(res.get_resname(), res.get_id()[1]) for res in 
                               set(atom.get_parent() for atom in catalytic_atoms[:20])],
    }


def extract_pharmacophore_from_cocrystal(pdb_path: Path) -> dict | None:
    """
    Extrai farmacóforo de inibidor co-cristalizado (HETATM) se presente.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    model = structure[0]
    
    # Buscar HETATM (exceto água)
    het_residues = []
    for chain in model:
        for res in chain:
            if res.id[0].startswith("H_") and res.get_resname() not in ["HOH", "WAT"]:
                het_residues.append(res)
    
    if not het_residues:
        return None
    
    # Usar primeiro ligante encontrado
    ligand = het_residues[0]
    ligand_atoms = [atom.get_coord() for atom in ligand]
    
    if len(ligand_atoms) < 3:
        return None
    
    ligand_center = np.mean(ligand_atoms, axis=0)
    ligand_size = np.max([np.linalg.norm(coord - ligand_center) for coord in ligand_atoms])
    
    return {
        "name": ligand.get_resname(),
        "center": ligand_center.tolist(),
        "size": float(ligand_size),
        "num_atoms": len(ligand_atoms),
    }


# ── Função de fitness baseada na proteína ─────────────────────────────────────

def fitness_protein_based(
    mol: Chem.Mol,
    site_info: dict,
    ref_fps: list,
    pharmacophore: dict | None = None
) -> float:
    """
    Fitness otimizada baseada nas características do sítio ativo da proteína.
    
    Componentes:
    1. Tamanho/forma compatível com sítio (30%)
    2. Complementaridade química (carga, hidrofobicidade) (25%)
    3. Características estruturais de inibidores (20%)
    4. Drug-likeness (QED) (15%)
    5. Similaridade com inibidores conhecidos (10%)
    """
    try:
        Chem.SanitizeMol(mol)
        
        # 1. Compatibilidade de tamanho (30%)
        mw = Descriptors.MolWt(mol)
        volume_score = 0.0
        
        # Volume ideal baseado no sítio (ajustado)
        site_volume = site_info["volume"]
        ideal_mw_min = max(200, site_volume / 5)  # Heurística
        ideal_mw_max = min(600, site_volume / 2)
        
        if ideal_mw_min <= mw <= ideal_mw_max:
            volume_score = 1.0
        elif mw < ideal_mw_min:
            volume_score = mw / ideal_mw_min
        else:
            volume_score = max(0.0, 1.0 - (mw - ideal_mw_max) / 200)
        
        # 2. Complementaridade química (25%)
        logp = Crippen.MolLogP(mol)
        site_hydrophobic = site_info["hydrophobic_ratio"]
        
        # Se sítio é hidrofóbico, favorece LogP maior
        if site_hydrophobic > 0.5:
            hydro_score = min(1.0, (logp + 2) / 5)  # Favorece 0-3
        else:
            hydro_score = max(0.0, 1.0 - abs(logp) / 3)  # Favorece próximo de 0
        
        # Carga: se sítio tem muitos resíduos carregados, favorece grupos carregados
        charged_site = site_info["charged_positive"] + site_info["charged_negative"]
        
        # Detectar grupos carregados na molécula
        carboxyl = Chem.MolFromSmarts("C(=O)[O;H1,H0-]")
        amine = Chem.MolFromSmarts("[N;H2,H3+]")
        
        has_carboxyl = mol.HasSubstructMatch(carboxyl)
        has_amine = mol.HasSubstructMatch(amine)
        
        charge_score = 0.5  # Neutro
        if charged_site > 3:
            if has_carboxyl or has_amine:
                charge_score = 1.0
        
        chem_complementarity = (hydro_score * 0.6 + charge_score * 0.4)
        
        # 3. Características estruturais de inibidores (20%)
        struct_score = 0.0
        
        # Beta-lactam ou análogos
        beta_lactam = Chem.MolFromSmarts("[C;R1]=1[C;R1][N;R1][C;R1]=1=O")
        if mol.HasSubstructMatch(beta_lactam):
            struct_score += 0.5
        
        # Carboxilato (crítico para beta-lactamases)
        if has_carboxyl:
            struct_score += 0.3
        else:
            struct_score -= 0.1  # Penaliza ausência
        
        # Heterociclos
        heterocycle = Chem.MolFromSmarts("[#7,#8,#16;R]")
        if mol.HasSubstructMatch(heterocycle):
            struct_score += 0.2
        
        struct_score = max(0.0, min(1.0, struct_score))
        
        # 4. Drug-likeness (15%)
        qed_score = QED.qed(mol)
        
        # 5. Similaridade com inibidores conhecidos (10%)
        similarity = 0.0
        if ref_fps:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
            similarity = max(sims) if sims else 0.0
        
        # 6. Bônus de farmacóforo (se ligante co-cristalizado disponível)
        pharmacophore_bonus = 0.0
        if pharmacophore:
            # Se tamanho similar ao ligante co-cristalizado
            ligand_size = pharmacophore["size"]
            mol_size = Descriptors.MolWt(mol) / 50  # Aproximação de raio
            size_diff = abs(mol_size - ligand_size)
            if size_diff < 2.0:
                pharmacophore_bonus = 0.1
        
        # Fitness final
        fitness = (
            volume_score * 0.30 +
            chem_complementarity * 0.25 +
            struct_score * 0.20 +
            qed_score * 0.15 +
            similarity * 0.10 +
            pharmacophore_bonus
        )
        
        return max(0.0, min(1.0, fitness))
    
    except:
        return 0.0


# ── Operadores genéticos (idênticos à versão standalone) ──────────────────────

def mol_to_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol)


def mutate(mol: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    """Mutação química."""
    TRANSFORMS = [
        ("[c:1][F]", "[c:1][Cl]"),
        ("[c:1][Cl]", "[c:1][F]"),
        ("[c:1][Br]", "[c:1][Cl]"),
        ("[c:1][CH3]", "[c:1][OH]"),
        ("[c:1][OH]", "[c:1][NH2]"),
        ("[c:1][NH2]", "[c:1][CH3]"),
        ("[C:1](=O)[OH]", "[C:1](=O)[NH2]"),
        ("[C:1](=O)[NH2]", "[C:1](=O)[OH]"),
        ("[S:1](=O)(=O)[OH]", "[S:1](=O)(=O)[NH2]"),
        ("[c:1][H]", "[c:1][F]"),
    ]
    
    rng.shuffle(TRANSFORMS)
    
    for reactant_smarts, product_smarts in TRANSFORMS:
        try:
            rxn = AllChem.ReactionFromSmarts(f"{reactant_smarts}>>{product_smarts}")
            products = rxn.RunReactants((mol,))
            if products:
                product = products[0][0]
                Chem.SanitizeMol(product)
                return product
        except:
            continue
    
    return None


def crossover(mol_a: Chem.Mol, mol_b: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    """Crossover via BRICS."""
    from rdkit.Chem import BRICS
    
    try:
        frags_a = list(BRICS.BRICSDecompose(mol_a))
        frags_b = list(BRICS.BRICSDecompose(mol_b))
        
        if not frags_a or not frags_b:
            return None
        
        frag_from_b = rng.choice(frags_b)
        combined = frags_a[: len(frags_a) // 2] + [frag_from_b]
        
        frag_mols = [Chem.MolFromSmiles(f) for f in combined if Chem.MolFromSmiles(f)]
        if not frag_mols:
            return None
        
        built = list(BRICS.BRICSBuild(frag_mols))
        if not built:
            return None
        
        candidate = rng.choice(built[:5])
        Chem.SanitizeMol(candidate)
        return candidate
    except:
        return None


# ── Algoritmo Genético ────────────────────────────────────────────────────────

def run_genetic_algorithm(
    seeds: list[tuple[str, Chem.Mol]],
    site_info: dict,
    pharmacophore: dict | None,
    generations: int,
    population_size: int,
    mutation_rate: float,
    rng: random.Random,
) -> list[tuple[float, str, str]]:
    """Executa GA com fitness baseada na proteína."""
    
    # Fingerprints de referência
    ref_fps = []
    for _, mol in seeds:
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            ref_fps.append(fp)
        except:
            pass
    
    # População inicial
    population: list[tuple[float, str, str]] = []
    for name, mol in seeds:
        smi = mol_to_smiles(mol)
        fit = fitness_protein_based(mol, site_info, ref_fps, pharmacophore)
        population.append((fit, smi, f"seed:{name}"))
    
    print(f"\n  População inicial: {len(population)} compostos")
    print(f"  Fitness médio inicial: {sum(f for f, _, _ in population) / len(population):.4f}")
    
    best_seen: dict[str, float] = {smi: fit for fit, smi, _ in population}
    best_ever = max(population, key=lambda x: x[0])
    
    for gen in range(1, generations + 1):
        new_individuals: list[tuple[float, str, str]] = []
        
        population.sort(key=lambda x: x[0], reverse=True)
        elite_size = max(2, population_size // 5)
        elite = population[:elite_size]
        
        attempts = 0
        max_attempts = population_size * 4
        
        while len(new_individuals) < population_size and attempts < max_attempts:
            attempts += 1
            
            if rng.random() < mutation_rate:
                _, smi, _ = rng.choice(elite)
                parent = Chem.MolFromSmiles(smi)
                if parent is None:
                    continue
                child = mutate(parent, rng)
                origin_tag = f"mut_gen{gen}"
            else:
                if len(elite) < 2:
                    continue
                (_, smi_a, _), (_, smi_b, _) = rng.sample(elite, 2)
                mol_a = Chem.MolFromSmiles(smi_a)
                mol_b = Chem.MolFromSmiles(smi_b)
                if mol_a is None or mol_b is None:
                    continue
                child = crossover(mol_a, mol_b, rng)
                origin_tag = f"cross_gen{gen}"
            
            if child is None:
                continue
            
            child_smi = mol_to_smiles(child)
            if child_smi in best_seen:
                continue
            
            fit = fitness_protein_based(child, site_info, ref_fps, pharmacophore)
            best_seen[child_smi] = fit
            new_individuals.append((fit, child_smi, origin_tag))
            
            if fit > best_ever[0]:
                best_ever = (fit, child_smi, origin_tag)
        
        population = elite + new_individuals
        population.sort(key=lambda x: x[0], reverse=True)
        population = population[:population_size]
        
        avg_fit = sum(f for f, _, _ in population) / len(population) if population else 0.0
        best_fit = population[0][0] if population else 0.0
        
        print(f"  Gen {gen:>3}/{generations}  |  Pop: {len(population):>4}  |  "
              f"Melhor: {best_fit:.4f}  |  Média: {avg_fit:.4f}  |  "
              f"Best ever: {best_ever[0]:.4f}", end="\r", flush=True)
    
    print()
    return population


# ── Salvar resultados ─────────────────────────────────────────────────────────

def smiles_to_3d_mol(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    res = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if res != 0:
        res = AllChem.EmbedMolecule(mol, randomSeed=42)
    if res != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def save_results(
    candidates: list[tuple[float, str, str]],
    top_n: int,
    out_sdf: Path,
    out_csv: Path,
) -> int:
    saved = 0
    rows = []
    
    writer = SDWriter(str(out_sdf))
    
    for rank, (fit, smi, origin) in enumerate(candidates[:top_n * 3], 1):
        mol = smiles_to_3d_mol(smi)
        if mol is None:
            continue
        
        mol.SetProp("_Name", f"PROTEIN_GUIDED_{rank:03d}")
        mol.SetProp("SMILES", smi)
        mol.SetProp("FITNESS", f"{fit:.4f}")
        mol.SetProp("ORIGIN", origin)
        mol.SetProp("MW", f"{Descriptors.MolWt(mol):.2f}")
        mol.SetProp("LOGP", f"{Crippen.MolLogP(mol):.2f}")
        mol.SetProp("HBD", str(Lipinski.NumHDonors(mol)))
        mol.SetProp("HBA", str(Lipinski.NumHAcceptors(mol)))
        mol.SetProp("QED", f"{QED.qed(mol):.3f}")
        
        writer.write(mol)
        rows.append({
            "rank": saved + 1,
            "id": f"PROTEIN_GUIDED_{rank:03d}",
            "smiles": smi,
            "fitness": f"{fit:.4f}",
            "mw": f"{Descriptors.MolWt(mol):.2f}",
            "logp": f"{Crippen.MolLogP(mol):.2f}",
            "hbd": Lipinski.NumHDonors(mol),
            "hba": Lipinski.NumHAcceptors(mol),
            "qed": f"{QED.qed(mol):.3f}",
            "origin": origin,
        })
        saved += 1
        if saved >= top_n:
            break
    
    writer.close()
    
    with open(out_csv, "w", newline="") as f:
        if rows:
            writer_csv = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer_csv.writeheader()
            writer_csv.writerows(rows)
    
    return saved


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="GA Guiado por Estrutura da Beta-Lactamase"
    )
    parser.add_argument("--pdb", type=str, required=True,
                        help="Arquivo PDB da beta-lactamase")
    parser.add_argument("--generations", type=int, default=50,
                        help="Número de gerações (padrão: 50)")
    parser.add_argument("--population", type=int, default=100,
                        help="Tamanho da população (padrão: 100)")
    parser.add_argument("--mutation-rate", type=float, default=0.5,
                        help="Taxa de mutação (padrão: 0.5)")
    parser.add_argument("--top-out", type=int, default=20, dest="top_out",
                        help="Quantos candidatos salvar (padrão: 20)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Semente aleatória (padrão: 42)")
    args = parser.parse_args()
    
    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        print(f"{RED}Erro: Arquivo {pdb_path} não encontrado.{RESET}")
        sys.exit(1)
    
    print(f"\n{CYAN}{'═' * 70}{RESET}")
    print(f"{CYAN}  GA Guiado por Estrutura de Beta-Lactamase{RESET}")
    print(f"{CYAN}{'═' * 70}{RESET}\n")
    
    rng = random.Random(args.seed)
    
    # 1. Escolher fallback de resíduos catalíticos
    fallback_mode, fallback_resnames = prompt_catalytic_selection()

    # 2. Analisar sítio ativo
    print(f"Analisando estrutura: {pdb_path.name}")
    site_info = analyze_binding_site(
        pdb_path,
        radius=10.0,
        fallback_mode=fallback_mode,
        fallback_resnames=fallback_resnames,
    )
    
    print(f"\n{BOLD}Características do sítio ativo:{RESET}")
    print(f"  Centro: ({site_info['center'][0]:.2f}, {site_info['center'][1]:.2f}, {site_info['center'][2]:.2f})")
    print(f"  Volume: {site_info['volume']:.1f} Å³")
    print(f"  Resíduos no sítio: {site_info['num_residues']}")
    print(f"  Hidrofobicidade: {site_info['hydrophobic_ratio']:.2%}")
    print(f"  Resíduos carregados: +{site_info['charged_positive']} / -{site_info['charged_negative']}")
    print(f"  Doadores H-bond: {site_info['h_bond_donors']}")
    print(f"  Aceptores H-bond: {site_info['h_bond_acceptors']}")
    
    if site_info['catalytic_residues']:
        print(f"\n{BOLD}Resíduos catalíticos identificados:{RESET}")
        for resname, resnum in site_info['catalytic_residues'][:5]:
            print(f"    {resname}{resnum}")
    
    # 3. Buscar farmacóforo (ligante co-cristalizado)
    pharmacophore = extract_pharmacophore_from_cocrystal(pdb_path)
    if pharmacophore:
        print(f"\n{GREEN}✓ Ligante co-cristalizado encontrado: {pharmacophore['name']}{RESET}")
        print(f"  Tamanho: {pharmacophore['size']:.2f} Å")
        print(f"  Átomos: {pharmacophore['num_atoms']}")
    else:
        print(f"\n{YELLOW}! Nenhum ligante co-cristalizado encontrado.{RESET}")
    
    # 4. Carregar seeds
    print(f"\n{BOLD}Carregando inibidores conhecidos como seeds:{RESET}")
    seeds = []
    for name, smiles in KNOWN_INHIBITORS.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            seeds.append((name, mol))
            print(f"  ✓ {name}")
    
    # 5. Executar GA
    print(f"\n{BOLD}Parâmetros do GA:{RESET}")
    print(f"  Gerações: {args.generations}")
    print(f"  População: {args.population}")
    print(f"  Taxa de mutação: {args.mutation_rate}")
    
    candidates = run_genetic_algorithm(
        seeds,
        site_info,
        pharmacophore,
        generations=args.generations,
        population_size=args.population,
        mutation_rate=args.mutation_rate,
        rng=rng,
    )
    
    # 6. Salvar resultados
    results_dir = PROJECT_DIR / "results"
    results_dir.mkdir(exist_ok=True)
    out_sdf = results_dir / f"inhibitors_{pdb_path.stem}.sdf"
    out_csv = results_dir / f"inhibitors_{pdb_path.stem}.csv"
    
    print(f"\nGerando coordenadas 3D e salvando top {args.top_out} candidatos...")
    saved = save_results(candidates, args.top_out, out_sdf, out_csv)
    
    # 7. Resumo
    print(f"\n{GREEN}{'═' * 70}{RESET}")
    print(f"{GREEN}Concluído!{RESET}")
    print(f"  {saved} candidatos salvos:")
    print(f"    {out_sdf}")
    print(f"    {out_csv}")
    
    if saved > 0:
        print(f"\n{BOLD}Top 10 candidatos gerados:{RESET}")
        print(f"  {'#':<4} {'Fitness':<10} {'MW':<8} {'LogP':<7} {'QED':<7} {'SMILES':<60}")
        print(f"  {'─'*4} {'─'*10} {'─'*8} {'─'*7} {'─'*7} {'─'*60}")
        with open(out_csv) as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader, 1):
                if i > 10:
                    break
                smi_short = row["smiles"][:57] + "..." if len(row["smiles"]) > 60 else row["smiles"]
                print(f"  {row['rank']:<4} {row['fitness']:<10} {row['mw']:<8} "
                      f"{row['logp']:<7} {row['qed']:<7} {smi_short}")
    
    print(f"\n{CYAN}Próximos passos:{RESET}")
    print(f"  1. Revisar candidatos em {out_csv}")
    print(f"  2. Validar com docking: python scripts/screen_and_rank.py --ligands {out_sdf}")
    print(f"  3. Sintetizar e testar experimentalmente\n")


if __name__ == "__main__":
    main()
