#!/usr/bin/env python3
"""
generate_candidate.py — Etapa 7: Geração de molécula inibidora de novo

Usa os top compostos do docking (results/ranking.csv) como população inicial
e aplica um algoritmo genético simples via RDKit para gerar moléculas novas
com melhor potencial de ligação à beta-lactamase.

Uso:
    python scripts/generate_candidate.py
    python scripts/generate_candidate.py --generations 50 --population 100 --top-seed 20

Saída:
    results/generated_candidates.sdf   — moléculas geradas (top-5)
    results/generated_candidates.csv   — resumo com SMILES, fitness e origem
"""

import argparse
import csv
import random
import sys
import warnings
from pathlib import Path

import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import SDWriter
from rdkit import RDLogger

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


# ── Config ────────────────────────────────────────────────────────────────────

def load_config() -> dict:
    p = PROJECT_DIR / "config.yaml"
    with open(p) as f:
        return yaml.safe_load(f)


# ── Load seed compounds from ranking.csv + SDF ───────────────────────────────

def load_seed_compounds(top_n: int) -> list[tuple[str, Chem.Mol, float]]:
    """
    Lê os top_n compostos do ranking.csv e busca seus SMILES no SDF.
    Retorna lista de (compound_id, mol, dg_kcal_mol).
    """
    ranking_csv = PROJECT_DIR / "results" / "ranking.csv"
    sdf_path    = PROJECT_DIR / "data"    / "compounds" / "compounds.sdf"

    if not ranking_csv.exists():
        print(f"{RED}Erro:{RESET} results/ranking.csv não encontrado. Execute screen_and_rank.py primeiro.")
        sys.exit(1)

    if not sdf_path.exists():
        print(f"{RED}Erro:{RESET} data/compounds/compounds.sdf não encontrado. Execute fetch_compounds.py primeiro.")
        sys.exit(1)

    # Lê os top_n IDs do ranking
    top_ids: dict[str, float] = {}
    with open(ranking_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if len(top_ids) >= top_n:
                break
            top_ids[row["compound_id"]] = float(row["dg_kcal_mol"])

    # Carrega moléculas do SDF
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    seeds = []
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        if name in top_ids:
            dg = top_ids[name]
            seeds.append((name, mol, dg))
        if len(seeds) >= top_n:
            break

    print(f"  {len(seeds)} compostos semente carregados do ranking.")
    return seeds


# ── Fitness function ──────────────────────────────────────────────────────────

def lipinski_score(mol: Chem.Mol) -> float:
    """
    Retorna 1.0 se passa em todas as regras de Lipinski, penaliza por violações.
    """
    violations = 0
    if Descriptors.MolWt(mol) > 500:
        violations += 1
    if Descriptors.MolLogP(mol) > 5:
        violations += 1
    if Descriptors.NumHDonors(mol) > 5:
        violations += 1
    if Descriptors.NumHAcceptors(mol) > 10:
        violations += 1
    return max(0.0, 1.0 - violations * 0.25)


def tanimoto_to_reference(mol: Chem.Mol, ref_fps: list) -> float:
    """
    Retorna a similaridade de Tanimoto máxima em relação aos compostos de referência.
    """
    if not ref_fps:
        return 0.0
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    from rdkit import DataStructs
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    return max(sims) if sims else 0.0


def fitness(mol: Chem.Mol, ref_fps: list) -> float:
    """
    Função de fitness combinada:
      - Similaridade com compostos de referência (alto = estrutura promissora)
      - Penalidade por violações de Lipinski
      - Penalidade por tamanho excessivo
    """
    lip   = lipinski_score(mol)
    sim   = tanimoto_to_reference(mol, ref_fps)
    mw    = Descriptors.MolWt(mol)
    size_pen = max(0.0, (mw - 400) / 200)  # penaliza MW > 400

    return sim * lip - size_pen * 0.1


# ── Genetic operators ─────────────────────────────────────────────────────────

def mol_to_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol)


def smiles_to_3d_mol(smiles: str) -> Chem.Mol | None:
    """Converte SMILES para mol com coordenadas 3D."""
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


def mutate(mol: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    """
    Aplica uma mutação simples: troca um átomo de halogênio ou grupo metil
    usando transformações SMARTS do RDKit.
    Retorna nova molécula ou None se a mutação falhou.
    """
    from rdkit.Chem import AllChem

    # Lista de transformações (SMARTS de substituição simples)
    TRANSFORMS = [
        ("[c:1][F]", "[c:1][Cl]"),
        ("[c:1][Cl]", "[c:1][F]"),
        ("[c:1][Br]", "[c:1][Cl]"),
        ("[c:1][CH3]", "[c:1][OH]"),
        ("[c:1][OH]", "[c:1][NH2]"),
        ("[c:1][NH2]", "[c:1][CH3]"),
        ("[C:1](=O)[OH]", "[C:1](=O)[NH2]"),
    ]

    rng.shuffle(TRANSFORMS)
    smiles = mol_to_smiles(mol)

    for reactant_smarts, product_smarts in TRANSFORMS:
        try:
            rxn = AllChem.ReactionFromSmarts(f"{reactant_smarts}>>{product_smarts}")
            products = rxn.RunReactants((mol,))
            if products:
                product = products[0][0]
                try:
                    Chem.SanitizeMol(product)
                    return product
                except Exception:
                    continue
        except Exception:
            continue

    return None


def crossover(mol_a: Chem.Mol, mol_b: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    """
    Crossover via fragmentação BRICS: troca fragmentos entre dois compostos.
    Retorna nova molécula ou None se falhou.
    """
    from rdkit.Chem import BRICS

    try:
        frags_a = list(BRICS.BRICSDecompose(mol_a))
        frags_b = list(BRICS.BRICSDecompose(mol_b))

        if not frags_a or not frags_b:
            return None

        # Troca um fragmento aleatório
        frag_from_b = rng.choice(frags_b)
        combined = frags_a[: len(frags_a) // 2] + [frag_from_b]

        # Reconstrói via BRICS.BRICSBuild
        frag_mols = [Chem.MolFromSmiles(f) for f in combined if Chem.MolFromSmiles(f)]
        if not frag_mols:
            return None

        built = list(BRICS.BRICSBuild(frag_mols))
        if not built:
            return None

        candidate = rng.choice(built[:5])  # pega um dos primeiros
        Chem.SanitizeMol(candidate)
        return candidate
    except Exception:
        return None


# ── Genetic algorithm ─────────────────────────────────────────────────────────

def run_genetic_algorithm(
    seeds: list[tuple[str, Chem.Mol, float]],
    generations: int,
    population_size: int,
    rng: random.Random,
) -> list[tuple[float, str, str]]:
    """
    Executa o algoritmo genético.
    Retorna lista de (fitness, smiles, origin) ordenada por fitness decrescente.
    """
    # Fingerprints de referência para a função de fitness
    ref_fps = []
    for _, mol, _ in seeds:
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            ref_fps.append(fp)
        except Exception:
            pass

    # População inicial = SMILES dos seeds
    population: list[tuple[float, str, str]] = []
    for name, mol, dg in seeds:
        smi = mol_to_smiles(mol)
        fit = fitness(mol, ref_fps)
        population.append((fit, smi, f"seed:{name}"))

    print(f"  População inicial: {len(population)} compostos")

    best_seen: dict[str, float] = {smi: fit for fit, smi, _ in population}

    for gen in range(1, generations + 1):
        new_individuals: list[tuple[float, str, str]] = []

        # Ordenar por fitness
        population.sort(key=lambda x: x[0], reverse=True)
        elite = population[: max(2, population_size // 5)]

        attempts = 0
        while len(new_individuals) < population_size and attempts < population_size * 4:
            attempts += 1
            op = rng.choice(["mutate", "crossover"])

            try:
                if op == "mutate" and elite:
                    _, smi, origin = rng.choice(elite)
                    parent = Chem.MolFromSmiles(smi)
                    if parent is None:
                        continue
                    child = mutate(parent, rng)
                    origin_tag = f"mutate(gen{gen})"

                else:  # crossover
                    if len(elite) < 2:
                        continue
                    (_, smi_a, _), (_, smi_b, _) = rng.sample(elite, 2)
                    mol_a = Chem.MolFromSmiles(smi_a)
                    mol_b = Chem.MolFromSmiles(smi_b)
                    if mol_a is None or mol_b is None:
                        continue
                    child = crossover(mol_a, mol_b, rng)
                    origin_tag = f"crossover(gen{gen})"

                if child is None:
                    continue

                child_smi = mol_to_smiles(child)
                if child_smi in best_seen:
                    continue  # já visto

                fit = fitness(child, ref_fps)
                best_seen[child_smi] = fit
                new_individuals.append((fit, child_smi, origin_tag))

            except Exception:
                continue

        population = elite + new_individuals
        population.sort(key=lambda x: x[0], reverse=True)
        population = population[:population_size]

        best_fit = population[0][0] if population else 0.0
        print(f"  Gen {gen:>3}/{generations}  |  tamanho: {len(population)}  |  "
              f"melhor fitness: {best_fit:.4f}", end="\r", flush=True)

    print()
    return population


# ── Build 3D and save ─────────────────────────────────────────────────────────

def save_results(
    candidates: list[tuple[float, str, str]],
    top_n: int,
    out_sdf: Path,
    out_csv: Path,
) -> int:
    """Salva os top_n candidatos em SDF e CSV. Retorna quantos foram salvos."""
    saved = 0
    rows  = []

    writer = SDWriter(str(out_sdf))

    for rank, (fit, smi, origin) in enumerate(candidates[:top_n * 3], 1):
        mol = smiles_to_3d_mol(smi)
        if mol is None:
            continue

        mol.SetProp("_Name",         f"GEN_CAND_{rank:03d}")
        mol.SetProp("SMILES",        smi)
        mol.SetProp("FITNESS",       f"{fit:.4f}")
        mol.SetProp("ORIGIN",        origin)
        mol.SetProp("MW",            f"{Descriptors.MolWt(mol):.2f}")
        mol.SetProp("LOGP",          f"{Descriptors.MolLogP(mol):.2f}")
        mol.SetProp("HBD",           str(Descriptors.NumHDonors(mol)))
        mol.SetProp("HBA",           str(Descriptors.NumHAcceptors(mol)))

        writer.write(mol)
        rows.append({
            "rank":    saved + 1,
            "id":      f"GEN_CAND_{rank:03d}",
            "smiles":  smi,
            "fitness": f"{fit:.4f}",
            "mw":      f"{Descriptors.MolWt(mol):.2f}",
            "logp":    f"{Descriptors.MolLogP(mol):.2f}",
            "hbd":     Descriptors.NumHDonors(mol),
            "hba":     Descriptors.NumHAcceptors(mol),
            "origin":  origin,
        })
        saved += 1
        if saved >= top_n:
            break

    writer.close()

    with open(out_csv, "w", newline="") as f:
        writer_csv = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        writer_csv.writeheader()
        writer_csv.writerows(rows)

    return saved


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="generate_candidate.py — Gera moléculas inibidoras de novo (Etapa 7)"
    )
    parser.add_argument("--generations", type=int, default=30,
                        help="Número de gerações do algoritmo genético (padrão: 30)")
    parser.add_argument("--population",  type=int, default=50,
                        help="Tamanho da população por geração (padrão: 50)")
    parser.add_argument("--top-seed",    type=int, default=20, dest="top_seed",
                        help="Quantos compostos do ranking usar como semente (padrão: 20)")
    parser.add_argument("--top-out",     type=int, default=5,  dest="top_out",
                        help="Quantos candidatos gerados salvar (padrão: 5)")
    parser.add_argument("--seed",        type=int, default=42,
                        help="Semente aleatória para reprodutibilidade (padrão: 42)")
    args = parser.parse_args()

    print(f"\n{CYAN}{'═' * 60}{RESET}")
    print(f"{CYAN}  Geração de Candidatos — Beta-Lactamase Screening{RESET}")
    print(f"{CYAN}{'═' * 60}{RESET}\n")

    rng = random.Random(args.seed)

    # 1. Carregar sementes
    print(f"Carregando top {args.top_seed} compostos do ranking...")
    seeds = load_seed_compounds(args.top_seed)

    if not seeds:
        print(f"{RED}Nenhum composto semente encontrado.{RESET}")
        sys.exit(1)

    # 2. Executar algoritmo genético
    print(f"\nExecutando algoritmo genético:")
    print(f"  Gerações: {args.generations}  |  População: {args.population}")
    candidates = run_genetic_algorithm(
        seeds,
        generations=args.generations,
        population_size=args.population,
        rng=rng,
    )

    if not candidates:
        print(f"{RED}Nenhum candidato gerado.{RESET}")
        sys.exit(1)

    # 3. Salvar resultados
    results_dir = PROJECT_DIR / "results"
    results_dir.mkdir(exist_ok=True)
    out_sdf = results_dir / "generated_candidates.sdf"
    out_csv = results_dir / "generated_candidates.csv"

    print(f"\nGerando coordenadas 3D e salvando top {args.top_out} candidatos...")
    saved = save_results(candidates, args.top_out, out_sdf, out_csv)

    # 4. Resumo
    print(f"\n{GREEN}Concluído!{RESET}")
    print(f"  {saved} candidatos salvos:")
    print(f"    {out_sdf}")
    print(f"    {out_csv}")

    if saved > 0:
        print(f"\nTop candidatos gerados:")
        print(f"  {'#':<4} {'Fitness':<10} {'MW':<8} {'LogP':<7} {'SMILES':<50} Origem")
        print(f"  {'─'*4} {'─'*10} {'─'*8} {'─'*7} {'─'*50} {'─'*20}")
        with open(out_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                smi_short = row["smiles"][:48] + ".." if len(row["smiles"]) > 50 else row["smiles"]
                print(f"  {row['rank']:<4} {row['fitness']:<10} {row['mw']:<8} "
                      f"{row['logp']:<7} {smi_short:<50} {row['origin']}")

    print(f"\nPróxima etapa sugerida:")
    print(f"  Realizar docking dos candidatos gerados com screen_and_rank.py")
    print(f"  Ou visualizar com: {CYAN}python scripts/save_and_visualize.py{RESET}\n")


if __name__ == "__main__":
    main()
