"""
fetch_compounds.py
Baixa biblioteca de compostos candidatos a inibidores de beta-lactamase.

Fontes suportadas:
  - ChEMBL: compostos com atividade biológica documentada (IC50/Ki) contra beta-lactamases
  - ZINC20:  compostos drug-like filtrados por critérios de Lipinski via RDKit
  - both:    combina ChEMBL + ZINC20, desduplicando por SMILES canônico

Uso:
    python scripts/fetch_compounds.py --source chembl --max 5000
    python scripts/fetch_compounds.py --source zinc   --max 5000
    python scripts/fetch_compounds.py --source both   --max 5000

Saída:
    data/compounds/compounds.sdf      — biblioteca de compostos em formato SDF
    data/compounds_summary.json       — metadados da coleta
"""

import argparse
import json
import os
import sys
import time
import warnings

import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, SDWriter
from rdkit import RDLogger


# Evita warnings do RDKit no terminal
RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

# ── Caminhos ──────────────────────────────────────────────────────────────────

COMPOUNDS_DIR  = os.path.join("data", "compounds")
OUTPUT_SDF     = os.path.join(COMPOUNDS_DIR, "compounds.sdf")
SUMMARY_FILE   = os.path.join("data", "compounds_summary.json")

# ── URLs das APIs ─────────────────────────────────────────────────────────────

CHEMBL_TARGET_SEARCH = "https://www.ebi.ac.uk/chembl/api/data/target/search"
CHEMBL_ACTIVITY      = "https://www.ebi.ac.uk/chembl/api/data/activity"
CHEMBL_MOLECULE      = "https://www.ebi.ac.uk/chembl/api/data/molecule"

# Alvos ChEMBL de beta-lactamase pré-selecionados (todas as classes relevantes)
# Gerado via: /target/search?q=beta-lactamase&target_type=SINGLE PROTEIN
CHEMBL_BETALACTAMASE_TARGETS = [
    ("CHEMBL2026",    "Beta-lactamase (TEM-1, Classe A)"),
    ("CHEMBL4114",    "Beta-lactamase (SHV, Classe A)"),
    ("CHEMBL5437",    "Class A carbapenemase KPC-2"),
    ("CHEMBL3562180", "Carbapenem-hydrolyzing beta-lactamase KPC"),
    ("CHEMBL6053",    "Class A carbapenemase"),
    ("CHEMBL1667679", "Beta-lactamase CTX-M-16"),
    ("CHEMBL3499",    "Beta-lactamase class C (AmpC)"),
    ("CHEMBL1255145", "Class D beta-lactamase (OXA)"),
    ("CHEMBL4840",    "Metallo-beta-lactamase type 2 (NDM)"),
    ("CHEMBL5798",    "Beta-lactamase VIM-2 (Classe B)"),
    ("CHEMBL3562178", "Beta-lactamase IMP-1 (Classe B)"),
    ("CHEMBL3326",    "Beta-lactamase L1 (Classe B)"),
    ("CHEMBL2725",    "Beta-lactamase"),
    ("CHEMBL5031",    "Beta-lactamase"),
    ("CHEMBL1641358", "Beta-lactamase"),
]

# ── Cores ANSI ────────────────────────────────────────────────────────────────

GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"


# ── Filtro de Lipinski (Ro5) ──────────────────────────────────────────────────

def passes_lipinski(mol) -> bool:
    """Retorna True se a molécula passa na Regra dos 5 de Lipinski."""
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = Descriptors.NumHDonors(mol)
    hba  = Descriptors.NumHAcceptors(mol)
    return mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    """
    Retorna o maior fragmento de uma molécula (remove sais / contra-íons).
    Se a molécula já é uma única peça, retorna ela mesma.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) == 1:
        return mol
    return max(frags, key=lambda f: f.GetNumHeavyAtoms())


def smiles_to_mol(smiles: str, compound_id: str, name: str = "") -> Chem.Mol | None:
    """
    Converte SMILES para objeto RDKit Mol com coordenadas 3D.
    Sais e contra-íons são removidos; apenas o maior fragmento é mantido.
    Retorna None se inválido.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Remove sais/contra-íons — Meeko exige molécula com 1 único fragmento
    mol = largest_fragment(mol)
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        # fallback: tenta sem ETKDGv3
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if result != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)
    mol.SetProp("_Name", compound_id)
    mol.SetProp("COMPOUND_ID", compound_id)
    mol.SetProp("COMPOUND_NAME", name or compound_id)
    mol.SetProp("SMILES", smiles)
    return mol


# ── Fonte ChEMBL ──────────────────────────────────────────────────────────────

def collect_chembl(max_compounds: int) -> list[dict]:
    """
    Coleta compostos do ChEMBL com atividade (IC50/Ki) contra beta-lactamases.
    Usa lista pré-definida de alvos para evitar chamadas lentas de descoberta.
    O endpoint /activity já retorna canonical_smiles — sem requisição extra por molécula.
    """
    print(f"\n{CYAN}Fonte: ChEMBL{RESET}")
    print(f"  Buscando compostos com IC50/Ki em {len(CHEMBL_BETALACTAMASE_TARGETS)} alvos...")

    seen_ids  = set()
    compounds = []
    limit     = 100

    for tid, tname in CHEMBL_BETALACTAMASE_TARGETS:
        if len(compounds) >= max_compounds:
            break

        offset = 0
        while len(compounds) < max_compounds:
            try:
                resp = requests.get(
                    CHEMBL_ACTIVITY,
                    params={
                        "target_chembl_id":  tid,
                        "standard_type__in": "IC50,Ki",
                        "format":            "json",
                        "limit":             limit,
                        "offset":            offset,
                    },
                    timeout=20
                )
                resp.raise_for_status()
                data = resp.json()
            except Exception as e:
                print(f"\n  {YELLOW}Aviso [{tid}]: {e}{RESET}")
                break

            activities = data.get("activities", [])
            if not activities:
                break

            for act in activities:
                mol_id = act.get("molecule_chembl_id")
                smiles = act.get("canonical_smiles")

                if not mol_id or not smiles or mol_id in seen_ids:
                    continue
                seen_ids.add(mol_id)

                compounds.append({
                    "id":             mol_id,
                    "name":           act.get("molecule_pref_name") or mol_id,
                    "smiles":         smiles,
                    "activity_type":  act.get("standard_type"),
                    "activity_value": act.get("standard_value"),
                    "activity_units": act.get("standard_units"),
                    "target_id":      tid,
                    "target_name":    tname,
                    "source":         "chembl",
                })

                if len(compounds) >= max_compounds:
                    break

            total  = data.get("page_meta", {}).get("total_count", 0)
            offset += limit
            if offset >= total:
                break

            progress = min(len(compounds), max_compounds)
            print(f"  [{progress:>5}/{max_compounds}] compostos  |  alvo: {tname[:40]}...",
                  end="\r", flush=True)
            time.sleep(0.3)

    print()
    print(f"  {len(compounds)} compostos coletados do ChEMBL.")
    return compounds


# ── Fonte ZINC20 ──────────────────────────────────────────────────────────────

def collect_zinc(max_compounds: int) -> list[dict]:
    """
    Coleta compostos drug-like do ChEMBL usando filtros de propriedades
    físico-químicas equivalentes ao filtro drug-like do ZINC20:
      - Peso molecular <= 500 Da
      - AlogP <= 5
      - HBD <= 5, HBA <= 10
      - Apenas small molecules com estrutura definida

    Usa paginação por offset para garantir diversidade real entre chamadas.
    Nota: zinc20.docking.org foi substituído por esta fonte por instabilidade SSL.
    """
    print(f"\n{CYAN}Fonte: ChEMBL (drug-like — equivalente ao ZINC){RESET}")
    print(f"  Buscando moléculas com filtro MW≤500, AlogP≤5, HBD≤5, HBA≤10...")

    compounds = []
    seen_ids  = set()
    limit     = 100
    offset    = 0

    while len(compounds) < max_compounds:
        try:
            resp = requests.get(
                CHEMBL_MOLECULE,
                params={
                    "format":                               "json",
                    "limit":                                limit,
                    "offset":                               offset,
                    "molecule_properties__mw_freebase__lte": 500,
                    "molecule_properties__alogp__lte":       5,
                    "molecule_properties__hbd__lte":         5,
                    "molecule_properties__hba__lte":         10,
                    "molecule_type":                        "Small molecule",
                    "structure_type":                       "MOL",
                },
                timeout=20,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            print(f"\n  {YELLOW}Aviso ao buscar lote ChEMBL drug-like: {e}{RESET}")
            time.sleep(2)
            break

        molecules = data.get("molecules", [])
        if not molecules:
            break

        for mol_data in molecules:
            mol_id = mol_data.get("molecule_chembl_id")
            structs = mol_data.get("molecule_structures") or {}
            smiles  = structs.get("canonical_smiles")

            if not mol_id or not smiles or mol_id in seen_ids:
                continue
            seen_ids.add(mol_id)

            # Validação extra com RDKit + filtro de Lipinski
            mol_check = Chem.MolFromSmiles(smiles)
            if mol_check is None:
                continue
            if not passes_lipinski(mol_check):
                continue

            compounds.append({
                "id":             mol_id,
                "name":           mol_data.get("pref_name") or mol_id,
                "smiles":         smiles,
                "source":         "zinc",   # mantém "zinc" para compatibilidade
                "activity_type":  None,
                "activity_value": None,
                "activity_units": None,
                "target_id":      None,
            })

            if len(compounds) >= max_compounds:
                break

        total = data.get("page_meta", {}).get("total_count", 0)
        offset += limit
        if offset >= total:
            break

        progress = min(len(compounds), max_compounds)
        print(f"  [{progress:>5}/{max_compounds}] compostos coletados...",
              end="\r", flush=True)
        time.sleep(0.3)

    print()
    print(f"  {len(compounds)} compostos drug-like coletados.")
    return compounds


# ── Conversão para SDF ────────────────────────────────────────────────────────

def compounds_to_sdf(compounds: list[dict], output_path: str) -> int:
    """
    Converte lista de compostos (com SMILES) para arquivo SDF com coordenadas 3D.
    Retorna número de moléculas gravadas com sucesso.
    """
    print(f"\n  Gerando coordenadas 3D e salvando em SDF...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    writer  = SDWriter(output_path)
    success = 0
    failed  = 0

    for i, cmpd in enumerate(compounds, start=1):
        mol = smiles_to_mol(cmpd["smiles"], cmpd["id"], cmpd["name"])
        if mol is None:
            failed += 1
            continue

        # Adiciona propriedades extras ao SDF
        if cmpd.get("activity_type"):
            mol.SetProp("ACTIVITY_TYPE",  str(cmpd["activity_type"]))
            mol.SetProp("ACTIVITY_VALUE", str(cmpd["activity_value"] or ""))
            mol.SetProp("ACTIVITY_UNITS", str(cmpd["activity_units"] or ""))
        if cmpd.get("target_id"):
            mol.SetProp("TARGET_ID", str(cmpd["target_id"]))
        mol.SetProp("SOURCE", cmpd["source"])

        writer.write(mol)
        success += 1

        if i % 50 == 0:
            print(f"  [{i:>5}/{len(compounds)}] convertidos...", end="\r", flush=True)

    writer.close()
    print(f"  {success} moléculas gravadas no SDF ({failed} descartadas por erro de geometria).")
    return success


# ── Resumo ────────────────────────────────────────────────────────────────────

def deduplicate_compounds(compounds: list[dict]) -> list[dict]:
    """
    Remove compostos duplicados com base no SMILES canônico do RDKit.
    Em caso de duplicata entre ChEMBL e ZINC, mantém o registro do ChEMBL
    (que carrega dados de atividade biológica).
    """
    seen_smiles: dict[str, dict] = {}  # canonical_smiles → compound

    for cmpd in compounds:
        mol = Chem.MolFromSmiles(cmpd["smiles"])
        if mol is None:
            continue
        canonical = Chem.MolToSmiles(mol)
        if canonical not in seen_smiles:
            seen_smiles[canonical] = cmpd
        else:
            # Mantém ChEMBL sobre ZINC (tem dados de atividade)
            existing = seen_smiles[canonical]
            if existing["source"] == "zinc" and cmpd["source"] == "chembl":
                seen_smiles[canonical] = cmpd

    return list(seen_smiles.values())


def save_summary(source: str, requested: int, collected: int, saved: int,
                 breakdown: dict | None = None) -> None:
    summary = {
        "source":           source,
        "requested":        requested,
        "collected":        collected,
        "saved_to_sdf":     saved,
        "output_sdf":       OUTPUT_SDF,
        "lipinski_filter":  source in ("zinc", "both"),
    }
    if breakdown:
        summary["breakdown"] = breakdown
    with open(SUMMARY_FILE, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Resumo salvo em '{SUMMARY_FILE}'")


# ── Entrada principal ─────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Baixa biblioteca de compostos candidatos a inibidores de beta-lactamase."
    )
    parser.add_argument(
        "--source", "-s",
        choices=["chembl", "zinc", "both"],
        default="chembl",
        help="Fonte dos compostos: 'chembl', 'zinc' ou 'both' (padrão: chembl)"
    )
    parser.add_argument(
        "--max", "-m",
        type=int,
        default=5000,
        dest="max_compounds",
        help="Número máximo de compostos a coletar (padrão: 5000)"
    )
    args = parser.parse_args()

    print(f"\n{CYAN}{'═'*60}{RESET}")
    print(f"{CYAN}  Coleta de Compostos — Beta-Lactamase Virtual Screening{RESET}")
    print(f"{CYAN}{'═'*60}{RESET}")

    os.makedirs(COMPOUNDS_DIR, exist_ok=True)

    # Coleta compostos da(s) fonte(s) escolhida(s)
    if args.source == "chembl":
        compounds = collect_chembl(args.max_compounds)
        breakdown = None
    elif args.source == "zinc":
        compounds = collect_zinc(args.max_compounds)
        breakdown = None
    else:  # both
        half = args.max_compounds // 2
        print(f"\n{CYAN}Modo 'both': coletando até {half} compostos de cada fonte.{RESET}")
        chembl_compounds = collect_chembl(half)
        zinc_compounds   = collect_zinc(half)

        n_chembl = len(chembl_compounds)
        n_zinc   = len(zinc_compounds)
        combined = chembl_compounds + zinc_compounds

        print(f"\n  Desduplicando por SMILES canônico...")
        compounds = deduplicate_compounds(combined)
        removed   = len(combined) - len(compounds)
        print(f"  {len(combined)} compostos combinados → {len(compounds)} únicos "
              f"({removed} duplicatas removidas).")

        breakdown = {
            "chembl": n_chembl,
            "zinc":   n_zinc,
            "duplicates_removed": removed,
        }

    if not compounds:
        print(f"{RED}Nenhum composto coletado. Verifique a conexão e tente novamente.{RESET}")
        sys.exit(1)

    # Converte para SDF com coordenadas 3D
    saved = compounds_to_sdf(compounds, OUTPUT_SDF)

    # Salva resumo
    save_summary(args.source, args.max_compounds, len(compounds), saved,
                 breakdown=breakdown if args.source == "both" else None)

    # Resultado final
    print(f"\n{GREEN}Concluído!{RESET}")
    print(f"  {saved} compostos salvos em '{OUTPUT_SDF}'")
    print(f"\nPróxima etapa:")
    print(f"  {CYAN}python scripts/screen_and_rank.py{RESET}\n")


if __name__ == "__main__":
    main()
