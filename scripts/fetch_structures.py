"""
fetch_structures.py
Busca e baixa estruturas de beta-lactamases do RCSB PDB por nome/classe.

Uso:
    python scripts/fetch_structures.py --query "beta-lactamase class A" --max 20
    python scripts/fetch_structures.py --query "NDM-1 beta-lactamase" --max 10

Saída:
    data/structures/<ID>.pdb        — arquivos de estrutura
    data/structures_index.json      — índice com metadados de cada estrutura
"""

import argparse
import json
import os
import sys
import time

import requests
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

# ── Constantes ────────────────────────────────────────────────────────────────

RCSB_SEARCH_URL   = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY_URL    = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_URL   = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
RCSB_PDB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"

STRUCTURES_DIR    = os.path.join("data", "structures")
INDEX_FILE        = os.path.join("data", "structures_index.json")

RESOLUTION_WARNING_THRESHOLD = 3.0   # Å

# Cores ANSI
GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
RESET  = "\033[0m"


# ── Busca na API do RCSB ──────────────────────────────────────────────────────

def search_pdb(query: str, max_results: int) -> list[str]:
    """
    Busca estruturas no RCSB PDB por texto livre.
    Retorna lista de PDB IDs ordenados por score de relevância.
    """
    print(f"\n{CYAN}Buscando no RCSB PDB:{RESET} '{query}' (máx: {max_results})")

    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query}
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max_results},
            "sort": [{"sort_by": "score", "direction": "desc"}]
        }
    }

    try:
        resp = requests.post(RCSB_SEARCH_URL, json=payload, timeout=30)
        resp.raise_for_status()
    except requests.exceptions.ConnectionError:
        print(f"{RED}Erro:{RESET} Sem conexão com a internet ou RCSB indisponível.")
        sys.exit(1)
    except requests.exceptions.Timeout:
        print(f"{RED}Erro:{RESET} Tempo de espera esgotado ao contactar o RCSB.")
        sys.exit(1)
    except requests.exceptions.HTTPError as e:
        print(f"{RED}Erro HTTP:{RESET} {e}")
        sys.exit(1)

    data = resp.json()
    print(f"\n {YELLOW}Data: {data}\n")
    print(f"{RESET}")
    total = data.get("total_count", 0)
    results = data.get("result_set", [])
    ids = [r["identifier"] for r in results]

    print(f"  Total encontrado no PDB: {total} estruturas")
    print(f"  Processando os {len(ids)} mais relevantes...")
    return ids


# ── Busca de metadados ────────────────────────────────────────────────────────

def fetch_metadata(pdb_id: str) -> dict:
    """
    Busca metadados de uma estrutura: título, resolução, método e organismo.
    Retorna dicionário com os campos ou valores 'N/A' em caso de ausência.
    """
    metadata = {
        "pdb_id":     pdb_id,
        "title":      "N/A",
        "organism":   "N/A",
        "resolution": None,
        "method":     "N/A",
        "warning":    False,
        "warning_msg": "",
    }

    # Metadados da entrada (título, resolução, método)
    try:
        resp = requests.get(RCSB_ENTRY_URL.format(pdb_id=pdb_id), timeout=15)
        resp.raise_for_status()
        entry = resp.json()

        metadata["title"] = entry.get("struct", {}).get("title", "N/A")

        info = entry.get("rcsb_entry_info", {})
        metadata["method"] = info.get("experimental_method", "N/A")

        resolution_list = info.get("resolution_combined", None)
        if resolution_list:
            metadata["resolution"] = float(resolution_list[0])

    except Exception:
        pass  # metadados parciais são aceitáveis

    # Organismo via entidade polímero
    try:
        resp = requests.get(RCSB_ENTITY_URL.format(pdb_id=pdb_id), timeout=15)
        resp.raise_for_status()
        entity = resp.json()

        src = entity.get("rcsb_entity_source_organism", [])
        if src:
            metadata["organism"] = src[0].get("scientific_name", "N/A")

        # Título mais descritivo a partir da entidade, se disponível
        entity_name = entity.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        if entity_name and metadata["title"] == "N/A":
            metadata["title"] = entity_name

    except Exception:
        pass

    # Aviso de resolução baixa
    if metadata["resolution"] is not None:
        if metadata["resolution"] > RESOLUTION_WARNING_THRESHOLD:
            metadata["warning"] = True
            metadata["warning_msg"] = (
                f"Resolução {metadata['resolution']:.2f} Å > {RESOLUTION_WARNING_THRESHOLD} Å "
                f"— resultados de docking menos confiáveis"
            )

    return metadata


# ── Download do arquivo PDB ───────────────────────────────────────────────────

def download_pdb(pdb_id: str, output_dir: str) -> str | None:
    """
    Baixa o arquivo .pdb de uma estrutura.
    Retorna o caminho do arquivo salvo, ou None em caso de falha.
    """
    filepath = os.path.join(output_dir, f"{pdb_id}.pdb")

    if os.path.isfile(filepath):
        return filepath  # já baixado anteriormente

    url = RCSB_PDB_DOWNLOAD.format(pdb_id=pdb_id)
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        with open(filepath, "w") as f:
            f.write(resp.text)
        return filepath
    except requests.exceptions.HTTPError:
        return None  # estrutura sem arquivo PDB legacy (ex: só mmCIF)
    except Exception:
        return None


# ── Validação com Biopython ───────────────────────────────────────────────────

def validate_pdb(filepath: str) -> bool:
    """
    Verifica se o arquivo .pdb é parseável pelo Biopython.
    Retorna True se válido.
    """
    parser = PDBParser(QUIET=True)
    warnings.filterwarnings("ignore", category=PDBConstructionWarning)
    try:
        structure = parser.get_structure("test", filepath)
        # Verifica se há pelo menos um átomo
        atoms = list(structure.get_atoms())
        return len(atoms) > 0
    except Exception:
        return False


# ── Exibição de resultados ────────────────────────────────────────────────────

def print_summary(index: list[dict]) -> None:
    """Imprime tabela resumo das estruturas baixadas."""
    print(f"\n{'─'*80}")
    print(f"{'ID':<8} {'Título':<35} {'Organismo':<22} {'Res.(Å)':<9} {'Método'}")
    print(f"{'─'*80}")
    for entry in index:
        res_str = f"{entry['resolution']:.2f}" if entry['resolution'] else "N/A"
        warn    = f" {YELLOW}[!]{RESET}" if entry["warning"] else ""
        title   = entry["title"][:34] if entry["title"] != "N/A" else "N/A"
        org     = entry["organism"][:21] if entry["organism"] != "N/A" else "N/A"
        print(f"{entry['pdb_id']:<8} {title:<35} {org:<22} {res_str:<9} {entry['method']}{warn}")
    print(f"{'─'*80}")

    warnings_found = [e for e in index if e["warning"]]
    if warnings_found:
        print(f"\n{YELLOW}[!] Aviso de resolução baixa (> {RESOLUTION_WARNING_THRESHOLD} Å):{RESET}")
        for e in warnings_found:
            print(f"    {e['pdb_id']}: {e['warning_msg']}")


# ── Entrada principal ─────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Busca e baixa estruturas de beta-lactamases do RCSB PDB."
    )
    parser.add_argument(
        "--query", "-q",
        type=str,
        default="beta-lactamase class A",
        help="Texto de busca (padrão: 'beta-lactamase class A')"
    )
    parser.add_argument(
        "--max", "-m",
        type=int,
        default=20,
        dest="max_results",
        help="Número máximo de estruturas a baixar (padrão: 20)"
    )
    args = parser.parse_args()

    os.makedirs(STRUCTURES_DIR, exist_ok=True)
    print(f"{GREEN}Args: {args}\n{RESET}")
    # 1. Busca IDs no PDB
    pdb_ids = search_pdb(args.query, args.max_results)
    print(f"\n{GREEN}pdb_ids: {pdb_ids}\n{RESET}")

    if not pdb_ids:
        print(f"{RED}Nenhuma estrutura encontrada para a query '{args.query}'.{RESET}")
        sys.exit(1)

    # 2. Para cada ID: baixar + buscar metadados
    index = []
    success = 0
    failed  = []

    for i, pdb_id in enumerate(pdb_ids, start=1):
        print(f"  [{i:>3}/{len(pdb_ids)}] {pdb_id} ", end="", flush=True)

        # Download
        filepath = download_pdb(pdb_id, STRUCTURES_DIR)
        if filepath is None:
            print(f"{RED}falhou (sem arquivo PDB disponível){RESET}")
            failed.append(pdb_id)
            time.sleep(0.2)
            continue

        # Validação
        if not validate_pdb(filepath):
            print(f"{YELLOW}inválido (arquivo PDB corrompido){RESET}")
            os.remove(filepath)
            failed.append(pdb_id)
            time.sleep(0.2)
            continue

        # Metadados
        meta = fetch_metadata(pdb_id)
        meta["filepath"] = filepath
        index.append(meta)
        success += 1

        res_str = f"{meta['resolution']:.2f} Å" if meta["resolution"] else "N/A"
        warn    = f" {YELLOW}[resolução baixa]{RESET}" if meta["warning"] else ""
        print(f"{GREEN}OK{RESET}  {res_str}{warn}")

        time.sleep(0.15)  # respeita rate limit da API

    # 3. Salva índice
    with open(INDEX_FILE, "w") as f:
        json.dump(index, f, indent=2, ensure_ascii=False)

    # 4. Resumo final
    print_summary(index)
    print(f"\n{GREEN}Concluído:{RESET} {success} estruturas salvas em '{STRUCTURES_DIR}/'")
    print(f"Índice salvo em '{INDEX_FILE}'")

    if failed:
        print(f"{YELLOW}Não foi possível baixar {len(failed)} estrutura(s):{RESET} {', '.join(failed)}")


if __name__ == "__main__":
    main()
