"""
select_targets.py
Apresenta as estruturas de beta-lactamases disponíveis e permite ao usuário
selecionar qual(is) usar como alvo(s) do docking molecular.

Uso:
    python scripts/select_targets.py

Lê:
    data/structures_index.json   — índice gerado pelo fetch_structures.py

Escreve:
    config.yaml                  — campo 'selected_targets' atualizado
"""

import argparse
import json
import re
import shutil
import sys
import warnings
from pathlib import Path

import yaml
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# ── Caminhos ──────────────────────────────────────────────────────────────────

STRUCTURES_DIR = Path("data") / "structures"
INDEX_FILE     = Path("data") / "structures_index.json"
CONFIG_FILE    = Path("config.yaml")

# ── Cores ANSI ────────────────────────────────────────────────────────────────

GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"


# ── Leitura do índice ─────────────────────────────────────────────────────────

def load_index() -> list[dict]:
    if not INDEX_FILE.is_file():
        print(f"{RED}Erro:{RESET} Índice '{INDEX_FILE}' não encontrado.")
        print("Execute primeiro: python scripts/fetch_structures.py")
        sys.exit(1)

    with open(INDEX_FILE) as f:
        index = json.load(f)

    if not index:
        print(f"{RED}Erro:{RESET} Nenhuma estrutura no índice.")
        print("Execute primeiro: python scripts/fetch_structures.py")
        sys.exit(1)

    return index


def load_index_optional() -> list[dict]:
    """Carrega o índice se existir; retorna lista vazia para importação local."""
    if not INDEX_FILE.is_file():
        return []
    with open(INDEX_FILE) as f:
        return json.load(f)


def save_index(index: list[dict]) -> None:
    INDEX_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(INDEX_FILE, "w") as f:
        json.dump(index, f, indent=2, ensure_ascii=False)


# ── Exibição da tabela ────────────────────────────────────────────────────────

def print_table(index: list[dict]) -> None:
    """Imprime tabela formatada com todas as estruturas disponíveis."""
    print(f"\n{BOLD}Estruturas disponíveis:{RESET}\n")
    print(f"{'#':<5} {'ID':<8} {'Título':<38} {'Organismo':<24} {'Res.(Å)':<9} {'Método':<12} {'Aviso'}")
    print("─" * 110)

    for i, entry in enumerate(index, start=1):
        res_str  = f"{entry['resolution']:.2f}" if entry.get("resolution") else "N/A"
        warn_col = f"{YELLOW}[resolução baixa]{RESET}" if entry.get("warning") else ""
        title    = (entry.get("title") or "N/A")[:37]
        org      = (entry.get("organism") or "N/A")[:23]
        method   = (entry.get("method") or "N/A")[:11]

        print(f"{i:<5} {entry['pdb_id']:<8} {title:<38} {org:<24} {res_str:<9} {method:<12} {warn_col}")

    print("─" * 110)


# ── Input e validação ─────────────────────────────────────────────────────────

def prompt_selection(index: list[dict]) -> list[str]:
    """
    Solicita seleção interativa ao usuário.
    Aceita IDs (ex: 1ZG4,1BTL) ou números da tabela (ex: 1,3).
    Retorna lista de PDB IDs válidos selecionados.
    """
    valid_ids     = {e["pdb_id"].upper() for e in index}
    valid_numbers = {str(i): e["pdb_id"] for i, e in enumerate(index, start=1)}

    print(f"\nDigite os {BOLD}IDs{RESET} ou {BOLD}números{RESET} das estruturas desejadas, separados por vírgula.")
    print(f"  Exemplos: {CYAN}1AXB,3BLM{RESET}   ou   {CYAN}1,3{RESET}   ou   {CYAN}1AXB{RESET} (apenas uma)")
    print(f"  Digite {CYAN}all{RESET} para selecionar todas.\n")

    while True:
        raw = input("Seleção: ").strip()

        if not raw:
            print(f"{YELLOW}Nenhuma entrada. Tente novamente.{RESET}")
            continue

        # Selecionar todas
        if raw.lower() == "all":
            selected = [e["pdb_id"] for e in index]
            break

        tokens = [t.strip().upper() for t in raw.split(",") if t.strip()]
        selected = []
        errors   = []

        for token in tokens:
            if token in valid_ids:
                if token not in selected:
                    selected.append(token)
            elif token in valid_numbers:
                pdb_id = valid_numbers[token].upper()
                if pdb_id not in selected:
                    selected.append(pdb_id)
            else:
                errors.append(token)

        if errors:
            print(f"{RED}Não reconhecido:{RESET} {', '.join(errors)}")
            print(f"Use IDs do PDB (ex: 1AXB) ou números da tabela (ex: 1). Tente novamente.\n")
            continue

        if not selected:
            print(f"{YELLOW}Nenhuma estrutura selecionada. Tente novamente.{RESET}\n")
            continue

        break

    return selected


def prompt_source() -> str:
    """Pergunta se o usuário quer usar alvos baixados ou importar um PDB local."""
    print(f"\nEscolha a origem da beta-lactamase:")
    print(f"  {CYAN}1{RESET} — Selecionar estruturas já baixadas em data/structures/")
    print(f"  {CYAN}2{RESET} — Importar arquivo PDB local do computador")

    while True:
        raw = input("Origem [1]: ").strip().lower()
        if raw in ("", "1", "baixadas", "downloaded"):
            return "downloaded"
        if raw in ("2", "local", "pdb"):
            return "local"
        print(f"{YELLOW}Opção inválida. Digite 1 ou 2.{RESET}")


def prompt_local_pdb_path() -> Path:
    while True:
        raw = input("Caminho do arquivo PDB local: ").strip()
        if raw:
            return Path(raw).expanduser()
        print(f"{YELLOW}Informe um caminho para o arquivo .pdb.{RESET}")


def sanitize_target_id(path: Path) -> str:
    target_id = re.sub(r"[^A-Za-z0-9_-]+", "_", path.stem).strip("_")
    return target_id.upper() or "LOCAL_TARGET"


def validate_local_pdb(path: Path) -> None:
    if not path.is_file():
        raise ValueError(f"arquivo não encontrado: {path}")
    if path.suffix.lower() != ".pdb":
        raise ValueError("o arquivo local deve ter extensão .pdb")
    if path.stat().st_size == 0:
        raise ValueError("o arquivo PDB está vazio")

    parser = PDBParser(QUIET=True)
    warnings.filterwarnings("ignore", category=PDBConstructionWarning)
    try:
        structure = parser.get_structure("local_target", str(path))
        if not any(structure.get_atoms()):
            raise ValueError("nenhum átomo encontrado no arquivo PDB")
    except Exception as exc:
        raise ValueError(f"Biopython não conseguiu ler o PDB: {exc}") from exc


def copy_local_pdb(source: Path, target_id: str, force: bool) -> Path:
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    destination = STRUCTURES_DIR / f"{target_id}.pdb"

    try:
        same_file = destination.exists() and source.resolve() == destination.resolve()
    except OSError:
        same_file = False

    if destination.exists() and not same_file and not force:
        raise FileExistsError(
            f"já existe '{destination}'. Use --force ou escolha um arquivo com outro nome."
        )

    if not same_file:
        shutil.copy2(source, destination)

    return destination


def upsert_local_index_entry(index: list[dict], target_id: str, source: Path, destination: Path) -> list[dict]:
    entry = {
        "pdb_id": target_id,
        "title": f"Local PDB: {source.name}",
        "organism": "N/A",
        "resolution": None,
        "method": "local file",
        "warning": False,
        "warning_msg": "",
        "filepath": str(destination),
        "source": "local",
        "source_path": str(source),
    }

    filtered = [e for e in index if e.get("pdb_id", "").upper() != target_id.upper()]
    filtered.append(entry)
    return filtered


def import_local_target(local_pdb: Path, force: bool) -> str:
    validate_local_pdb(local_pdb)
    target_id = sanitize_target_id(local_pdb)
    destination = copy_local_pdb(local_pdb, target_id, force=force)
    index = upsert_local_index_entry(load_index_optional(), target_id, local_pdb, destination)
    save_index(index)
    return target_id


def import_local_target_interactive() -> str:
    while True:
        local_pdb = prompt_local_pdb_path()
        target_id = sanitize_target_id(local_pdb)
        destination = STRUCTURES_DIR / f"{target_id}.pdb"
        force = False

        if destination.exists():
            answer = input(f"'{destination}' já existe. Sobrescrever? [s/N]: ").strip().lower()
            force = answer in ("s", "sim", "y", "yes")
            if not force:
                print(f"{YELLOW}Importação cancelada para esse arquivo. Tente outro caminho.{RESET}\n")
                continue

        try:
            return import_local_target(local_pdb, force=force)
        except (ValueError, FileExistsError) as exc:
            print(f"{RED}Erro:{RESET} {exc}\n")


# ── Confirmação ───────────────────────────────────────────────────────────────

def confirm_selection(selected: list[str], index: list[dict]) -> bool:
    """Exibe as estruturas selecionadas e pede confirmação."""
    index_map = {e["pdb_id"].upper(): e for e in index}

    print(f"\n{BOLD}Estruturas selecionadas:{RESET}")
    for pdb_id in selected:
        entry   = index_map.get(pdb_id, {})
        title   = (entry.get("title") or "N/A")[:50]
        org     = entry.get("organism") or "N/A"
        res_str = f"{entry['resolution']:.2f} Å" if entry.get("resolution") else "N/A"
        warn    = f" {YELLOW}[resolução baixa]{RESET}" if entry.get("warning") else ""
        print(f"  {GREEN}✓{RESET}  {pdb_id}  —  {title}  ({org}, {res_str}){warn}")

    print()
    answer = input("Confirmar seleção? [S/n]: ").strip().lower()
    return answer in ("", "s", "sim", "y", "yes")


# ── Atualização do config.yaml ────────────────────────────────────────────────

def update_config(selected: list[str], clear_binding_sites: bool = False) -> None:
    """Atualiza o campo 'selected_targets' no config.yaml."""
    if not CONFIG_FILE.is_file():
        print(f"{RED}Erro:{RESET} '{CONFIG_FILE}' não encontrado. Execute a partir do diretório raiz do projeto.")
        sys.exit(1)

    with open(CONFIG_FILE) as f:
        config = yaml.safe_load(f)

    config["selected_targets"] = selected
    if clear_binding_sites and isinstance(config.get("binding_site"), dict):
        config["binding_site"]["targets"] = {}

    with open(CONFIG_FILE, "w") as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

    print(f"{GREEN}Seleção salva em '{CONFIG_FILE}':{RESET} {', '.join(selected)}")


def print_next_step(target_id: str | None = None) -> None:
    print(f"\n{GREEN}Pronto!{RESET} Execute a próxima etapa:")
    print(f"  {CYAN}python scripts/prepare_protein.py --force{RESET}")
    print(f"\nPara o GA guiado pela estrutura, use:")
    pdb_path = STRUCTURES_DIR / f"{target_id}.pdb" if target_id else STRUCTURES_DIR / "<ALVO>.pdb"
    print(f"  {CYAN}python scripts/generate_inhibitors_from_protein.py --pdb {pdb_path}{RESET}\n")


# ── Entrada principal ─────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Seleciona alvos de beta-lactamase baixados ou importa um PDB local."
    )
    parser.add_argument(
        "--local-pdb",
        type=Path,
        help="Caminho para um arquivo .pdb local. Ao usar esta opção, só este alvo será selecionado."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Sobrescreve data/structures/<ID>.pdb ao importar um PDB local com o mesmo nome."
    )
    args = parser.parse_args()

    print(f"\n{CYAN}{'═' * 60}{RESET}")
    print(f"{CYAN}  Seleção de Alvos — Beta-Lactamase Virtual Screening{RESET}")
    print(f"{CYAN}{'═' * 60}{RESET}")

    if args.local_pdb:
        try:
            target_id = import_local_target(args.local_pdb.expanduser(), force=args.force)
        except (ValueError, FileExistsError) as exc:
            print(f"{RED}Erro:{RESET} {exc}")
            sys.exit(1)
        update_config([target_id], clear_binding_sites=True)
        print(f"{GREEN}Arquivo local importado:{RESET} {STRUCTURES_DIR / f'{target_id}.pdb'}")
        print_next_step(target_id)
        return

    source = prompt_source()
    if source == "local":
        target_id = import_local_target_interactive()
        update_config([target_id], clear_binding_sites=True)
        print(f"{GREEN}Arquivo local importado:{RESET} {STRUCTURES_DIR / f'{target_id}.pdb'}")
        print_next_step(target_id)
        return

    index = load_index()
    print_table(index)

    while True:
        selected = prompt_selection(index)
        if confirm_selection(selected, index):
            break
        print(f"\n{YELLOW}Seleção cancelada. Vamos tentar novamente.{RESET}\n")
        print_table(index)

    update_config(selected)
    print_next_step()


if __name__ == "__main__":
    main()
