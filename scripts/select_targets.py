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

import json
import os
import sys

import yaml

# ── Caminhos ──────────────────────────────────────────────────────────────────

INDEX_FILE  = os.path.join("data", "structures_index.json")
CONFIG_FILE = "config.yaml"

# ── Cores ANSI ────────────────────────────────────────────────────────────────

GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"


# ── Leitura do índice ─────────────────────────────────────────────────────────

def load_index() -> list[dict]:
    if not os.path.isfile(INDEX_FILE):
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

def update_config(selected: list[str]) -> None:
    """Atualiza o campo 'selected_targets' no config.yaml."""
    if not os.path.isfile(CONFIG_FILE):
        print(f"{RED}Erro:{RESET} '{CONFIG_FILE}' não encontrado. Execute a partir do diretório raiz do projeto.")
        sys.exit(1)

    with open(CONFIG_FILE) as f:
        config = yaml.safe_load(f)

    config["selected_targets"] = selected

    with open(CONFIG_FILE, "w") as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)

    print(f"{GREEN}Seleção salva em '{CONFIG_FILE}':{RESET} {', '.join(selected)}")


# ── Entrada principal ─────────────────────────────────────────────────────────

def main():
    print(f"\n{CYAN}{'═' * 60}{RESET}")
    print(f"{CYAN}  Seleção de Alvos — Beta-Lactamase Virtual Screening{RESET}")
    print(f"{CYAN}{'═' * 60}{RESET}")

    index    = load_index()
    print(f"\n{YELLOW}index: {index}{RESET}\n")
    print_table(index)

    while True:
        selected = prompt_selection(index)
        if confirm_selection(selected, index):
            break
        print(f"\n{YELLOW}Seleção cancelada. Vamos tentar novamente.{RESET}\n")
        print_table(index)

    update_config(selected)
    print(f"\n{GREEN}Pronto!{RESET} Execute a próxima etapa:")
    print(f"  {CYAN}python scripts/fetch_compounds.py --source chembl --max 5000{RESET}\n")


if __name__ == "__main__":
    main()
