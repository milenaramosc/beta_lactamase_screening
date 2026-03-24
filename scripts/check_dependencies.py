"""
check_dependencies.py
Verifica se todas as dependências necessárias estão instaladas corretamente.
"""

import sys
import subprocess
import importlib

PYTHON_PACKAGES = [
    ("biopython",                 "Bio"),
    ("requests",                  "requests"),
    ("pyyaml",                    "yaml"),
    ("pandas",                    "pandas"),
    ("rdkit",                     "rdkit"),
    ("meeko",                     "meeko"),
    ("chembl_webresource_client", "chembl_webresource_client"),
    ("py3Dmol",                   "py3Dmol"),
    ("openbabel-wheel",           "openbabel.openbabel"),
]

# (binary, nome, obrigatório)
SYSTEM_TOOLS = [
    ("vina",   "AutoDock Vina", True),   # obrigatório — motor de docking
    ("pymol",  "PyMOL",         False),  # opcional — apenas para imagens PNG estáticas
]

GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
RESET  = "\033[0m"

ok      = f"{GREEN}[OK]{RESET}"
missing = f"{RED}[FALTANDO]{RESET}"
warning = f"{YELLOW}[AVISO]{RESET}"


def check_python_packages():
    print("\n=== Pacotes Python ===")
    all_ok = True
    for pkg_name, import_name in PYTHON_PACKAGES:
        try:
            mod = importlib.import_module(import_name)
            version = getattr(mod, "__version__", "versão desconhecida")
            print(f"  {ok}  {pkg_name} ({version})")
        except ImportError:
            print(f"  {missing}  {pkg_name}  →  pip install {pkg_name}")
            all_ok = False
    return all_ok


def check_system_tools():
    print("\n=== Ferramentas do Sistema ===")
    has_required = True
    for binary, name, required in SYSTEM_TOOLS:
        result = subprocess.run(
            ["which", binary],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            path = result.stdout.strip()
            ver_result = subprocess.run(
                [binary, "--version"],
                capture_output=True,
                text=True
            )
            version_line = (ver_result.stdout or ver_result.stderr).splitlines()
            version = version_line[0] if version_line else "versão desconhecida"
            print(f"  {ok}  {name}: {path}")
            print(f"        {version}")
        else:
            if required:
                print(f"  {missing}  {name} ({binary}) não encontrado no PATH  [OBRIGATÓRIO]")
                has_required = False
            else:
                print(f"  {warning}  {name} ({binary}) não encontrado  [opcional — imagens PNG desativadas]")
    return has_required


def check_vina_path():
    """Verifica também ~/bin/vina que pode não estar no PATH padrão."""
    import os
    local_vina = os.path.expanduser("~/bin/vina")
    if os.path.isfile(local_vina) and os.access(local_vina, os.X_OK):
        result = subprocess.run([local_vina, "--version"], capture_output=True, text=True)
        version = result.stdout.strip() or result.stderr.strip()
        print(f"\n  {warning}  AutoDock Vina encontrado em {local_vina} mas não está no PATH.")
        print(f"        {version}")
        print(f"        Adicione ao PATH com: export PATH=$PATH:~/bin")
        return True
    return False


def check_python_version():
    print("\n=== Versão do Python ===")
    v = sys.version_info
    if v.major >= 3 and v.minor >= 10:
        print(f"  {ok}  Python {v.major}.{v.minor}.{v.micro}")
        return True
    else:
        print(f"  {warning}  Python {v.major}.{v.minor}.{v.micro} — recomendado >= 3.10")
        return False


if __name__ == "__main__":
    print("Verificando dependências do pipeline de triagem virtual...")

    py_ok    = check_python_version()
    pkgs_ok  = check_python_packages()
    tools_ok = check_system_tools()

    # Verifica vina em ~/bin se não foi encontrado no PATH
    if not tools_ok:
        if check_vina_path():
            tools_ok = True

    print("\n=== Resumo ===")
    if py_ok and pkgs_ok and tools_ok:
        print(f"  {ok}  Todas as dependências obrigatórias estão instaladas. Pronto para executar o pipeline.")
        sys.exit(0)
    else:
        if not pkgs_ok:
            print(f"  {missing}  Instale os pacotes Python faltando com:")
            print("         pip install -r requirements.txt")
        if not tools_ok:
            print(f"  {missing}  Instale as ferramentas do sistema faltando (ver README.md).")
        sys.exit(1)
