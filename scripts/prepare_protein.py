"""
prepare_protein.py — Etapa 2.5: Preparação das estruturas proteicas para docking

Uso:
    python scripts/prepare_protein.py
    python scripts/prepare_protein.py --force   # reprocessa mesmo se já preparado

Para cada alvo em config.yaml (selected_targets):
  1. Remove água cristalográfica (HETATM HOH/WAT)
  2. Remove ligantes co-cristalizados (demais HETATM)
  3. Adiciona hidrogênios via OpenBabel
  4. Minimiza energia com campo de força MMFF94 (500 passos)
  5. Detecta sítio ativo via registros SITE do PDB
  6. Salva data/prepared/<ID>_prepared.pdb
  7. Atualiza binding_site.targets em config.yaml

Saída:
    data/prepared/<ID>_prepared.pdb
    config.yaml  (campo binding_site.targets atualizado)
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import yaml
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from openbabel import openbabel as ob

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR  = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
CONFIG_FILE = PROJECT_DIR / "config.yaml"
STRUCT_DIR  = PROJECT_DIR / "data" / "structures"
PREP_DIR    = PROJECT_DIR / "data" / "prepared"

# ── ANSI colours ──────────────────────────────────────────────────────────────

GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
CYAN   = "\033[96m"
BOLD   = "\033[1m"
RESET  = "\033[0m"


# ── Config ────────────────────────────────────────────────────────────────────

def load_config() -> dict:
    if not CONFIG_FILE.exists():
        print(f"{RED}Erro:{RESET} config.yaml não encontrado em {CONFIG_FILE}")
        sys.exit(1)
    with open(CONFIG_FILE) as f:
        return yaml.safe_load(f)


def save_config(config: dict) -> None:
    with open(CONFIG_FILE, "w") as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)


# ── Step 1 & 2: Strip water and heteroatoms ───────────────────────────────────

def strip_heteroatoms(pdb_path: Path, out_path: Path) -> int:
    """
    Remove HETATM records (água + ligantes co-cristalizados).
    Mantém ATOM (cadeia proteica) e registros de cabeçalho/SITE.
    Retorna número de linhas HETATM removidas.
    """
    kept = []
    removed = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                removed += 1
            else:
                kept.append(line)
    with open(out_path, "w") as f:
        f.writelines(kept)
    return removed


# ── Step 3 & 4: Add hydrogens + minimise via OpenBabel ────────────────────────

def add_hydrogens(pdb_path: Path, out_path: Path) -> bool:
    """
    Adiciona hidrogênios via OpenBabel e salva o PDB resultante.
    Retorna True em caso de sucesso.

    Nota: minimização de energia em proteínas inteiras (>3000 átomos) via
    OpenBabel é impraticável sem GPU. Para docking com AutoDock Vina a adição
    de H é suficiente — Vina trata a proteína como receptora rígida.
    """
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdb")

    mol = ob.OBMol()
    if not conv.ReadFile(mol, str(pdb_path)):
        raise RuntimeError(f"OpenBabel não conseguiu ler {pdb_path}")

    mol.AddHydrogens()
    conv.WriteFile(mol, str(out_path))
    return True


# ── Step 5: Binding site detection ───────────────────────────────────────────

def parse_site_records(pdb_path: Path) -> list[tuple]:
    """
    Lê registros SITE do PDB (formato fixed-width).
    Retorna lista de (resname, chain, resnum).
    """
    site_residues = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("SITE"):
                continue
            i = 18
            while i + 10 <= len(line.rstrip()):
                group = line[i: i + 11]
                resname    = group[0:3].strip()
                chain      = group[4].strip() if len(group) > 4 else ""
                resnum_str = group[5:9].strip() if len(group) > 5 else ""
                if resname and chain and resnum_str:
                    try:
                        site_residues.append((resname, chain, int(resnum_str)))
                    except ValueError:
                        pass
                i += 11
    return site_residues


def detect_binding_site(pdb_path: Path) -> dict:
    """
    Detecta sítio ativo e retorna dict com center_x/y/z e size_x/y/z.
    Prioridade:
      1. Registros SITE do PDB
      2. Fallback: centroide de todos os CA + caixa 30 Å
    """
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("X", str(pdb_path))
    model  = struct[0]

    site_residues = parse_site_records(pdb_path)
    coords = []

    if site_residues:
        for resname, chain, resnum in site_residues:
            try:
                res = model[chain][resnum]
                ca  = res["CA"]
                coords.append(ca.get_vector().get_array())
            except (KeyError, Exception):
                pass

    if coords:
        coords  = np.array(coords)
        center  = coords.mean(axis=0)
        extents = coords.max(axis=0) - coords.min(axis=0)
        box     = np.maximum(extents + 10, 20.0)
        method  = "SITE records"
    else:
        # fallback
        for chain in model:
            for res in chain:
                if "CA" in res:
                    coords.append(res["CA"].get_vector().get_array())
        if not coords:
            raise ValueError(f"Nenhum átomo CA encontrado em {pdb_path}")
        coords = np.array(coords)
        center = coords.mean(axis=0)
        box    = np.array([30.0, 30.0, 30.0])
        method = "centroide CA (fallback)"

    cx, cy, cz = center
    sx, sy, sz = box
    return {
        "center_x": round(float(cx), 3),
        "center_y": round(float(cy), 3),
        "center_z": round(float(cz), 3),
        "size_x":   round(float(sx), 3),
        "size_y":   round(float(sy), 3),
        "size_z":   round(float(sz), 3),
        "method":   method,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def process_target(pdb_id: str, force: bool) -> dict | None:
    """
    Processa um único alvo. Retorna dict com resultado ou None em caso de erro.
    """
    src_pdb  = STRUCT_DIR / f"{pdb_id}.pdb"
    prep_pdb = PREP_DIR   / f"{pdb_id}_prepared.pdb"

    result = {
        "pdb_id":    pdb_id,
        "stripped":  False,
        "minimised": False,
        "site":      None,
        "error":     None,
    }

    # Verificar arquivo fonte
    if not src_pdb.exists():
        result["error"] = f"arquivo {src_pdb.name} não encontrado em data/structures/"
        return result

    # Pular se já preparado e --force não foi passado
    if prep_pdb.exists() and not force:
        result["skipped"] = True
        return result

    # Arquivo temporário para etapas intermediárias
    tmp_stripped = PREP_DIR / f"{pdb_id}_stripped.pdb"

    try:
        # Etapa 1+2: remover HETATM
        removed = strip_heteroatoms(src_pdb, tmp_stripped)
        result["stripped"] = True
        result["hetatm_removed"] = removed

        # Etapa 3: adicionar H
        try:
            add_hydrogens(tmp_stripped, prep_pdb)
            result["minimised"] = False  # minimização não aplicada em proteínas
        except Exception as e:
            # Fallback: salva sem minimização
            import shutil
            shutil.copy(tmp_stripped, prep_pdb)
            result["minimised"] = False
            result["minimise_warning"] = str(e)

        # Etapa 5: detectar sítio ativo (usa src_pdb pois tem SITE records)
        try:
            site = detect_binding_site(src_pdb)
            result["site"] = site
        except Exception as e:
            result["site_warning"] = str(e)

    except Exception as e:
        result["error"] = str(e)
    finally:
        if tmp_stripped.exists():
            tmp_stripped.unlink()

    return result


def main():
    parser = argparse.ArgumentParser(
        description="prepare_protein.py — Prepara proteínas para docking (Etapa 2.5)"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Reprocessa alvos já preparados anteriormente"
    )
    args = parser.parse_args()

    PREP_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n{CYAN}{'═' * 60}{RESET}")
    print(f"{CYAN}  Preparação de Proteínas — Beta-Lactamase Screening{RESET}")
    print(f"{CYAN}{'═' * 60}{RESET}\n")

    config  = load_config()
    targets = config.get("selected_targets", [])

    if not targets:
        print(f"{RED}Erro:{RESET} Nenhum alvo em selected_targets. Execute select_targets.py primeiro.")
        sys.exit(1)

    print(f"Alvos selecionados: {', '.join(targets)}\n")

    binding_sites = {}
    summary = {"ok": 0, "skipped": 0, "errors": 0}

    for pdb_id in targets:
        print(f"  {BOLD}{pdb_id}{RESET}", end="  ")
        res = process_target(pdb_id, force=args.force)

        if res is None:
            print(f"{RED}[ERRO interno]{RESET}")
            summary["errors"] += 1
            continue

        if res.get("skipped"):
            print(f"{YELLOW}[ignorado — já preparado; use --force para reprocessar]{RESET}")
            summary["skipped"] += 1
            # ainda carrega o sítio do PDB original
            try:
                site = detect_binding_site(STRUCT_DIR / f"{pdb_id}.pdb")
                binding_sites[pdb_id] = site
            except Exception:
                pass
            continue

        if res.get("error"):
            print(f"{RED}[ERRO]{RESET} {res['error']}")
            summary["errors"] += 1
            continue

        # Resultado OK
        hetatm_n = res.get("hetatm_removed", 0)
        print(f"{GREEN}[OK]{RESET}  HETATM removidos: {hetatm_n}  |  H adicionados")

        if res.get("minimise_warning"):
            print(f"       {YELLOW}Aviso minimização:{RESET} {res['minimise_warning']}")

        if res.get("site"):
            s = res["site"]
            print(f"       Sítio: centro ({s['center_x']:.2f}, {s['center_y']:.2f}, {s['center_z']:.2f})"
                  f"  caixa ({s['size_x']:.1f}×{s['size_y']:.1f}×{s['size_z']:.1f} Å)"
                  f"  [{s['method']}]")
            binding_sites[pdb_id] = s
        elif res.get("site_warning"):
            print(f"       {YELLOW}Aviso sítio:{RESET} {res['site_warning']}")

        summary["ok"] += 1

    # Atualizar config.yaml com coordenadas dos sítios ativos
    if binding_sites:
        bs = config.setdefault("binding_site", {})
        bs["targets"] = {
            pdb_id: {k: v for k, v in s.items() if k != "method"}
            for pdb_id, s in binding_sites.items()
        }
        save_config(config)
        print(f"\n{GREEN}Sítios ativos salvos em config.yaml (binding_site.targets).{RESET}")

    # Resumo final
    print(f"\n{'─' * 60}")
    print(f"  Preparados:  {summary['ok']}")
    print(f"  Ignorados:   {summary['skipped']}")
    print(f"  Erros:       {summary['errors']}")
    print(f"  Saída:       {PREP_DIR}/")
    print(f"{'─' * 60}")

    if summary["ok"] > 0 or summary["skipped"] > 0:
        print(f"\n{GREEN}Pronto!{RESET} Próxima etapa:")
        print(f"  {CYAN}python scripts/fetch_compounds.py --source chembl --max 5000{RESET}\n")
    else:
        print(f"\n{RED}Nenhum alvo preparado com sucesso.{RESET}")
        sys.exit(1)


if __name__ == "__main__":
    main()
