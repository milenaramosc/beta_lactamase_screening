# Catalytic Residue Selection Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Permitir ao usuario escolher em tempo de execucao entre residuos cataliticos conhecidos (informados pelo usuario) ou buscar em todo o PDB, removendo o array fixo atual.

**Architecture:** Um prompt interativo coleta a escolha do fallback antes da analise do sitio. A funcao `identify_catalytic_residues` recebe o modo e lista de residuos, mantendo prioridade para SITE records e usando `all` como fallback em execucoes nao interativas.

**Tech Stack:** Python 3, Bio.PDB

---

## Chunk 1: Prompt e controle de fallback

### Task 1: Adicionar prompt e fluxo de fallback

**Files:**
- Modify: `scripts/generate_inhibitors_from_protein.py`

- [ ] **Step 1: Escrever teste manual esperado (nao automatizado)**

Entrada esperada no terminal:

```
Usar residuos cataliticos conhecidos? (s/n): s
Informe residuos (ex.: SER,LYS,GLU): SER,LYS
```

Resultado esperado:
- O fallback usa apenas SER e LYS quando nao ha SITE records.

- [ ] **Step 2: Implementar funcao de prompt**

Adicionar funcao:

```python
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
                    return "known", valid
                print(f"{YELLOW}Entrada invalida. Tente novamente.{RESET}")
            print(f"{YELLOW}Limite de tentativas. Usando fallback all.{RESET}")
            return "all", []
        if decision == "no":
            return "all", []
        print(f"{YELLOW}Entrada invalida. Tente novamente.{RESET}")

    print(f"{YELLOW}Limite de tentativas. Usando fallback all.{RESET}")
    return "all", []
```

- [ ] **Step 2.1: Garantir importacao de sys**

Confirmar que `import sys` existe no topo do arquivo (se nao existir, adicionar).

- [ ] **Step 3: Conectar prompt ao fluxo principal**

Em `main()` antes de `analyze_binding_site`:

```python
fallback_mode, fallback_resnames = prompt_catalytic_selection()
site_info = analyze_binding_site(pdb_path, radius=10.0, fallback_mode=fallback_mode, fallback_resnames=fallback_resnames)
```

- [ ] **Step 4: Atualizar assinatura de `analyze_binding_site`**

```python
def analyze_binding_site(pdb_path: Path, radius: float = 10.0, fallback_mode: str = "all", fallback_resnames: list[str] | None = None) -> dict:
    if fallback_resnames is None:
        fallback_resnames = []
```

E repassar para `identify_catalytic_residues`.

- [ ] **Step 5: Atualizar `identify_catalytic_residues`**

```python
def identify_catalytic_residues(structure, site_residues: list, fallback_mode: str = "all", fallback_resnames: list[str] | None = None) -> list:
    if fallback_resnames is None:
        fallback_resnames = []
```

No fallback:
- `known`: filtrar por `fallback_resnames`.
- `all`: incluir apenas aminoacidos padrao (excluir agua/HOH/WAT e hetero).
- Tratar `site_residues` vazio/invalidado como ausente, aplicando fallback.

- [ ] **Step 6: Ajustar mensagens**

Adicionar mensagens de aviso quando o fallback for ativado e qual modo foi escolhido.

- [ ] **Step 7: Teste manual**

Run:
```
python scripts/generate_inhibitors_from_protein.py --pdb data/structures/1ZG4.pdb
```

Verificar:
- Prompt aparece em TTY.
- Se SITE records existem, eles prevalecem.
- Em entrada invalida, aplica retry e fallback `all`.

Teste adicional (nao interativo):
```
printf "" | python scripts/generate_inhibitors_from_protein.py --pdb data/structures/1ZG4.pdb
```
Verificar aviso e uso de `all`.

---

## Chunk 2: Ajustes de tratamento de residuos

### Task 2: Definir filtro de aminoacidos padrao

**Files:**
- Modify: `scripts/generate_inhibitors_from_protein.py`

- [ ] **Step 1: Adicionar lista de aminoacidos padrao**

```python
STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
}
```

- [ ] **Step 2: Usar filtro no fallback `all`**

```python
if res.get_resname() in STANDARD_RESIDUES:
    for atom in res:
        catalytic_atoms.append(atom)
```

- [ ] **Step 3: Teste manual**

Run:
```
python scripts/generate_inhibitors_from_protein.py --pdb data/structures/1ZG4.pdb
```

Verificar que agua e ligantes hetero nao entram no fallback `all`.
