# Pipeline de Triagem Virtual de Inibidores de Beta-Lactamase

Pipeline computacional em Python para identificar candidatos a inibidores de beta-lactamase
usando triagem virtual (virtual screening) com docking molecular.

O objetivo é encontrar compostos que se liguem à beta-lactamase e impeçam sua ação
(como o clavulanato faz em relação à amoxicilina), gerando candidatos para validação
experimental em bancada.

---

## Requisitos do Sistema

- Python >= 3.10
- AutoDock Vina >= 1.2
- OpenBabel >= 3.1
- PyMOL (versão open-source)

---

## Instalação

### 1. Clonar o repositório

```bash
git clone https://github.com/milenaramosc/beta_lactamase_screening
cd beta_lactamase_screening
```

### 2. Instalar dependências Python

```bash
pip install -r requirements.txt
```

### 3. Instalar ferramentas do sistema

**Ubuntu/Debian:**
```bash
sudo apt install openbabel python3-openbabel pymol
```

**Conda (recomendado):**
```bash
conda install -c conda-forge openbabel pymol-open-source autodock-vina
```

**AutoDock Vina (instalação manual):**
```bash
# Baixar binário em https://vina.scripps.edu/downloads/
# Adicionar ao PATH:
export PATH=$PATH:/caminho/para/vina
```

### 4. Verificar instalação

```bash
python scripts/check_dependencies.py
```

---

## Uso

Execute os scripts em ordem, um por vez:

### Etapa 1 — Buscar estruturas de beta-lactamases no PDB

```bash
python scripts/fetch_structures.py --query "beta-lactamase class A" --max 20
```

Saída: `data/structures/*.pdb` e `data/structures_index.json`

---

### Etapa 2 — Selecionar alvos para o docking

```bash
python scripts/select_targets.py
```

Exibe tabela interativa. Selecione os IDs desejados (ex: `1ZG4,1BTL`).
A seleção é salva automaticamente em `config.yaml`.

---

### Etapa 2.5 — Preparar proteínas para o docking

```bash
python scripts/prepare_protein.py
```

Remove moléculas de água e ligantes co-cristalizados (HETATM), adiciona hidrogênios via OpenBabel,
detecta automaticamente o sítio de ligação (registros SITE do PDB) e salva as coordenadas em `config.yaml`.

Saída: `data/prepared/<ID>_prepared.pdb` para cada alvo selecionado.

Use `--force` para reprocessar estruturas já preparadas.

---

### Etapa 3 — Baixar compostos candidatos

```bash
# Usando ChEMBL (recomendado):
python scripts/fetch_compounds.py --source chembl --max 5000

# Usando ZINC20 (redireciona internamente para ChEMBL — ZINC20 inacessível neste ambiente):
python scripts/fetch_compounds.py --source zinc --max 5000

# Usando Ambos
python scripts/fetch_compounds.py --source both --max 5000
```

Saída: `data/compounds/compounds.sdf` e `data/compounds_summary.json`

> **Nota:** `--source zinc` usa a API do ChEMBL internamente com filtros drug-like equivalentes
> (Lipinski: MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10). O ZINC20 está inacessível neste ambiente
> por erro SSL.

---

### Etapa 4 — Docking molecular e ranqueamento

```bash
python scripts/screen_and_rank.py
```

Saída: `results/ranking.csv` e `results/top50_poses/`

> Pode levar várias horas dependendo do número de compostos e CPUs disponíveis.
> O progresso é salvo incrementalmente; interrompa e retome sem perda de dados.

---

### Etapa 5 — Análise ADMET

```bash
python scripts/admet_analysis.py
```

Saída: `results/admet_report.csv` e `results/final_candidates.csv`

---

### Etapa 6 — Salvar complexos e gerar visualizações

```bash
python scripts/save_and_visualize.py
```

Saída:
- `results/complexes/` — arquivos `.pdb` proteína + inibidor
- `results/visualizations/` — visualizações 3D (`.html`) e imagens (`.png`)
- `results/final_report.html` — relatório consolidado com todos os candidatos

---

### Etapa 7 — Geração de candidatos de novo (algoritmo genético)

```bash
python scripts/generate_candidate.py
```

Lê os top compostos de `results/ranking.csv` como população-semente e executa um algoritmo genético
para gerar novas moléculas via fragmentação BRICS (RDKit) + operadores de crossover e mutação.
A função de fitness combina similaridade de Tanimoto com compostos de referência, escore de Lipinski
e penalidade de peso molecular.

Parâmetros opcionais:
```bash
python scripts/generate_candidate.py \
  --generations 20 \
  --population 50 \
  --top-seed 10 \
  --top-out 5 \
  --seed 42
```

Saída: `results/generated_candidates.sdf` e `results/generated_candidates.csv`

---

## Estrutura do Projeto

```
beta_lactamase_screening/
├── config.yaml                  # configuração central do pipeline
├── requirements.txt
├── README.md
├── data/
│   ├── structures/              # estruturas .pdb das beta-lactamases (brutas)
│   ├── prepared/                # estruturas .pdb preparadas (sem HETATM, com H)
│   ├── structures_index.json    # metadados das estruturas baixadas
│   └── compounds/
│       └── compounds.sdf        # biblioteca de compostos candidatos
├── scripts/
│   ├── fetch_structures.py      # Etapa 1
│   ├── select_targets.py        # Etapa 2
│   ├── prepare_protein.py       # Etapa 2.5
│   ├── fetch_compounds.py       # Etapa 3
│   ├── screen_and_rank.py       # Etapa 4
│   ├── admet_analysis.py        # Etapa 5
│   ├── save_and_visualize.py    # Etapa 6
│   └── generate_candidate.py   # Etapa 7
└── results/
    ├── ranking.csv
    ├── admet_report.csv
    ├── final_candidates.csv
    ├── generated_candidates.csv
    ├── generated_candidates.sdf
    ├── top50_poses/
    ├── complexes/
    ├── visualizations/
    └── final_report.html
```

---

## Configuração

Edite `config.yaml` para ajustar parâmetros como:
- `target_query` — termo de busca para as beta-lactamases
- `compound_source` — `"chembl"`, `"zinc"` ou `"both"`
- `top_n_docking` — quantos compostos enviar para análise ADMET
- `vina_exhaustiveness` — precisão do docking (padrão: 8)
- `docking_workers` — número de CPUs para paralelismo

---

## Referências

- RCSB Protein Data Bank: https://www.rcsb.org
- ChEMBL: https://www.ebi.ac.uk/chembl
- ZINC20: https://zinc20.docking.org
- AutoDock Vina: https://vina.scripps.edu
- SwissADME: http://www.swissadme.ch
