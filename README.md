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

Para escolher entre estruturas já baixadas ou importar um arquivo PDB local:

```bash
python scripts/select_targets.py
```

Para importar diretamente uma beta-lactamase do seu computador:

```bash
python scripts/select_targets.py --local-pdb /caminho/para/minha_beta_lactamase.pdb --force
```

Ao importar um PDB local, o arquivo é copiado para `data/structures/` e passa a ser o único alvo em `config.yaml`.
Isso faz com que a preparação, o docking e o fluxo do algoritmo genético usem apenas essa beta-lactamase.

Ao selecionar estruturas baixadas, a tabela interativa permite escolher os IDs desejados (ex: `1ZG4,1BTL`).

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

### Etapa 4.5 — Análise pós-docking de interação (opcional)

```bash
python scripts/score_interactions.py \
  --target-profile results/experiments/teste/target_profile.json \
  --protein data/structures/betalac13.pdb \
  --ligand-pose results/top50_poses/BETALAC13.3_CYS34-CHARGED__CHEMBL122450.pdbqt \
  --out results/experiments/teste/interaction_score.json \
  --verbose
```

Saída: `results/experiments/teste/interaction_score.json`

---

### Etapa 5 — Ranking multiobjetivo (consensus score)

```bash
python scripts/consensus_score.py \
  --ranking results/ranking.csv \
  --interactions results/experiments/teste/interaction_score.json \
  --out results/experiments/teste/final_candidates.csv \
  --verbose
```

Saída: `results/experiments/teste/final_candidates.csv`

> Para usar um CSV consolidado com múltiplas poses, forneça `--interactions` apontando
> para o arquivo CSV contendo várias linhas com `compound_id` e `interaction_score`.

---

### Etapa 6 — Análise ADMET

```bash
python scripts/admet_analysis.py
```

Saída: `results/admet_report.csv` e `results/final_candidates.csv`

> Esse `results/final_candidates.csv` pertence ao fluxo ADMET antigo. Para a Etapa 9,
> use `results/experiments/<experimento>/final_candidates.csv`, gerado pela Etapa 5
> (`scripts/consensus_score.py`).

---

### Etapa 7 — Salvar complexos e gerar visualizações

```bash
python scripts/save_and_visualize.py
```

Saída:
- `results/complexes/` — arquivos `.pdb` proteína + inibidor
- `results/visualizations/` — visualizações 3D (`.html`) e imagens (`.png`)
- `results/final_report.html` — relatório consolidado com todos os candidatos

---

### Etapa 8 — Geração de candidatos de novo (algoritmo genético)

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

### Etapa 9 — Otimização genética guiada pelo consensus score

```bash
python scripts/genetic_optimize.py \
  --final-candidates results/experiments/teste/final_candidates.csv \
  --compounds data/compounds/compounds.sdf \
  --out results/experiments/teste/generated_candidates.sdf \
  --out-csv results/experiments/teste/generated_candidates_summary.csv \
  --top-n-seeds 10 \
  --generations 10 \
  --population-size 30 \
  --seed 42 \
  --verbose
```

Com a mesma seed, os mesmos inputs e os mesmos parâmetros, a saída deve ser reprodutível.
`--population-size` é o alvo de moléculas únicas após deduplicação; se o algoritmo não conseguir
atingir esse total, ele mostra um aviso.

Filtros químicos mínimos são aplicados antes da saída principal:
`--max-molecular-weight` (650), `--max-tpsa` (250), `--max-hbd` (8),
`--max-hba` (15) e `--min-qed` (0.05).

Saída:
- `generated_candidates.sdf` — apenas moléculas válidas/aprovadas
- `generated_candidates_summary.csv` — apenas moléculas válidas/aprovadas
- `generated_candidates_filtered.csv` — auditoria de moléculas filtradas, criado automaticamente quando `--out-csv` é usado

Os candidatos gerados devem retornar ao ciclo:
- preparação de ligantes
- docking
- interaction scoring
- consensus scoring

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
│   ├── score_interactions.py    # Etapa 4.5
│   ├── consensus_score.py       # Etapa 5
│   ├── admet_analysis.py        # Etapa 6
│   ├── save_and_visualize.py    # Etapa 7
│   ├── generate_candidate.py    # Etapa 8
│   └── genetic_optimize.py      # Etapa 9
└── results/
    ├── ranking.csv
    ├── admet_report.csv
    ├── final_candidates.csv
    ├── experiments/
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
