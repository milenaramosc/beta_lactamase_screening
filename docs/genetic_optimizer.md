# Genetic optimizer (fitness guiado por consensus score)

## Objetivo

Esta etapa usa o `final_score` do `final_candidates.csv` para selecionar sementes e
guiar um algoritmo genetico simples. O resultado sao candidatos preliminares, que
devem retornar ao ciclo de docking e consensus scoring.

Este documento substitui o design inicial de 2026-05-03 para esta etapa.

## Como o final_score e usado

- Seleciona sementes prioritarias (high > medium > low, evitando deprioritized).
- Pondera a selecao de pais na evolucao.
- Mantem o `final_score` como referencia principal do lineage.

## Fitness real vs fitness estimado

- Moleculas sementes herdam o `final_score` como fitness.
- Moleculas geradas recebem um fitness estimado (similaridade com sementes, QED
  e penalizacao por tamanho, TPSA alto, QED baixo e SMILES desconectado).
- O fitness estimado e marcado como preliminar e nao substitui docking.

## Reprodutibilidade

O argumento `--seed` controla as escolhas aleatorias do algoritmo. Com os mesmos
arquivos de entrada, os mesmos parametros e a mesma seed, o SDF e os CSVs
gerados devem ser identicos.

O `population-size` e tratado como alvo do conjunto final apos deduplicacao.
Se o algoritmo nao conseguir recompor moleculas unicas suficientes, ele gera um
warning informando quantas moleculas unicas foram obtidas.

## Entradas

Obrigatorias:

- `--final-candidates`: CSV gerado pelo consensus score.
- `--compounds`: SDF original de compostos.
- `--out`: caminho do SDF de saida.

Opcionais:

- `--out-csv`: resumo em CSV.
- `--top-n-seeds`, `--min-final-score`, `--generations`, `--population-size`.
- `--mutation-rate`, `--crossover-rate`, `--elite-size`.
- `--seed`, `--run-id`, `--verbose`.
- `--max-molecular-weight` (650), `--max-tpsa` (250), `--max-hbd` (8),
  `--max-hba` (15), `--min-qed` (0.05).

Use o CSV gerado por:

```bash
python scripts/consensus_score.py
```

O arquivo esperado e `results/experiments/<experimento>/final_candidates.csv`,
contendo `compound_id`, `final_score` e `final_classification`. Nao use o
`results/final_candidates.csv` antigo da etapa ADMET para esta otimizacao.

## Filtros quimicos minimos

Moleculas geradas fora dos limites configurados nao entram no SDF nem no CSV
principal. Elas sao registradas no CSV de auditoria quando `--out-csv` e usado.

Motivos possiveis em `filter_reason`:

- `mw>{limite}`
- `tpsa>{limite}`
- `hbd>{limite}`
- `hba>{limite}`
- `qed<{limite}`
- `invalid_molecule`

## Saidas

### SDF

Cada molecula inclui:

- `generated_id`
- `parent_1`
- `parent_2`
- `generation`
- `operation` (seed|mutation|crossover|elite|clone)
- `inherited_or_estimated_fitness`
- `source_final_score`
- `source_final_classification`
- `canonical_smiles`

### CSV (opcional)

Colunas minimas:

- `generated_id`
- `canonical_smiles`
- `generation`
- `operation`
- `parent_1`
- `parent_2`
- `inherited_or_estimated_fitness`
- `source_final_score`
- `source_final_classification`
- `molecular_weight`
- `logp`
- `hbd`
- `hba`
- `tpsa`
- `qed`
- `valid_molecule`
- `filter_reason`
- `notes`

`valid_molecule` e sempre `true` no CSV principal, e `filter_reason` fica vazio.

### CSV de auditoria

Quando `--out-csv` e informado, tambem e criado automaticamente:

```text
generated_candidates_filtered.csv
```

Esse arquivo fica no mesmo diretorio do CSV principal e registra moleculas
filtradas ou invalidas. Se nenhuma molecula for filtrada, o arquivo e criado
apenas com cabecalho.

Colunas:

- `generated_id`
- `canonical_smiles`
- `generation`
- `operation`
- `parent_1`
- `parent_2`
- `filter_reason`
- `molecular_weight`
- `tpsa`
- `hbd`
- `hba`
- `qed`
- `notes`
- `valid_molecule`
- `filtered`
- `source_final_score`
- `source_final_classification`

## Limitações

- Fitness estimado e preliminar, nao confirma atividade biologica.
- Nao executa docking nem consenso automaticamente.
- Operadores geneticos sao simples e focados em robustez.

## Exemplo

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

Os candidatos gerados devem retornar ao ciclo de docking e consensus scoring.
