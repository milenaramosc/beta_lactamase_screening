# Genetic Optimizer Hardening (final_score-guided)

## Contexto
Precisamos finalizar a etapa de otimizacao genetica guiada por final_score.
No README essa etapa aparece como Etapa 9 e deve permanecer nesse ponto do
pipeline (apos Etapa 8). Quando mencionar o ciclo interno de otimizacao, usar
a expressao "etapa de otimizacao genetica" sem renumerar.
determinismo forte, filtros quimicos minimos, melhor diversidade, erro amigavel
no CLI e saidas limpas para o proximo ciclo de docking. Nao usar git worktree e
nao alterar scripts antigos.

## Objetivos
- Reprodutibilidade: mesmo comando com `--seed` gera o mesmo SDF/CSV.
- SDF limpo: somente moleculas validas e aprovadas em filtros.
- CSV principal limpo: somente moleculas validas/aprovadas.
- CSV de auditoria: registrar filtradas quando `--out-csv` for informado.
- Diversidade: tentar recompor moleculas unicas apos dedup final.
- Fitness estimado mais conservador para tamanhos/propiedades extremas.
- CLI com erro amigavel para CSV antigo sem `final_score`/`final_classification`.
- Documentacao clara: seed, filtros, alvo de population-size, uso do final_candidates correto.
- Limpeza do projeto: remover __pycache__ e *.pyc; ajustar .gitignore.

## Nao objetivos
- Nao rodar docking automatico.
- Nao implementar ADMET local.
- Nao implementar sitios alostericos.
- Nao implementar dinamica molecular.
- Nao remover scripts antigos. Nao modificar scripts existentes (exceto
  `scripts/genetic_optimize.py`, que e o script desta etapa).
  Alteracoes permitidas apenas em:
  - `src/betalactamase_engine/generation/genetic_optimizer.py`
  - `scripts/genetic_optimize.py`
  - `docs/genetic_optimizer.md`
  - `README.md`
  - `.gitignore`

## Abordagem (determinismo)
- Ordenar fragmentos BRICS com `sorted(BRICS.BRICSDecompose(mol))`.
- Ordenar listas derivadas de dict/set antes de amostragem.
- Ordenar populacao final antes de salvar por:
  - `inherited_or_estimated_fitness` desc
  - `canonical_smiles` asc
- IDs gerados baseados na ordenacao final.
- Selecionar parent_2 diferente de parent_1 quando possivel; caso contrario,
  permitir e registrar nota `crossover_same_parent`.

Determinismo:
- Fixar `random.seed(seed)` e, se numpy vier a ser usado, `numpy.random.seed(seed)`.
- RDKit nao recebe seed global aqui; o determinismo e garantido por:
  - ordenacao de fragmentos BRICS
  - ordenacao da populacao antes de salvar
  - ausencia de sets/dicts sem ordenacao antes de amostragem
  - ordenacao deterministica de `filter_reason`
  - inputs de mutacao/crossover ordenados e escolhas via RNG com seed fixa
  - o SDF de entrada (`--compounds`) sera ordenado de forma deterministica por
    um indice de propriedades (SDF_KEYS) antes do mapeamento; isso reduz
    variacao por ordem.
  - SDF_KEYS: ["compound_id", "chembl_id", "molecule_id", "ChEMBL ID", "ID", "name", "_Name"]
  - chave de ordenacao: tupla com valores de todas as chaves acima (string
    vazia quando ausente; valores convertidos para string e strip) + SMILES
    canonico como desempate. Moleculas invalidas sao descartadas antes da
    ordenacao.

## Diversidade apos deduplicacao
- Apos a deduplicacao final, tentar repor moleculas unicas via mutacao/crossover
  ate `population_size` ou `max_attempts`.
- Se nao atingir, registrar warning:
  `Only N unique molecules could be generated after deduplication.`

Escopo de `population_size`:
- O alvo principal e o conjunto final apos deduplicacao. As geracoes internas
  mantem `population_size`, mas o criterio de diversidade se aplica ao output.

Regras de `population_size`:
- Se `population_size` <= 1, nao tentar repor diversidade apos deduplicacao.

`max_attempts`:
- Valor default: `population_size * 8` (interno, sem novo flag).
- Conta tentativas de geracao apos a deduplicacao final (cada tentativa gera
  uma nova molecula candidata). Tentativas contam mesmo se a molecula falhar
  na sanitizacao ou for filtrada.

## Filtros quimicos minimos
Novos argumentos CLI (com defaults):
- `--max-molecular-weight` (650)
- `--max-tpsa` (250)
- `--max-hbd` (8)
- `--max-hba` (15)
- `--min-qed` (0.05)

Regras:
- Moleculas geradas fora dos limites sao filtradas.
- Seeds podem ser preservadas para evolucao, mas sao filtradas na saida se
  estiverem fora dos limites. Elas entram apenas no CSV de auditoria.
- `filter_reason` deve ser lista separada por ponto e virgula.
- `filter_reason` deve refletir os limites configurados (ex.: se
  `--max-molecular-weight=700`, usar `mw>700`).

Valores possiveis de `filter_reason`:
- `mw>{max_molecular_weight}`
- `tpsa>{max_tpsa}`
- `hbd>{max_hbd}`
- `hba>{max_hba}`
- `qed<{min_qed}`
- `invalid_molecule`

Ordem de `filter_reason`:
- Seguir esta ordem fixa para reprodutibilidade: mw, tpsa, hbd, hba, qed.

## Saidas
- `generated_candidates.sdf`: somente validas/aprovadas (nome definido por
  `--out`).
- `generated_candidates_summary.csv`: somente validas/aprovadas (nome definido
  por `--out-csv`) e com `valid_molecule=true` (coluna mantida para
  compatibilidade).
- `generated_candidates_filtered.csv`: gerado automaticamente quando `--out-csv`
  for informado, no mesmo diretorio do `--out-csv`, com nome fixo
  `generated_candidates_filtered.csv`.
  Campos:
  `generated_id, canonical_smiles, generation, operation, parent_1, parent_2,
  filter_reason, molecular_weight, tpsa, hbd, hba, qed, notes, valid_molecule,
  filtered, source_final_score, source_final_classification`.
  Se nao houver filtradas, gerar o arquivo apenas com cabecalho e registrar
  no verbose: "No filtered molecules were generated."

Observacao:
- `generated_candidates_filtered.csv` so e criado quando `--out-csv` e fornecido.

Campos no CSV de auditoria:
- `valid_molecule=false`
- `filtered=true`
- `filter_reason` conforme lista acima

Moleculas invalidas:
- Falhas de sanitizacao nao entram no SDF nem no CSV principal.
- Registrar no CSV de auditoria com `filter_reason=invalid_molecule`.

Campos em linhas validas:
- `filter_reason`: string vazia (""), tanto no CSV quanto no SDF
- `filtered`: coluna omitida no CSV principal (somente no CSV de auditoria)
 - `filter_reason` deve ser propriedade presente no SDF para validas (vazia).

Mapeamento de colunas:
- `final_candidates.csv.final_score` -> `source_final_score`
- `final_candidates.csv.final_classification` -> `source_final_classification`

CSV principal:
- manter `logp` como coluna do resumo (compatibilidade com doc existente),
  calculado via `Descriptors.MolLogP`.

`crossover_same_parent`:
- Registrar em `notes` quando parent_1 == parent_2 em crossover por falta de
  alternativas.

## Fitness estimado
- Penalizar mais fortemente MW > 500 e TPSA > 180.
- Penalizar QED baixo.
- Penalizar SMILES desconectado (`.` no SMILES canonico).
- Penalizar numero de fragmentos > 3 (usando `Chem.GetMolFrags`).
  Usar `Chem.GetMolFrags(mol, asMols=False)` e contar fragmentos desconectados.

Pesos propostos (ajuste direto no modulo):
- `score = 0.60*similaridade + 0.20*qed - 0.12*penalidade_mw - 0.06*penalidade_tpsa`
- penalidade_mw = max(0, (mw-500)/200)
- penalidade_tpsa = max(0, (tpsa-180)/80)
- bonus/penalidade adicional: -0.10 se SMILES desconectado, -0.05 se fragments>3

Similaridade:
- usar fingerprints das seeds mapeadas (Morgan, radius=2, 2048 bits), com
  Tanimoto maximo como referencia.

Regra de limite:
- Clampar `inherited_or_estimated_fitness` para [0.0, 1.0] apos penalidades.

## Erros amigaveis no CLI
- Capturar `ValueError` e `FileNotFoundError` e imprimir:
  "Invalid final candidates file. This step requires the CSV generated by
  scripts/consensus_score.py containing compound_id, final_score and
  final_classification."

Detectar colunas ausentes com validacao explicita e gerar `ValueError`.

`compound_id`:
- Obrigatorio para mapear sementes ao SDF e manter lineage (parent_1/parent_2).

## Documentacao
- Atualizar `docs/genetic_optimizer.md` e `README.md` com:
  - reprodutibilidade e significado de `--seed`
  - filtros quimicos minimos
  - `population-size` como alvo
  - uso de `results/experiments/<experimento>/final_candidates.csv`
  - nao usar `results/final_candidates.csv` antigo (ADMET)
- ajustar a secao ADMET para evitar conflito com consensus
- CSV de auditoria gerado automaticamente
- incluir esquema do `generated_candidates_filtered.csv` em docs
- corrigir inconsistencias de numeracao de etapas (Etapa 9 deve ser unica)
- alinhar numeracao tanto no fluxo de uso quanto na secao "Estrutura do Projeto"
  - numeracao proposta do pipeline:
    - Etapa 1: fetch_structures.py
    - Etapa 2: select_targets.py
    - Etapa 2.5: prepare_protein.py
    - Etapa 3: fetch_compounds.py
    - Etapa 4: screen_and_rank.py
    - Etapa 4.5: score_interactions.py
    - Etapa 5: consensus_score.py
    - Etapa 6: admet_analysis.py
    - Etapa 7: save_and_visualize.py
    - Etapa 8: generate_candidate.py
    - Etapa 9: genetic_optimize.py
- manter estilo de linguagem consistente com README (portugues claro, sem
  trocar para ingles)
- substituir frases em ingles (ex.: "Generated candidates...") por portugues
- diferenciar Etapa 8 (generate_candidate.py) e Etapa 9 (genetic_optimize.py)
  e manter ambas no README
- adicionar `consensus_score.py` e `genetic_optimize.py` na "Estrutura do Projeto"
 - registrar que `generated_candidates_filtered.csv` so existe quando `--out-csv` e usado
  - registrar que este documento substitui o design 2026-05-03 para esta etapa

## Limpeza do projeto
- Remover todos `__pycache__/` e arquivos `*.pyc`.
- Manter a convencao atual do repo: continuar ignorando `results/*` e
  `data/structures/*` no `.gitignore`.
- Ajustar `.gitignore` apenas para remover entradas redundantes (ex.: duplicatas
  de `results/` e `results/*`) e manter as regras existentes para resultados
  e estruturas.

## Decisoes obrigatorias
- Nenhuma decisao pendente de `.gitignore` (manter ignorados `results/*` e
  `data/structures/*`).

## Notas de implementacao
- Regras de processo (nao usar git worktree, nao alterar scripts antigos) sao
  restricoes do ticket e devem ser seguidas na implementacao.

IDs gerados:
- Formato `GEN_<run_id>_<contador>` com contador apos ordenacao final para
  garantir estabilidade entre execucoes com a mesma seed.
- `run_id`: se `--run-id` for informado, usar esse valor; caso contrario,
  derivar de `seed` via hash (sha256 do seed, 8 chars) como ja implementado.

Ordenacao final:
- `inherited_or_estimated_fitness` e o campo interno usado para ordenar a
  populacao antes de gerar IDs e escrever SDF/CSV.
- SMILES canonico usado no desempate deve ser calculado apos sanitizacao.

## Testes manuais
- Rodar duas vezes com `--seed 42` e confirmar SDF/CSV identicos.
- Verificar que o SDF/CSV principais nao incluem filtradas.
- Confirmar warnings quando diversidade insuficiente.
- Validar erro amigavel com CSV antigo sem `final_score`.
Definicao de "valida/aprovada":
- RDKit sanitiza com sucesso.
- Atende aos filtros quimicos minimos (mw/tpsa/hbd/hba/qed).
- Nao e duplicata apos deduplicacao final.
