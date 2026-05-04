Você é um agente de desenvolvimento especializado em Python, quimioinformática, docking molecular, algoritmos genéticos, RDKit e pipelines científicos reprodutíveis.

IMPORTANTE:
Não use git worktree.
Não crie git worktree.
Não altere a estratégia de branches usando git worktree.
Trabalhe diretamente na estrutura atual do projeto.

Contexto do projeto:
O projeto é um motor em terminal para triagem, ranqueamento e geração de candidatos a inibidores de beta-lactamase.

O projeto já possui:

1. Target Profiler:
- scripts/profile_target.py
- src/betalactamase_engine/target_profiler.py
- gera target_profile.json

2. Interaction Scoring:
- scripts/score_interactions.py
- src/betalactamase_engine/scoring/interaction_scoring.py
- gera interaction_score.json ou interaction_score.csv

3. Consensus Scoring:
- scripts/consensus_score.py
- src/betalactamase_engine/scoring/consensus_score.py
- gera final_candidates.csv
- final_candidates.csv contém final_score e final_classification

Objetivo desta tarefa:
Implementar a Etapa 4 — fazer o algoritmo genético usar o `final_score` como função de seleção/fitness.

A geração de candidatos de novo não deve mais depender apenas de heurísticas internas, QED, similaridade ou docking isolado. Ela deve usar o ranking multiobjetivo gerado pelo `consensus_score.py`.

O objetivo é que os melhores candidatos do `final_candidates.csv` sejam usados como sementes, e que o algoritmo genético selecione, cruze, mute e priorize novas moléculas com base no `final_score`.

Importante:
Não criar interface gráfica.
Não implementar sítios alostéricos.
Não implementar dinâmica molecular.
Não implementar ADMET local agora.
Não refazer toda a arquitetura do projeto.
Não quebrar o pipeline antigo.
Não remover scripts existentes.
Não afirmar que moléculas geradas são inibidores reais.
Usar termos conservadores como candidato, preliminar, priorizado e requer validação experimental.

Estrutura desejada:

Criar ou adaptar:

src/
  betalactamase_engine/
    generation/
      __init__.py
      genetic_optimizer.py

Criar também:

scripts/
  genetic_optimize.py

Se já existir script semelhante, como:
- scripts/generate_candidate.py
- scripts/generate_inhibitors_from_protein.py

não remover. Pode reaproveitar funções existentes, mas preservar compatibilidade.

Objetivo do novo módulo:
O novo módulo deve receber um `final_candidates.csv` e usar os compostos mais bem classificados como população inicial ou sementes do algoritmo genético.

Entradas obrigatórias do script:

--final-candidates
Caminho para o CSV gerado pelo consensus_score.py.
Exemplo:
results/experiments/teste/final_candidates.csv

--compounds
Caminho para o arquivo SDF original com compostos/sementes.
Exemplo:
data/compounds/compounds.sdf

--out
Caminho para salvar os candidatos gerados.
Exemplo:
results/experiments/teste/generated_candidates.sdf

Argumentos opcionais:

--top-n-seeds
Número de compostos do final_candidates.csv usados como sementes.
Padrão: 10

--generations
Número de gerações do algoritmo genético.
Padrão: 10

--population-size
Tamanho da população por geração.
Padrão: 30

--mutation-rate
Taxa de mutação.
Padrão: 0.25

--crossover-rate
Taxa de crossover.
Padrão: 0.50

--elite-size
Número de melhores candidatos preservados por geração.
Padrão: 5

--min-final-score
Score mínimo para usar um composto como semente.
Padrão: 0.0

--seed
Seed aleatória para reprodutibilidade.
Padrão: 42

--out-csv
Caminho opcional para salvar resumo dos candidatos gerados em CSV.
Exemplo:
results/experiments/teste/generated_candidates_summary.csv

--verbose
Exibir resumo no terminal.

Fluxo desejado:

1. Ler final_candidates.csv.

2. Validar que ele contém pelo menos:
- compound_id
- final_score
- final_classification

3. Ordenar candidatos por final_score decrescente.

4. Selecionar sementes:
- priorizar high_priority_candidate;
- depois medium_priority_candidate;
- depois low_priority_candidate;
- evitar deprioritized_candidate, salvo se não houver candidatos suficientes.

5. Aplicar filtro:
- usar apenas candidatos com final_score >= --min-final-score;
- limitar pelo --top-n-seeds.

6. Ler compounds.sdf.

7. Mapear compound_id do final_candidates.csv para moléculas reais no compounds.sdf.

O código deve tentar casar compound_id com propriedades do SDF, como:
- compound_id
- molecule_id
- chembl_id
- ChEMBL ID
- ID
- name
- _Name

Se não encontrar correspondência:
- registrar warning;
- ignorar essa semente;
- continuar com as demais.

8. Criar população inicial.

A população inicial deve conter:
- moléculas sementes encontradas no SDF;
- variações simples dessas moléculas, se necessário para completar population-size.

9. Implementar fitness inicial baseada em final_score.

Cada molécula semente deve herdar o final_score vindo do final_candidates.csv.

Para novas moléculas geradas, como ainda não passaram por docking e consensus, calcular um fitness provisório usando heurísticas químicas leves, por exemplo:
- similaridade com sementes de alto final_score;
- QED, se RDKit estiver disponível;
- penalização por moléculas inválidas;
- penalização por tamanho molecular extremo;
- penalização por SMILES inválido;
- penalização por moléculas desconectadas ou absurdas.

Importante:
O fitness provisório das novas moléculas deve ser marcado como aproximado/preliminar.

10. Preservar o papel do final_score.

O final_score deve ser a referência principal para seleção das sementes e para ponderar a probabilidade de seleção dos pais.

Exemplo:
- compostos com final_score alto têm maior chance de serem selecionados como pais;
- elite deve preservar os melhores por fitness;
- generated_candidates_summary.csv deve registrar de qual semente ou geração o candidato veio.

11. Operadores genéticos.

Implementar operadores simples e robustos usando RDKit, se já estiver no projeto:

a) Mutação:
- pequenas alterações em fragmentos;
- substituição simples de fragmentos;
- remoção ou adição conservadora de pequenos grupos, se implementável;
- evitar gerar moléculas inválidas.

b) Crossover:
- combinar fragmentos de duas moléculas usando BRICS, se disponível;
- se BRICS falhar, usar fallback simples;
- validar molécula final.

Não precisa criar química perfeita nesta etapa. Priorizar robustez e rastreabilidade.

12. Validação química mínima.

Para cada molécula gerada:
- validar com RDKit;
- sanitizar;
- remover duplicatas por canonical SMILES;
- calcular propriedades básicas:
  - molecular_weight
  - logp
  - hbd
  - hba
  - tpsa
  - qed, se possível

13. Saída SDF.

Salvar moléculas geradas em SDF no caminho de --out.

Cada molécula deve conter propriedades:
- generated_id
- parent_1
- parent_2, se houver
- generation
- operation: seed|mutation|crossover|elite
- inherited_or_estimated_fitness
- source_final_score, quando aplicável
- source_final_classification, quando aplicável
- canonical_smiles

14. Saída CSV.

Se --out-csv for fornecido, salvar CSV com pelo menos:

generated_id,
canonical_smiles,
generation,
operation,
parent_1,
parent_2,
inherited_or_estimated_fitness,
source_final_score,
source_final_classification,
molecular_weight,
logp,
hbd,
hba,
tpsa,
qed,
valid_molecule,
notes

15. Não rodar docking nesta etapa.

Esta etapa deve gerar candidatos novos e ranqueá-los de forma preliminar. O docking dos candidatos gerados será uma etapa posterior.

Adicionar no final do terminal uma orientação como:

"Generated candidates should be submitted to docking and consensus scoring in the next validation cycle."

16. Reprodutibilidade.

O script deve aceitar --seed e usar essa seed para:
- random;
- numpy, se usado;
- qualquer amostragem interna.

17. Tratamento de erros obrigatório.

Lidar de forma amigável com:
- final_candidates.csv inexistente;
- final_candidates.csv vazio;
- ausência de final_score;
- ausência de compound_id;
- compounds.sdf inexistente;
- compounds.sdf sem moléculas válidas;
- nenhuma semente encontrada no SDF;
- RDKit ausente;
- população insuficiente;
- falhas de mutação/crossover.

18. Logs e warnings.

Registrar warnings para:
- compound_id não encontrado no SDF;
- molécula inválida descartada;
- mutação/crossover que falhou;
- população preenchida com sementes repetidas ou variações simples;
- fitness estimado usado em moléculas não redockadas.

19. Verbose.

Quando --verbose for usado, imprimir resumo como:

Genetic optimization completed.
Final candidates input: results/experiments/teste/final_candidates.csv
Seed molecules selected: 10
Seed molecules found in SDF: 8
Generations: 10
Population size: 30
Generated valid molecules: 120
Unique generated molecules: 85
Output SDF: results/experiments/teste/generated_candidates.sdf
Output CSV: results/experiments/teste/generated_candidates_summary.csv

20. Documentação.

Criar:

docs/genetic_optimizer.md

Explicar:
- objetivo do algoritmo genético;
- como o final_score é usado;
- diferença entre fitness real e fitness estimado;
- por que moléculas geradas precisam passar por docking novamente;
- como interpretar generated_candidates.sdf;
- como interpretar generated_candidates_summary.csv;
- limitações da abordagem;
- exemplos de comando.

21. Atualizar README.

Adicionar uma nova etapa após o ranking multiobjetivo:

Etapa 4 — Otimização genética guiada pelo consensus score

Exemplo:

python scripts/genetic_optimize.py \
  --final-candidates results/experiments/teste/final_candidates.csv \
  --compounds data/compounds/compounds.sdf \
  --out results/experiments/teste/generated_candidates.sdf \
  --out-csv results/experiments/teste/generated_candidates_summary.csv \
  --top-n-seeds 10 \
  --generations 10 \
  --population-size 30 \
  --verbose

Explicar que os candidatos gerados devem retornar ao ciclo:
- preparação de ligantes;
- docking;
- interaction scoring;
- consensus scoring.

22. Testes manuais mínimos.

Após implementar, testar:

A) Rodar com final_candidates.csv real:

python scripts/genetic_optimize.py \
  --final-candidates results/experiments/teste/final_candidates.csv \
  --compounds data/compounds/compounds.sdf \
  --out results/experiments/teste/generated_candidates.sdf \
  --out-csv results/experiments/teste/generated_candidates_summary.csv \
  --generations 3 \
  --population-size 10 \
  --verbose

B) Validar que os arquivos foram criados:
- generated_candidates.sdf
- generated_candidates_summary.csv

C) Validar que o CSV contém:
- generated_id
- canonical_smiles
- generation
- operation
- inherited_or_estimated_fitness
- valid_molecule

D) Validar que o SDF contém moléculas válidas.

E) Rodar com --seed duas vezes e conferir que os resultados são reprodutíveis ou pelo menos consistentes.

23. Critérios de aceite.

A tarefa estará concluída quando:

1. O comando abaixo funcionar:

python scripts/genetic_optimize.py \
  --final-candidates results/experiments/teste/final_candidates.csv \
  --compounds data/compounds/compounds.sdf \
  --out results/experiments/teste/generated_candidates.sdf \
  --out-csv results/experiments/teste/generated_candidates_summary.csv \
  --generations 3 \
  --population-size 10 \
  --verbose

2. O algoritmo usar final_score para selecionar sementes e ponderar seleção dos pais.

3. O script gerar SDF com candidatos novos ou variações válidas.

4. O script gerar CSV resumo.

5. O código não rodar docking nesta etapa.

6. O código deixar claro que o fitness de moléculas novas é preliminar até passarem por docking e consensus scoring.

7. O pipeline antigo continuar funcionando.

8. Não usar git worktree.

9. Não implementar alostéricos, dinâmica molecular ou ADMET local.

10. Não remover scripts antigos.

Observação final:
Essa etapa fecha o ciclo de inovação inicial:

Target Profiler
→ Docking
→ Interaction Scoring
→ Consensus Scoring
→ Genetic Optimization guiado por final_score
→ novos candidatos para redocking

Não implemente redocking automático nesta tarefa.