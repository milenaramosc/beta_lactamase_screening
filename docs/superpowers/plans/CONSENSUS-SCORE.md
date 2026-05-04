Você é um agente de desenvolvimento especializado em Python, docking molecular, bioinformática estrutural, quimioinformática e pipelines científicos reprodutíveis.

IMPORTANTE:
Não use git worktree.
Não crie git worktree.
Não altere a estratégia de branches usando git worktree.
Trabalhe diretamente na estrutura atual do projeto.

Contexto do projeto:
O projeto é um motor em terminal para triagem e geração de candidatos a inibidores de beta-lactamase.

O projeto já possui:

1. Target Profiler:
- scripts/profile_target.py
- src/betalactamase_engine/target_profiler.py
- gera target_profile.json

2. Interaction Scoring:
- scripts/score_interactions.py
- src/betalactamase_engine/scoring/interaction_scoring.py
- gera interaction_score.json ou interaction_score.csv

3. Docking molecular:
- o pipeline já gera um ranking de docking, geralmente em results/ranking.csv
- também pode gerar poses em results/top50_poses/

Objetivo desta tarefa:
Criar a Etapa 3 — ranking multiobjetivo — através de um módulo chamado `consensus_score.py`.

Essa etapa deve combinar:
- energia de docking do AutoDock Vina;
- interaction_score pós-docking;
- pose_classification;
- dados ADMET existentes, se disponíveis;
- penalizações por pose ruim ou dados ausentes;

e gerar um ranking final mais robusto chamado, por exemplo:

results/experiments/teste/final_candidates.csv

Importante:
Não criar interface gráfica.
Não implementar sítios alostéricos.
Não implementar dinâmica molecular.
Não implementar algoritmo genético.
Não refazer o docking.
Não implementar ADMET local agora.
Não depender obrigatoriamente de SwissADME.
Não quebrar o pipeline antigo.
Não remover scripts existentes.

Estrutura desejada:

Criar:

src/
  betalactamase_engine/
    scoring/
      consensus_score.py

Criar também:

scripts/
  consensus_score.py

Se necessário, atualizar:

docs/
  consensus_score.md

E adicionar uma seção curta no README.

Entradas esperadas do script:

Obrigatórias:

--ranking
Caminho para o CSV de ranking do docking.
Exemplo:
results/ranking.csv

--out
Caminho para salvar o ranking final.
Exemplo:
results/experiments/teste/final_candidates.csv

Opcionais:

--interactions
Caminho para arquivo CSV, JSON ou diretório contendo resultados de interaction scoring.
Pode aceitar:
- um CSV com múltiplos compostos;
- um JSON único;
- um diretório com vários arquivos .json gerados pelo score_interactions.py.

--admet
Caminho opcional para CSV de ADMET já existente, caso o projeto tenha gerado.
Exemplo:
results/admet_report.csv
ou
results/final_candidates.csv

--compound-column
Nome da coluna de identificador do composto no ranking de docking.
Se não informado, tentar detectar automaticamente entre:
compound_id, ligand_id, molecule_id, chembl_id, name, title

--vina-column
Nome da coluna de energia/score de docking.
Se não informado, tentar detectar automaticamente entre:
vina_score, binding_affinity, affinity, docking_score, score, best_score

--top-n
Número máximo de candidatos a salvar.
Padrão: salvar todos.

--min-interaction-score
Filtrar candidatos abaixo de determinado interaction_score.
Padrão: não filtrar, apenas penalizar.

--verbose
Exibir resumo no terminal.

Requisitos principais:

1. Ler o ranking de docking.

O arquivo de ranking pode ter nomes de colunas diferentes.
O código deve tentar detectar automaticamente:
- coluna do composto;
- coluna do score Vina.

Se não conseguir detectar, retornar erro amigável explicando quais colunas foram encontradas e como usar:
--compound-column
--vina-column

2. Normalizar o score de docking.

Como no Vina valores mais negativos são melhores, converter para uma escala entre 0.0 e 1.0:

- melhor energia, mais negativa: docking_norm próximo de 1.0
- pior energia: docking_norm próximo de 0.0

Se todos os valores forem iguais, usar docking_norm = 0.5 para todos e registrar warning.

3. Ler interaction scoring.

O script deve aceitar:

a) CSV com colunas:
compound_id,
interaction_score,
pose_classification,
distance_to_site_center,
critical_residue_contact_count,
nearby_residue_contact_count,
metal_contact_count

b) JSON único com estrutura do score_interactions.py

c) Diretório com vários JSONs

Se --interactions não for fornecido:
- continuar funcionando;
- registrar warning;
- calcular ranking apenas com docking e ADMET, se houver;
- não falhar.

4. Integrar interaction_score por compound_id.

Se algum composto do ranking não tiver interaction_score:
- não remover o composto automaticamente;
- preencher interaction_score como vazio ou 0.0;
- aplicar penalização por dado ausente;
- registrar em uma coluna:
missing_interaction_score = true

5. Usar pose_classification como penalização ou bônus.

Classificações esperadas:
- mechanistically_plausible
- weak_active_site_pose
- outside_active_site
- unknown

Sugestão de multiplicadores ou penalizações:

mechanistically_plausible:
- sem penalização;
- pose_bonus = 1.0

weak_active_site_pose:
- penalização leve;
- pose_bonus = 0.7

outside_active_site:
- penalização forte;
- pose_bonus = 0.2

unknown:
- penalização intermediária/forte;
- pose_bonus = 0.4

6. Ler ADMET, se disponível.

Não implementar ADMET local agora.
Apenas consumir CSV existente se o usuário fornecer --admet.

O arquivo ADMET pode ter colunas variadas.
Tentar detectar, se existirem:
- lipinski_violations
- lipinski
- qed
- tpsa
- logp
- molecular_weight
- mw
- hbd
- hba
- rotatable_bonds
- pains
- admet_score

Se houver uma coluna admet_score:
- usar diretamente, normalizando para 0.0 a 1.0 se necessário.

Se não houver admet_score, mas houver propriedades:
- criar um admet_component simples e conservador:
  - penalizar lipinski_violations;
  - penalizar PAINS se presente;
  - beneficiar QED alto, se presente;
  - não falhar se faltarem colunas.

Se --admet não for fornecido:
- continuar funcionando;
- redistribuir pesos;
- registrar warning.

7. Calcular final_score entre 0.0 e 1.0.

Pesos sugeridos quando todos os componentes existem:

- docking_norm: 0.40
- interaction_score: 0.35
- admet_component: 0.20
- pose_bonus: 0.05

Se ADMET não existir:
redistribuir o peso de ADMET proporcionalmente entre docking e interaction.

Exemplo:
- docking_norm: 0.50
- interaction_score: 0.45
- pose_bonus: 0.05

Se interaction_score não existir:
- docking_norm: 0.75
- admet_component: 0.20
- pose_bonus/missing penalty: 0.05

Se nem ADMET nem interaction_score existirem:
- final_score deve refletir principalmente docking_norm;
- registrar warning forte de que o ranking está baseado quase apenas no docking.

8. Aplicar penalizações.

Criar colunas explícitas:
- docking_norm
- interaction_component
- admet_component
- pose_bonus
- missing_interaction_score
- missing_admet
- penalty_notes
- final_score

Penalizações sugeridas:
- outside_active_site: reduzir final_score;
- unknown pose: reduzir final_score;
- interaction_score ausente: penalizar;
- ADMET ausente: não penalizar fortemente se o usuário não forneceu ADMET, mas registrar warning;
- lipinski_violations alto: penalizar;
- PAINS positivo: penalizar, se a coluna existir.

9. Criar classificação final.

Adicionar coluna:

final_classification

Valores possíveis:
- high_priority_candidate
- medium_priority_candidate
- low_priority_candidate
- deprioritized_candidate

Sugestão:
high_priority_candidate:
- final_score >= 0.75
- pose_classification == mechanistically_plausible ou interaction_score >= 0.65

medium_priority_candidate:
- final_score >= 0.55

low_priority_candidate:
- final_score >= 0.35

deprioritized_candidate:
- final_score < 0.35
- ou outside_active_site com interaction_score baixo

10. Não afirmar que um composto é inibidor real.

Usar linguagem conservadora:
- candidato;
- priorizado;
- plausível;
- preliminar;
- compatível com docking;
- requer validação experimental.

Evitar:
- “inibidor confirmado”
- “molécula ativa”
- “bloqueia a enzima”

11. Saída CSV obrigatória.

O arquivo final deve conter pelo menos:

compound_id,
vina_score,
docking_norm,
interaction_score,
pose_classification,
interaction_component,
admet_component,
pose_bonus,
final_score,
final_classification,
missing_interaction_score,
missing_admet,
penalty_notes

Se existirem colunas extras do ranking original, preservar quando possível.

12. Opcionalmente gerar JSON resumo.

Se for simples, além do CSV, gerar um arquivo:

final_candidates_summary.json

com:
- total_candidates;
- candidates_with_interaction_score;
- candidates_with_admet;
- high_priority_count;
- medium_priority_count;
- low_priority_count;
- deprioritized_count;
- weights_used;
- warnings.

Se isso complicar muito, deixar como opcional.

13. Criar CLI:

scripts/consensus_score.py

Exemplo de uso com interaction scoring:

python scripts/consensus_score.py \
  --ranking results/ranking.csv \
  --interactions results/experiments/teste/interaction_scores.csv \
  --out results/experiments/teste/final_candidates.csv \
  --verbose

Exemplo com diretório de JSONs:

python scripts/consensus_score.py \
  --ranking results/ranking.csv \
  --interactions results/experiments/teste/interaction_scores/ \
  --out results/experiments/teste/final_candidates.csv \
  --verbose

Exemplo com ADMET:

python scripts/consensus_score.py \
  --ranking results/ranking.csv \
  --interactions results/experiments/teste/interaction_scores.csv \
  --admet results/admet_report.csv \
  --out results/experiments/teste/final_candidates.csv \
  --verbose

14. Verbose no terminal.

Quando --verbose for usado, imprimir resumo como:

Consensus scoring completed.
Docking ranking: results/ranking.csv
Interactions: results/experiments/teste/interaction_scores.csv
ADMET: not provided
Total candidates: 50
With interaction scores: 42
With ADMET: 0
High priority candidates: 8
Medium priority candidates: 17
Low priority candidates: 20
Deprioritized candidates: 5
Output: results/experiments/teste/final_candidates.csv

15. Tratamento de erros obrigatório.

Lidar de forma amigável com:
- arquivo ranking inexistente;
- ranking CSV vazio;
- coluna de composto não encontrada;
- coluna de Vina não encontrada;
- valores de Vina não numéricos;
- arquivo interaction_score inexistente;
- diretório interactions vazio;
- JSON interaction inválido;
- ADMET CSV inválido;
- compound_id ausente ou duplicado.

16. Duplicatas.

Se houver múltiplos interaction_score para o mesmo compound_id:
- usar o maior interaction_score;
- registrar warning ou penalty_note informando duplicidade;
- preservar a pose_classification associada ao melhor score.

17. Atualizar documentação.

Criar:

docs/consensus_score.md

Explicar:
- objetivo do ranking multiobjetivo;
- por que docking score sozinho não basta;
- como docking_norm é calculado;
- como interaction_score entra no ranking;
- como pose_classification influencia o resultado;
- como ADMET é usado quando disponível;
- como interpretar final_score;
- como interpretar final_classification;
- limitações da abordagem;
- exemplos de comando.

18. Atualizar README.

Adicionar etapa após interaction scoring:

Etapa 3 — Ranking multiobjetivo:

python scripts/consensus_score.py \
  --ranking results/ranking.csv \
  --interactions results/experiments/teste/interaction_scores.csv \
  --out results/experiments/teste/final_candidates.csv \
  --verbose

19. Testes manuais mínimos.

Após implementar, testar pelo menos:

a) Com ranking.csv e sem interactions:
- deve funcionar;
- deve registrar warning;
- deve gerar final_candidates.csv.

b) Com ranking.csv e interactions:
- deve integrar interaction_score;
- deve gerar final_candidates.csv.

c) Se não houver arquivo de interactions com múltiplos compostos, criar um CSV sintético pequeno apenas para validar a integração, deixando claro que é teste técnico sem valor científico.

20. Critérios de aceite.

A tarefa estará concluída quando:

1. O comando abaixo funcionar:

python scripts/consensus_score.py \
  --ranking results/ranking.csv \
  --interactions results/experiments/teste/interaction_scores.csv \
  --out results/experiments/teste/final_candidates.csv \
  --verbose

2. O arquivo final_candidates.csv for gerado.

3. O CSV tiver pelo menos:
- compound_id
- vina_score
- docking_norm
- interaction_score
- pose_classification
- final_score
- final_classification
- missing_interaction_score
- missing_admet
- penalty_notes

4. O script também funcionar sem --interactions, apenas com --ranking.

5. O script não quebrar se --admet não for fornecido.

6. O pipeline antigo continuar funcionando.

7. Não usar git worktree.

8. Não implementar algoritmo genético ainda.

Observação final:
Essa etapa prepara o projeto para a próxima fase: fazer o algoritmo genético usar o ranking multiobjetivo como função de seleção/fitness.
Não implemente essa próxima fase agora.