Você é um agente de desenvolvimento especializado em Python, bioinformática estrutural, docking molecular, quimioinformática e pipelines científicos reprodutíveis.

IMPORTANTE:
Não use git worktree.
Não crie git worktree.
Não altere a estratégia de branches usando git worktree.
Trabalhe diretamente na estrutura atual do projeto.

Contexto do projeto:
O projeto é um motor em terminal para triagem e geração de candidatos a inibidores de beta-lactamase.

O projeto já possui uma etapa chamada Target Profiler, responsável por analisar arquivos PDB e gerar um arquivo `target_profile.json`.

Esse JSON contém informações como:
- target;
- selected_chain;
- classification;
- active_site;
- metals;
- heteroatoms;
- pocket_properties;
- quality_control;
- warnings.

O objetivo agora é implementar a Etapa 2: análise pós-docking da pose do ligante.

Objetivo da tarefa:
Criar um módulo chamado `interaction_scoring.py` para avaliar se a pose gerada pelo docking faz sentido biologicamente, em vez de depender apenas da energia do AutoDock Vina.

Essa etapa deve responder perguntas como:
- O ligante ficou próximo ao centro do sítio ativo?
- O ligante interage com resíduos catalíticos candidatos?
- O ligante está próximo de metais, caso existam?
- A pose parece compatível com uma interação biologicamente plausível?
- A molécula teve energia boa, mas ficou fora do sítio ativo?

Importante:
Não criar interface gráfica.
Não implementar sítios alostéricos.
Não implementar dinâmica molecular.
Não implementar ADMET.
Não implementar algoritmo genético.
Não refazer o docking.
Não alterar o funcionamento principal do pipeline antigo.
Não remover scripts existentes.

Estrutura desejada:

Criar:

src/
  betalactamase_engine/
    scoring/
      __init__.py
      interaction_scoring.py

Criar também:

scripts/
  score_interactions.py

Caso a pasta `src/betalactamase_engine/scoring/` ainda não exista, criá-la.

Entradas esperadas:

O script `scripts/score_interactions.py` deve receber:

--target-profile
Caminho para o arquivo target_profile.json.

--protein
Caminho para o arquivo PDB da proteína/receptor.

--ligand-pose
Caminho para a pose dockada do ligante.

--out
Caminho para salvar o resultado da análise.

Argumentos opcionais:

--compound-id
Identificador do composto. Se não for informado, usar o nome do arquivo do ligante.

--site
Nesta etapa, aceitar apenas:
active_site

Valor padrão:
active_site

--distance-threshold
Distância máxima, em Å, para considerar contato simples entre ligante e resíduo.
Padrão: 4.0

--center-threshold
Distância máxima, em Å, entre o centro do ligante e o centro do sítio ativo para considerar que o ligante está dentro do sítio.
Padrão: 8.0

--metal-threshold
Distância máxima, em Å, para considerar contato com metal.
Padrão: 3.0

--format
Formato de saída:
json ou csv.
Padrão: json.

--verbose
Exibir resumo da análise no terminal.

Formatos de ligante:
O módulo deve tentar suportar, nesta ordem:
- PDBQT;
- PDB;
- SDF;
- MOL2.

Se o suporte a algum formato for difícil, implementar pelo menos PDBQT e PDB nesta etapa, pois são comuns no pipeline de docking.

O parser de ligante deve:
- extrair coordenadas atômicas;
- ignorar linhas inválidas;
- lidar com arquivos vazios;
- retornar erro amigável se não encontrar átomos válidos.

Requisitos do módulo interaction_scoring.py:

1. Ler o `target_profile.json`.

2. Ler o PDB da proteína.

3. Ler a pose dockada do ligante.

4. Obter do `target_profile.json`:
   - centro do active_site;
   - resíduos catalíticos candidatos;
   - metais detectados;
   - cadeia selecionada, se houver;
   - raio do sítio ativo.

5. Calcular o centro geométrico do ligante dockado.

6. Calcular a distância entre:
   - centro do ligante;
   - centro do sítio ativo.

7. Identificar resíduos da proteína próximos ao ligante.

Um resíduo deve ser considerado próximo se qualquer átomo do resíduo estiver a uma distância menor ou igual ao valor de `--distance-threshold`, padrão 4.0 Å, de qualquer átomo do ligante.

8. Identificar contatos com resíduos catalíticos candidatos.

Um contato catalítico ocorre quando o resíduo próximo ao ligante também aparece em `active_site.candidate_catalytic_residues` do target_profile.json.

9. Identificar contatos com metais.

Um contato com metal ocorre quando qualquer átomo do ligante estiver a distância menor ou igual ao valor de `--metal-threshold`, padrão 3.0 Å, de um metal listado no target_profile.json ou detectado no PDB.

10. Calcular um `interaction_score` entre 0.0 e 1.0.

Sugestão de cálculo inicial:

- proximidade ao centro do sítio ativo: até 0.35
- contato com resíduos catalíticos: até 0.40
- contato com metais, se existirem: até 0.15
- número de contatos totais com o bolso: até 0.10

Regras sugeridas:

a) Proximidade ao centro:
Se distance_to_site_center <= center_threshold:
pontuação proporcional, maior quando mais próximo.
Se estiver acima do center_threshold:
pontuação baixa ou zero.

b) Contatos catalíticos:
Se houver pelo menos 1 contato com resíduo catalítico candidato:
dar bônus.
Se houver 2 ou mais:
dar bônus maior.
Se não houver:
pontuação zero nessa parte.

c) Metais:
Se existirem metais no target_profile e houver contato:
dar bônus.
Se existirem metais e não houver contato:
não dar bônus.
Se não existirem metais:
não penalizar fortemente; apenas não usar esse componente ou redistribuir peso.

d) Contatos totais:
Dar pequeno bônus se houver contatos com resíduos do bolso.

11. Classificar a pose.

Possíveis classificações:

- "mechanistically_plausible"
- "weak_active_site_pose"
- "outside_active_site"
- "unknown"

Regras sugeridas:

"mechanistically_plausible":
- ligante dentro do center_threshold;
- pelo menos um contato com resíduo catalítico candidato;
- interaction_score >= 0.65.

"weak_active_site_pose":
- ligante dentro do center_threshold;
- poucos ou nenhum contato catalítico;
- interaction_score entre 0.35 e 0.65.

"outside_active_site":
- distância do centro do ligante ao centro do sítio ativo maior que center_threshold.

"unknown":
- dados insuficientes;
- ausência de átomos;
- ausência de centro do sítio ativo;
- erro parcial de parsing.

12. A saída JSON deve seguir estrutura semelhante a:

{
  "compound_id": "ligand_001",
  "evaluated_site": "active_site",
  "site_id": "active_site",
  "ligand_pose_file": "results/poses/ligand_001.pdbqt",
  "protein_file": "data/structures/target.pdb",
  "target_profile_file": "results/experiments/teste/target_profile.json",
  "ligand_center": [0.0, 0.0, 0.0],
  "site_center": [0.0, 0.0, 0.0],
  "distance_to_site_center": 0.0,
  "distance_threshold": 4.0,
  "center_threshold": 8.0,
  "metal_threshold": 3.0,
  "nearby_residue_contacts": [
    {
      "chain": "A",
      "residue_name": "SER",
      "residue_id": 70,
      "min_distance": 3.1
    }
  ],
  "critical_residue_contacts": [
    {
      "chain": "A",
      "residue_name": "SER",
      "residue_id": 70,
      "min_distance": 3.1
    }
  ],
  "metal_contacts": [],
  "contact_summary": {
    "nearby_residue_contact_count": 0,
    "critical_residue_contact_count": 0,
    "metal_contact_count": 0
  },
  "interaction_score": 0.0,
  "pose_classification": "unknown",
  "warnings": []
}

13. Se o formato de saída for CSV, gerar pelo menos as colunas:

compound_id,
evaluated_site,
pose_classification,
interaction_score,
distance_to_site_center,
nearby_residue_contact_count,
critical_residue_contact_count,
metal_contact_count,
ligand_pose_file,
protein_file,
target_profile_file,
warnings

14. Criar script `scripts/score_interactions.py`.

Exemplo de uso:

python scripts/score_interactions.py \
  --target-profile results/experiments/teste/target_profile.json \
  --protein data/structures/betalac13.pdb \
  --ligand-pose results/poses/ligand_001.pdbqt \
  --compound-id ligand_001 \
  --out results/experiments/teste/interaction_score.json \
  --verbose

15. O resumo em terminal com `--verbose` deve exibir algo como:

Interaction scoring completed.
Compound: ligand_001
Evaluated site: active_site
Distance to active site center: 4.2 Å
Nearby residue contacts: 8
Critical residue contacts: 2
Metal contacts: 0
Interaction score: 0.78
Pose classification: mechanistically_plausible
Output: results/experiments/teste/interaction_score.json

16. Tratamento de erros obrigatório:

O código deve lidar de forma amigável com:
- target_profile.json inexistente;
- target_profile.json inválido;
- ausência de active_site.center;
- protein PDB inexistente;
- ligand pose inexistente;
- ligand pose sem átomos válidos;
- resíduos catalíticos ausentes no JSON;
- metais ausentes;
- formato de ligante não suportado.

17. Não afirmar que a molécula é inibidor real.

A documentação e nomes devem usar termos como:
- candidato;
- pose compatível;
- plausibilidade preliminar;
- análise geométrica;
- contato provável.

Evitar frases definitivas como:
- "a molécula inibe a enzima";
- "a molécula é ativa";
- "a molécula bloqueia a beta-lactamase".

18. Atualizar documentação.

Criar ou atualizar:

docs/interaction_scoring.md

A documentação deve explicar:

- objetivo do módulo;
- por que energia de docking sozinha não basta;
- como o score de interação complementa o AutoDock Vina;
- quais entradas são usadas;
- como interpretar `interaction_score`;
- como interpretar `pose_classification`;
- exemplo de comando;
- exemplo resumido de saída JSON;
- limitações da abordagem.

19. Atualizar README, se existir seção de pipeline.

Adicionar uma etapa opcional depois do docking:

Etapa 2 — Análise pós-docking de interação:

python scripts/score_interactions.py \
  --target-profile results/experiments/teste/target_profile.json \
  --protein data/structures/betalac13.pdb \
  --ligand-pose results/poses/ligand_001.pdbqt \
  --out results/experiments/teste/interaction_score.json \
  --verbose

20. Testes manuais mínimos:

Após implementar, rodar pelo menos um teste com um arquivo PDB real existente em `data/structures/`.

Se não houver ligand pose disponível, criar um pequeno arquivo de pose de teste sintético somente para validar parsing e execução, sem afirmar valor científico.

21. Critérios de aceite:

A tarefa estará concluída quando:

a) O comando abaixo funcionar sem erro:

python scripts/score_interactions.py \
  --target-profile results/experiments/teste/target_profile.json \
  --protein data/structures/betalac13.pdb \
  --ligand-pose CAMINHO_DE_UMA_POSE_EXISTENTE.pdbqt \
  --out results/experiments/teste/interaction_score.json \
  --verbose

b) O arquivo de saída for criado corretamente.

c) O JSON de saída contiver:
- compound_id;
- evaluated_site;
- ligand_center;
- site_center;
- distance_to_site_center;
- nearby_residue_contacts;
- critical_residue_contacts;
- metal_contacts;
- interaction_score;
- pose_classification;
- warnings.

d) O pipeline antigo continuar funcionando.

e) Nenhuma funcionalidade existente deve ser removida.

f) O projeto não deve usar git worktree.

22. Observação final:

Essa etapa é preparatória para o ranking multiobjetivo.

O próximo passo futuro será combinar:
- docking_score;
- interaction_score;
- ADMET;
- propriedades físico-químicas;
- penalizações por pose ruim;

em um futuro módulo chamado `consensus_score.py`.

Não implemente `consensus_score.py` agora.