Você é um agente de desenvolvimento especializado em Python, bioinformática estrutural, docking molecular e pipelines científicos reprodutíveis.

Contexto:
O projeto já possui uma primeira versão do Target Profiler em:

src/betalactamase_engine/target_profiler.py
src/betalactamase_engine/utils/pdb_utils.py
scripts/profile_target.py
docs/target_profiler.md

O profiler já lê arquivos PDB, identifica informações básicas do alvo, detecta heteroátomos, metais, estima o sítio ativo e gera target_profile.json.

Objetivo desta tarefa:
Refinar cientificamente o Target Profiler para torná-lo mais robusto, reprodutível e adequado para as próximas etapas do pipeline: análise pós-docking, score de interação, ranking multiobjetivo e algoritmo genético guiado pelo alvo.

Importante:
Não implementar interface gráfica.
Não implementar sítios alostéricos agora.
Não implementar bolsões crípticos agora.
Não implementar dinâmica molecular agora.
Não implementar docking novo agora.
Não implementar ADMET agora.
Não implementar algoritmo genético agora.

Tarefas obrigatórias:

1. Manter o JSON com as seções principais já existentes:

{
  "target": {},
  "classification": {},
  "active_site": {},
  "metals": [],
  "heteroatoms": {},
  "pocket_properties": {},
  "warnings": []
}

2. Adicionar uma seção de controle de qualidade:

"quality_control": {
  "has_site_records": true,
  "has_organic_ligand": true,
  "has_metals": false,
  "used_fallback": false,
  "notes": []
}

3. Adicionar suporte ao argumento opcional:

--chain

Esse argumento deve permitir escolher uma cadeia específica da proteína.

Exemplo:

python scripts/profile_target.py \
  --pdb data/structures/1BTL.pdb \
  --chain A \
  --out results/experiments/teste/target_profile.json \
  --verbose

4. Se --chain for informado, analisar preferencialmente essa cadeia.

5. Se --chain não for informado e o PDB tiver múltiplas cadeias, escolher uma cadeia de forma consistente e registrar no JSON:

"selected_chain": "A"

Também adicionar uma nota em quality_control explicando que a cadeia foi selecionada automaticamente.

6. Corrigir o parser de registros SITE.

O método que lê registros SITE deve extrair corretamente todos os resíduos presentes em linhas SITE do PDB, não apenas parte deles.

7. Melhorar a representação de águas no JSON.

Em vez de listar todas as águas individualmente, usar uma estrutura resumida:

"waters": {
  "count": 199,
  "near_active_site_count": 5,
  "near_active_site": []
}

A lista near_active_site deve conter apenas águas próximas ao sítio ativo, dentro do raio configurado.

8. Melhorar a detecção de metais.

Usar uma lista explícita de íons/metais conhecidos, por exemplo:

ZN, MG, MN, FE, FE2, FE3, CA, CO, NI, CU, NA, K

Não classificar qualquer HETATM curto como metal.

9. Melhorar classificação preliminar da beta-lactamase.

A classificação deve continuar conservadora:

- "possible_serine_beta_lactamase"
- "possible_class_B_metallo_beta_lactamase"
- "unknown_beta_lactamase_class"

Também incluir:

"confidence": "low|medium|high"
"evidence": []

10. Não afirmar resultados científicos definitivos.

Usar sempre termos como:
- provável;
- possível;
- candidato;
- preliminar;
- inferido.

11. Remover arquivos __pycache__ e arquivos .pyc do projeto.

12. Atualizar .gitignore para ignorar:

__pycache__/
*.pyc
results/experiments/
*.log

13. Não remover arquivos de dados existentes sem autorização.

14. Atualizar docs/target_profiler.md explicando:

- objetivo do Target Profiler;
- como ele identifica o sítio ativo;
- como ele detecta metais;
- como ele representa águas;
- como usar --chain;
- como interpretar quality_control;
- exemplo de comando;
- exemplo resumido de saída JSON.

Critérios de aceite:

1. O comando abaixo deve funcionar:

python scripts/profile_target.py \
  --pdb data/structures/1BTL.pdb \
  --chain A \
  --out results/experiments/teste/target_profile.json \
  --verbose

2. O JSON gerado deve conter:
- target
- classification
- active_site
- metals
- heteroatoms
- pocket_properties
- quality_control
- warnings
- selected_chain ou informação equivalente

3. O pipeline antigo não pode quebrar.

4. Caso o PDB não tenha SITE, ligante ou metal, o profiler deve continuar funcionando com fallback e registrar isso em warnings ou quality_control.

5. Não implementar sítios alostéricos nesta etapa.