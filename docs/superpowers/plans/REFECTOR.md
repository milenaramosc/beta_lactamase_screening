IMPORTANTE:
Não implemente interface.
Não implemente docking novo.
Não implemente ADMET.
Não implemente algoritmo genético agora.
Esta tarefa é somente a criação do Target Profiler e do target_profile.json.
O objetivo é preparar a base científica para as próximas etapas.

Você é um agente de desenvolvimento especializado em Python, bioinformática, quimioinformática, docking molecular e organização de pipelines científicos reprodutíveis.

Contexto do projeto:
O projeto atual é um motor em terminal para triagem e geração de candidatos a inibidores de beta-lactamase. Ele trabalha com arquivos .PDB de beta-lactamases, identifica sítio ativo, coleta compostos, prepara ligantes, realiza docking molecular com AutoDock Vina, ranqueia moléculas candidatas e possui uma etapa inicial de geração de candidatos de novo com algoritmo genético.

Objetivo geral:
Evoluir a base científica do projeto antes de evoluir interface. Não criar frontend, dashboard ou interface gráfica. A prioridade é transformar o projeto em um motor científico mais robusto, reprodutível e defensável para dissertação.

Objetivo desta tarefa:
Implementar a Etapa 1 — Corrigir base científica — criando um módulo de caracterização do alvo chamado Target Profiler e preparando o projeto para que a geração de candidatos de novo possa ser guiada por características bioquímicas e mecanísticas da beta-lactamase.

Requisitos principais:
1. Criar um módulo responsável por analisar o arquivo .PDB da beta-lactamase.
2. Detectar automaticamente informações estruturais relevantes do alvo.
3. Gerar um arquivo `target_profile.json` que será usado pelas próximas etapas do pipeline.
4. Não alterar o funcionamento principal do pipeline atual de docking.
5. Manter tudo funcionando via terminal.
6. Preservar compatibilidade com a estrutura atual do projeto.
7. Evitar dependência obrigatória de serviços web externos.

Estrutura desejada:
Criar preferencialmente uma estrutura como:

src/
  betalactamase_engine/
    __init__.py
    target_profiler.py
    receptor_preparation.py
    utils/
      __init__.py
      pdb_utils.py

scripts/
  profile_target.py

Caso o projeto atual não tenha a pasta `src/`, criar a estrutura de forma limpa sem quebrar os scripts existentes.

Implementação obrigatória do Target Profiler:

O módulo `target_profiler.py` deve receber um arquivo `.pdb` e extrair:

1. Identificação básica:
   - caminho do arquivo PDB;
   - nome/base do arquivo;
   - número de cadeias;
   - número de resíduos;
   - número de átomos.

2. Detecção de resíduos catalíticos prováveis:
   - procurar resíduos típicos de serina-beta-lactamases, especialmente:
     - SER
     - LYS
     - GLU
     - ASN
     - ARG
   - quando possível, identificar resíduos compatíveis com motivos catalíticos de beta-lactamases classe A, C ou D.
   - não assumir rigidamente numeração exata, pois diferentes PDBs podem ter numerações diferentes.
   - detectar padrões aproximados por proximidade espacial e tipo de resíduo.

3. Detecção de metais:
   - identificar íons metálicos no PDB, especialmente:
     - ZN
     - MG
     - MN
     - FE
     - CA
   - se houver ZN próximo ao sítio ativo, marcar o alvo como possível metallo-beta-lactamase/classe B.
   - preservar essa informação para etapas futuras de docking e preparo de receptor.

4. Detecção de HETATM:
   - listar ligantes, cofatores, metais e águas encontrados no PDB.
   - separar:
     - água;
     - metais;
     - ligantes orgânicos;
     - outros heteroátomos.
   - não remover nada nesta etapa; apenas analisar e reportar.

5. Estimativa do centro do sítio ativo:
   - se o PDB possuir registros SITE, usar essa informação.
   - se houver ligante co-cristalizado orgânico, usar o centro geométrico do ligante como referência.
   - se houver metais catalíticos, usar o centro ao redor dos metais como referência.
   - se não houver nada disso, usar um fallback baseado em resíduos candidatos do sítio ativo.
   - documentar no JSON qual método foi usado:
     - "SITE_RECORD"
     - "COCRYSTAL_LIGAND"
     - "METAL_CENTER"
     - "CATALYTIC_RESIDUES"
     - "GEOMETRIC_FALLBACK"

6. Caracterização simples do bolso:
   - calcular resíduos próximos ao centro do sítio ativo dentro de um raio configurável, por padrão 8 Å.
   - contar resíduos:
     - polares;
     - apolares;
     - carregados positivamente;
     - carregados negativamente;
     - aromáticos.
   - gerar uma descrição simples do ambiente:
     - hydrophobic_ratio;
     - polar_ratio;
     - charged_ratio;
     - aromatic_ratio.

7. Inferência preliminar da classe da beta-lactamase:
   - se houver ZN ou outro metal catalítico relevante no sítio ativo, classificar como:
     - "possible_class_B_metallo_beta_lactamase"
   - se houver agrupamento de SER/LYS/GLU/ASN compatível com mecanismo de serina-beta-lactamase, classificar como:
     - "possible_serine_beta_lactamase"
   - caso não seja possível inferir:
     - "unknown_beta_lactamase_class"
   - essa classificação deve ser conservadora, ou seja, não fazer afirmações fortes sem evidência suficiente.

8. Saída obrigatória:
Gerar um arquivo JSON com estrutura semelhante a:

{
  "target": {
    "pdb_file": "data/structures/example.pdb",
    "pdb_name": "example",
    "chains": ["A"],
    "atom_count": 0,
    "residue_count": 0
  },
  "classification": {
    "predicted_class": "possible_serine_beta_lactamase",
    "confidence": "low|medium|high",
    "evidence": []
  },
  "active_site": {
    "center": [0.0, 0.0, 0.0],
    "detection_method": "COCRYSTAL_LIGAND",
    "radius_angstrom": 8.0,
    "candidate_catalytic_residues": []
  },
  "metals": [],
  "heteroatoms": {
    "waters": [],
    "metals": [],
    "organic_ligands": [],
    "others": []
  },
  "pocket_properties": {
    "nearby_residue_count": 0,
    "hydrophobic_ratio": 0.0,
    "polar_ratio": 0.0,
    "charged_ratio": 0.0,
    "aromatic_ratio": 0.0,
    "positive_residue_count": 0,
    "negative_residue_count": 0,
    "polar_residue_count": 0,
    "hydrophobic_residue_count": 0,
    "aromatic_residue_count": 0
  },
  "warnings": []
}

Criar script de terminal:

Criar `scripts/profile_target.py` para executar a análise:

python scripts/profile_target.py \
  --pdb data/structures/alvo.pdb \
  --out results/experiments/teste_target_profile/target_profile.json

Argumentos obrigatórios:
- `--pdb`: caminho do arquivo PDB.
- `--out`: caminho do JSON de saída.

Argumentos opcionais:
- `--radius`: raio em Å para análise do sítio ativo, padrão 8.0.
- `--verbose`: exibir resumo no terminal.

Ao final, o script deve imprimir um resumo como:

Target profile generated successfully.
PDB: alvo.pdb
Predicted class: possible_serine_beta_lactamase
Active site method: COCRYSTAL_LIGAND
Active site center: x, y, z
Metals detected: 0
Candidate catalytic residues: N
Output: results/experiments/.../target_profile.json

Regras importantes:
1. Não criar interface gráfica.
2. Não remover funcionalidades existentes.
3. Não quebrar os scripts atuais.
4. Não fazer docking nesta etapa.
5. Não fazer ADMET nesta etapa.
6. Não implementar algoritmo genético nesta etapa.
7. Esta etapa serve apenas para corrigir e fortalecer a base científica do alvo.
8. Usar código limpo, modular e testável.
9. Adicionar tratamento de erro amigável para:
   - arquivo PDB inexistente;
   - PDB vazio;
   - PDB sem átomos válidos;
   - ausência de ligantes;
   - ausência de SITE;
   - ausência de metais.
10. O código deve continuar funcionando mesmo quando algumas informações não forem encontradas.

Dependências:
Preferir usar Biopython/Bio.PDB se já existir no projeto. Caso não exista, verificar o `requirements.txt` e adicionar `biopython` se necessário. Usar somente dependências leves e adequadas para ambiente terminal.

Também pode usar:
- pathlib
- json
- argparse
- math
- statistics
- dataclasses
- typing

Evitar dependências pesadas nesta etapa.

Critérios de aceite:
A tarefa estará concluída quando:

1. O comando abaixo funcionar:

python scripts/profile_target.py --pdb data/structures/algum_alvo.pdb --out results/experiments/teste/target_profile.json --verbose

2. O arquivo `target_profile.json` for criado corretamente.

3. O JSON conter:
   - classificação preliminar;
   - centro do sítio ativo;
   - método usado para detectar o sítio ativo;
   - metais detectados;
   - HETATM organizados;
   - resíduos candidatos;
   - propriedades simples do bolso;
   - warnings quando algo não for encontrado.

4. O pipeline antigo continuar funcionando.

5. O código estiver preparado para que, em uma próxima etapa, o algoritmo genético use o `target_profile.json` como entrada.

Após implementar, gerar também um pequeno exemplo de uso no README ou em um arquivo `docs/target_profiler.md`, contendo:

- objetivo do Target Profiler;
- comando de execução;
- exemplo resumido de saída;
- explicação de como esse JSON será usado futuramente pela geração de candidatos de novo.

Não inventar resultados científicos. Qualquer inferência deve ser marcada como provável, possível ou desconhecida.
