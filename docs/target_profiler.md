# Target Profiler

## Objetivo
O **Target Profiler** é um módulo de caracterização bioquímica e estrutural do alvo (beta-lactamase). Ele analisa arquivos PDB para extrair informações críticas que guiam o restante do pipeline, especialmente a geração de candidatos *de novo*.

Sua principal função é transformar um arquivo de coordenadas brutas (.pdb) em um perfil estruturado (.json) contendo:
- Classificação provável da beta-lactamase (Serina vs Metalo).
- Localização e método de detecção do sítio ativo.
- Identificação de resíduos catalíticos candidatos.
- Presença de metais e heteroátomos (ligantes, águas).
- Propriedades físico-químicas do bolso (pocket characterization).

## Como o profiler identifica o sitio ativo
O profiler estima o sitio ativo de forma conservadora e reprodutivel, seguindo uma hierarquia:
1. Registros `SITE` do PDB (quando presentes).
2. Ligante co-cristalizado (maior ligante organico detectado).
3. Centro de metais, caso haja ions compativeis com metalo-beta-lactamases.
4. Residuos cataliticos candidatos (ex.: Ser, Lys, Glu) em proximidade.
5. Fallback geometric (centroide de todos os atomos da cadeia selecionada).

## Como o profiler detecta metais
Metais sao detectados usando uma lista explicita de ions conhecidos (ex.: ZN, MG, MN, FE, FE2, FE3, CA, CO, NI, CU, NA, K). Isso evita classificar qualquer HETATM curto como metal.

## Como o profiler representa aguas
As aguas sao resumidas no JSON com contagem total e quantas estao proximas ao sitio ativo dentro do raio configurado. Apenas essas aguas proximas sao listadas explicitamente.

## Execucao
Para gerar o perfil de um alvo, utilize o script `scripts/profile_target.py`:

```bash
python scripts/profile_target.py \
  --pdb data/structures/alvo.pdb \
  --chain A \
  --out results/experiments/nome_experimento/target_profile.json \
  --verbose
```

### Argumentos:
- `--pdb`: Caminho para o arquivo PDB de entrada.
- `--out`: Caminho para o arquivo JSON de saida.
- `--radius`: (Opcional) Raio em A para analise do bolso (padrao: 8.0).
- `--chain`: (Opcional) Cadeia especifica a ser analisada. Se omitido e houver varias cadeias, uma cadeia e selecionada automaticamente.
- `--verbose`: Exibe um resumo da analise no terminal.

## Interpretando quality_control
O bloco `quality_control` informa se o profiler encontrou registros `SITE`, ligantes organicos e metais, alem de indicar quando foi necessario usar fallback. Notas explicam a selecao automatica de cadeia e outras decisoes conservadoras.

## Exemplo de Saída
O arquivo `target_profile.json` gerado possui a seguinte estrutura (resumida):

```json
{
  "target": {
    "pdb_file": "data/structures/1BTL.pdb",
    "pdb_name": "1BTL",
    "atom_count": 2234,
    "residue_count": 463,
    "selected_chain": "A"
  },
  "classification": {
    "predicted_class": "possible_serine_beta_lactamase",
    "confidence": "medium",
    "evidence": ["Active site shows candidate Serine and Lysine residues"]
  },
  "active_site": {
    "center": [8.958, 10.722, 35.256],
    "detection_method": "COCRYSTAL_LIGAND",
    "candidate_catalytic_residues": [...]
  },
  "metals": [],
  "heteroatoms": {
    "waters": {
      "count": 199,
      "near_active_site_count": 5,
      "near_active_site": []
    },
    "metals": [],
    "organic_ligands": [],
    "others": []
  },
  "quality_control": {
    "has_site_records": true,
    "has_organic_ligand": true,
    "has_metals": false,
    "used_fallback": false,
    "notes": []
  },
  "pocket_properties": {
    "hydrophobic_ratio": 0.27,
    "polar_ratio": 0.36,
    "charged_ratio": 0.18
  }
}
```

## Uso Futuro
Este perfil será utilizado pela Etapa de Geração de Candidatos (Algoritmo Genético) para:
1. **Definir a Box de Docking:** O centro e o raio detectados eliminam a necessidade de definição manual de coordenadas.
2. **Guiar a Pontuação:** O Algoritmo Genético poderá usar a `predicted_class` para priorizar fragmentos que interagem melhor com Serinas ou Metais (como Zinco).
3. **Filtros de Bioafinidade:** Candidatos que não possuam complementaridade com as `pocket_properties` (ex: alta aromaticidade em bolsos hidrofóbicos) podem ser penalizados precocemente.
