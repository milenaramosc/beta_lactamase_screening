# Etapa 2 — Análise pós-docking de interação

Esta etapa avalia a plausibilidade geométrica da pose do ligante no sítio ativo, complementando o escore do AutoDock Vina. A energia de docking sozinha não indica se a pose ficou próxima ao centro do sítio, se encostou em resíduos catalíticos ou se interagiu com metais relevantes. O `interaction_score` fornece uma avaliação adicional baseada em distâncias e contatos estruturais.

## Objetivo

- Verificar se a pose do ligante está dentro do sítio ativo.
- Detectar contatos com resíduos catalíticos candidatos.
- Detectar contatos com metais quando presentes.
- Gerar uma classificação conservadora da pose.

## Entradas

- `--target-profile`: arquivo `target_profile.json` gerado pelo Target Profiler.
- `--protein`: arquivo PDB da proteína (receptor).
- `--ligand-pose`: pose dockada do ligante (PDBQT, PDB, SDF ou MOL2).
- `--out`: caminho do arquivo de saída (JSON ou CSV).

Argumentos opcionais:

- `--compound-id`: identificador do composto (padrão: nome do arquivo do ligante).
- `--site`: aceita apenas `active_site` nesta etapa.
- `--distance-threshold`: contato simples ligante-resíduo (padrão 4.0 Å).
- `--center-threshold`: distância máxima do centro do sítio (padrão 8.0 Å).
- `--metal-threshold`: distância máxima para contato com metal (padrão 3.0 Å).
- `--format`: `json` ou `csv` (padrão `json`).
- `--verbose`: mostra resumo no terminal.

## Interpretação do interaction_score

O escore vai de 0.0 a 1.0 e considera:

- Proximidade ao centro do sítio (até 0.35).
- Contatos com resíduos catalíticos (até 0.40).
- Contatos com metais (até 0.15), apenas se metais existirem.
- Número total de contatos com resíduos do bolso (até 0.10).

Se não houver metais no alvo, essa parte não penaliza a pose.

## Interpretação da pose_classification

- `mechanistically_plausible`: pose dentro do sítio e contato catalítico, com escore alto.
- `weak_active_site_pose`: pose dentro do sítio, mas com poucos contatos catalíticos.
- `outside_active_site`: pose fora do centro do sítio ativo.
- `unknown`: dados insuficientes ou erro parcial de parsing.

## Exemplo de comando

```bash
python scripts/score_interactions.py \
  --target-profile results/experiments/teste/target_profile.json \
  --protein data/structures/betalac13.pdb \
  --ligand-pose results/top50_poses/BETALAC13.3_CYS34-CHARGED__CHEMBL122450.pdbqt \
  --compound-id ligand_001 \
  --out results/experiments/teste/interaction_score.json \
  --verbose
```

## Exemplo resumido de saída JSON

```json
{
  "compound_id": "ligand_001",
  "evaluated_site": "active_site",
  "site_id": "active_site",
  "ligand_center": [0.0, 0.0, 0.0],
  "site_center": [0.0, 0.0, 0.0],
  "distance_to_site_center": 4.2,
  "nearby_residue_contacts": [],
  "critical_residue_contacts": [],
  "metal_contacts": [],
  "interaction_score": 0.42,
  "pose_classification": "weak_active_site_pose",
  "warnings": []
}
```

## Limitações

- A análise é geométrica e não confirma atividade biológica.
- Não substitui validação experimental.
- Não inclui dinâmica molecular, ADMET ou sítios alostéricos.

Esta etapa deve ser usada como triagem preliminar de plausibilidade estrutural.
