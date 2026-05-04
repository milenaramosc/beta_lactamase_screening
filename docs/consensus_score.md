# Etapa 3 — Ranking multiobjetivo (consensus score)

Esta etapa combina o ranking de docking com informações pós-docking e, quando disponíveis, dados ADMET. O objetivo é priorizar candidatos de forma mais robusta do que o docking isolado.

## Objetivo

- Normalizar o score do AutoDock Vina.
- Incorporar interaction_score e pose_classification.
- Incluir ADMET, se houver arquivo disponível.
- Gerar um ranking final com classificação conservadora.

## Entradas

Obrigatórias:

- `--ranking`: CSV com ranking do docking (ex: `results/ranking.csv`).
- `--out`: caminho para o CSV final.

Opcionais:

- `--interactions`: CSV, JSON único ou diretório com JSONs do interaction scoring.
- `--admet`: CSV de ADMET existente.
- `--compound-column`: nome da coluna de identificador do composto.
- `--vina-column`: nome da coluna de score de docking.
- `--top-n`: limitar número de linhas no output.
- `--min-interaction-score`: marca compostos abaixo desse score como penalizados.
- `--verbose`: imprime resumo no terminal.

## Como o docking_norm é calculado

Como valores de docking são mais negativos quando melhores, o script converte o score para uma escala 0.0–1.0 onde:

- docking melhor -> próximo de 1.0
- docking pior -> próximo de 0.0

Se todos os valores forem iguais, o script usa 0.5 para todos e registra warning.

## Como interaction_score e pose_classification entram

- `interaction_score` entra diretamente na composição do `final_score`.
- `pose_classification` aplica um bônus/penalização conservador.

Classificações e bônus padrão:

- `mechanistically_plausible`: 1.0
- `weak_active_site_pose`: 0.7
- `outside_active_site`: 0.2
- `unknown`: 0.4

## Como o ADMET é usado

- Se existir `admet_score`, ele é normalizado para 0.0–1.0.
- Caso contrário, o script usa um componente simples com base em:
  - Lipinski violations
  - PAINS
  - QED (se disponível)

Se `--admet` não for fornecido, o peso é redistribuído e o ranking segue com docking + interação.

## Interpretação do final_score

`final_score` varia de 0.0 a 1.0 e é calculado com pesos conservadores. A classificação final é:

- `high_priority_candidate`: final_score >= 0.75 e pose plausível
- `medium_priority_candidate`: final_score >= 0.55
- `low_priority_candidate`: final_score >= 0.35
- `deprioritized_candidate`: final_score < 0.35

## Exemplo de comando

```bash
python scripts/consensus_score.py \
  --ranking results/ranking.csv \
  --interactions results/experiments/teste/interaction_scores.csv \
  --out results/experiments/teste/final_candidates.csv \
  --verbose
```

## Limitações

- A abordagem é preliminar e não confirma atividade biológica.
- ADMET é consumido apenas se fornecido pelo usuário.
- Não substitui validação experimental.
