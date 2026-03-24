# Design: Conclusão do Pipeline + Geração de Molécula Inibidora

**Data:** 2026-03-15
**Status:** Aprovado

---

## Contexto

Pipeline computacional em Python para triagem virtual de inibidores de beta-lactamase. O projeto já tem estrutura completa (etapas 1-6), mas faltam dois componentes críticos:

1. Etapa de **preparação da proteína** (entre etapas 2 e 4)
2. Etapa de **geração de molécula de novo** (após etapa 5 ADMET)
3. **Correção da fonte ZINC20** (inacessível por SSL neste ambiente)

---

## Componente 1: `scripts/prepare_protein.py` (Etapa 2.5)

### Posição no pipeline

```
Etapa 1: fetch_structures.py    → data/structures/*.pdb
Etapa 2: select_targets.py      → config.yaml (selected_targets)
Etapa 2.5: prepare_protein.py   → data/prepared/*.pdb + config.yaml (binding_site)
Etapa 4: screen_and_rank.py     → lê data/prepared/ com fallback para data/structures/
```

### Etapas internas (por alvo selecionado)

1. **Remover água cristalográfica** — elimina linhas `HETATM` com resíduo `HOH`/`WAT`
2. **Remover ligantes co-cristalizados** — remove todos os `HETATM` restantes
3. **Adicionar hidrogênios** — via OpenBabel (`-h`)
4. **Minimização de energia** — via OpenBabel com campo de força MMFF94, 500 passos
5. **Detecção do sítio ativo** — registros SITE do PDB + cálculo de centro/box; fallback ao centroide de CA atoms
6. **Salvar resultado** — `data/prepared/<ID>_prepared.pdb` + atualiza `config.yaml` com coordenadas do sítio

### Tratamento de erros

- Minimização falha → salva PDB sem minimização + aviso amarelo (não aborta)
- SITE records ausentes → usa centroide de CA + caixa 30 Å (comportamento atual)
- Arquivo `.pdb` não encontrado → erro + pula alvo

### Interface CLI

```bash
python scripts/prepare_protein.py
python scripts/prepare_protein.py --force   # reprocessa mesmo se já preparado
```

### Integração com `screen_and_rank.py`

- `prepare_receptor_pdbqt()` busca primeiro `data/prepared/<ID>_prepared.pdb`
- Fallback silencioso para `data/structures/<ID>.pdb` se não preparado
- Binding site usa `config.yaml binding_site.targets.<ID>` se disponível

### Dependências

- OpenBabel (já em `requirements.txt` via `openbabel-wheel`)
- BioPython (já no projeto)
- PyYAML (já no projeto)

---

## Componente 2: Correção de `collect_zinc()` em `fetch_compounds.py`

### Problema

`zinc20.docking.org` inacessível por erro de TLS neste ambiente. A função atual não tem paginação real (parâmetro `count` estava comentado), retornando sempre o mesmo lote.

### Solução

Substituir a fonte ZINC20 pela API do ChEMBL com filtros de propriedades físico-químicas equivalentes ao filtro drug-like do ZINC:

| Filtro | Valor |
|--------|-------|
| `mw_freebase` | ≤ 500 Da |
| `alogp` | ≤ 5 |
| `hbd` | ≤ 5 |
| `hba` | ≤ 10 |
| `molecule_type` | Small molecule |

Paginação por `offset` → diversidade real (1.3M+ moléculas disponíveis).

### O que não muda

- Interface: `--source zinc` continua funcionando
- Deduplicação por SMILES canônico
- Filtro de Lipinski com RDKit
- Formato de saída

---

## Componente 3: `scripts/generate_candidate.py` (Etapa 7)

### Objetivo

Após executar todo o pipeline (etapas 1-5), usar os resultados de docking como referência para **gerar uma nova molécula candidata** via abordagem genética.

### Abordagem: Algoritmo Genético com RDKit

1. **Semente inicial** — os top-N compostos do `results/ranking.csv` (melhor ΔG) são usados como população inicial
2. **Representação** — SMILES canônico + fragmentação por RDKit (BRICS/RECAP)
3. **Operadores genéticos:**
   - *Crossover*: recombinação de fragmentos entre dois compostos pai
   - *Mutação*: substituição de grupo funcional via transformações RDKit
4. **Função fitness** — composta por:
   - Score de similaridade com o sítio ativo (fingerprint Tanimoto)
   - Filtro de Lipinski (penalidade se falhar)
   - Estimativa de docking rápida (opcional, se Vina disponível)
5. **Saída** — top-5 candidatos gerados em `results/generated_candidates.sdf` + relatório

### Dependências adicionais

Nenhuma além das já listadas no projeto (RDKit já está em `requirements.txt`).

### Interface CLI

```bash
python scripts/generate_candidate.py
python scripts/generate_candidate.py --generations 50 --population 100
```

---

## Ordem de implementação

1. Corrigir `collect_zinc()` em `fetch_compounds.py`
2. Criar `scripts/prepare_protein.py`
3. Ajustar `screen_and_rank.py` para usar `data/prepared/`
4. Criar `scripts/generate_candidate.py`
5. Atualizar `README.md` com as novas etapas

---

## Critérios de sucesso

- `fetch_compounds.py --source zinc` coleta moléculas diversas com sucesso
- `prepare_protein.py` processa os 4 alvos selecionados em `config.yaml`
- `screen_and_rank.py` usa automaticamente os arquivos preparados
- `generate_candidate.py` produz pelo menos 5 candidatos novos com SMILES válido
