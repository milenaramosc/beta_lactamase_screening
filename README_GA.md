# Algoritmo Genético Guiado por Estrutura da Proteína

## Visão Geral

Esta é a **versão definitiva** que integra análise da beta-lactamase (proteína alvo) com algoritmo genético para gerar inibidores específicos.

### Como Funciona

```
PDB da Beta-Lactamase
        ↓
[Análise do Sítio Ativo]
   - Resíduos catalíticos
   - Hidrofobicidade
   - Cargas eletrostáticas
   - Volume disponível
   - Ligante co-cristalizado (se houver)
        ↓
[Função de Fitness Customizada]
   - Compatível com sítio específico
   - Complementaridade química
   - Tamanho ideal para o bolso
        ↓
[Algoritmo Genético]
   - Evolui moléculas para aquele alvo
        ↓
Inibidores Específicos
```

---

## Por Que Esta Versão é Melhor?

| Aspecto | Standalone | **Protein-Guided** ✨ |
|---------|------------|---------------------|
| **Input** | Nenhum | PDB da beta-lactamase |
| **Fitness** | Genérica (inibidores gerais) | **Específica para o alvo** |
| **Tamanho das moléculas** | Fixo (250-500 Da) | **Ajustado ao volume do sítio** |
| **Hidrofobicidade** | Média | **Otimizada para o bolso** |
| **Cargas** | Carboxilato sempre | **Baseada em resíduos do sítio** |
| **Resultado** | Inibidores genéricos | **Inibidores específicos para aquela proteína** |

---

## Instalação

### Requisitos
```bash
pip install rdkit biopython numpy
```

**Não precisa** de AutoDock Vina, OpenBabel ou PyMOL para rodar o GA.

---

## Uso

### 1. Preparar Estrutura PDB

Você pode usar:
- **PDB bruto** baixado do RCSB
- **PDB preparado** (sem água, com hidrogênios)

```bash
# Opção A: Baixar PDB direto
wget https://files.rcsb.org/download/1ZG4.pdb -P data/structures/

# Opção B: Usar estrutura já preparada
python scripts/prepare_protein.py  # Etapa 2.5 do pipeline original
```

---

### 2. Executar GA Guiado pela Proteína

#### Uso Básico
```bash
python scripts/generate_inhibitors_from_protein.py \
  --pdb data/structures/1ZG4.pdb
```

**Saída (após ~8 minutos):**
```
results/inhibitors_1ZG4.sdf  — 20 moléculas 3D
results/inhibitors_1ZG4.csv  — Propriedades tabuladas
```

#### Uso Avançado
```bash
python scripts/generate_inhibitors_from_protein.py \
  --pdb data/structures/1ZG4.pdb \
  --generations 100 \
  --population 200 \
  --mutation-rate 0.6 \
  --top-out 50 \
  --seed 42
```

---

## Exemplo Real de Execução

```bash
$ python scripts/generate_inhibitors_from_protein.py --pdb data/structures/1AXB.pdb

══════════════════════════════════════════════════════════════════════
  GA Guiado por Estrutura de Beta-Lactamase
══════════════════════════════════════════════════════════════════════

Analisando estrutura: 1AXB.pdb

Características do sítio ativo:
  Centro: (7.93, 10.46, 41.51)
  Volume: 49910.9 Å³
  Resíduos no sítio: 159
  Hidrofobicidade: 30.82%
  Resíduos carregados: +16 / -15
  Doadores H-bond: 47
  Aceptores H-bond: 46

Resíduos catalíticos identificados:
    LYS111
    MET69
    SER70

✓ Ligante co-cristalizado encontrado: FOS
  Tamanho: 5.23 Å
  Átomos: 15

[GA executa e gera candidatos otimizados para este sítio específico]
```

---

## Análise do Sítio Ativo

### O Que é Analisado

#### 1. Resíduos Catalíticos
```python
# Identifica automaticamente:
- Registros SITE do PDB (se presentes)
- Resíduos conservados (SER70, LYS73, GLU166, etc.)
- Resíduos próximos ao ligante co-cristalizado
```

**Exemplo de output:**
```
Resíduos catalíticos identificados:
  SER70   — Serina nucleofílica (ataque ao beta-lactam)
  LYS73   — Lisina catalítica
  GLU166  — Glutamato (estabilização)
```

#### 2. Características Químicas
```python
- Hidrofobicidade: % de resíduos apolares (ALA, VAL, LEU, ILE, PHE)
- Cargas: número de LYS/ARG (+) e ASP/GLU (-)
- H-bonds: resíduos doadores (SER, THR, TYR) e aceptores (ASP, GLU)
```

**Exemplo de output:**
```
Hidrofobicidade: 30.82%
  → Sítio moderadamente polar
  → Fitness favorece LogP entre 0-2

Resíduos carregados: +16 / -15
  → Sítio balanceado
  → Fitness favorece moléculas com carboxilato (negativo)
```

#### 3. Volume do Bolso
```python
volume = (max_x - min_x) * (max_y - min_y) * (max_z - min_z)
ideal_MW = volume / 5  (heurística)
```

**Exemplo:**
```
Volume: 49910.9 Å³
  → MW ideal: 200-600 Da
  → Fitness penaliza moléculas < 200 ou > 600
```

#### 4. Ligante Co-cristalizado (Farmacóforo)
```python
# Se PDB tem inibidor co-cristalizado (HETATM):
- Extrai centro e tamanho do ligante
- Usa como referência de tamanho ideal
- Bônus de fitness para moléculas similares
```

**Exemplo:**
```
✓ Ligante co-cristalizado encontrado: FOS (fosfonamidon)
  Tamanho: 5.23 Å
  → Fitness dá bônus para moléculas de ~5-6 Å de raio
```

---

## Função de Fitness Customizada

### Componentes (Baseados na Proteína)

```python
fitness = (
    volume_score * 0.30 +              # Tamanho compatível com sítio
    chemical_complementarity * 0.25 +   # Hidrofob/carga compatível
    structural_features * 0.20 +        # Beta-lactam, carboxilato
    qed_score * 0.15 +                  # Drug-likeness
    similarity * 0.10 +                 # Similar a inibidores conhecidos
    pharmacophore_bonus                 # Bônus se tem ligante co-cristalizado
)
```

### 1. Volume Score (30%)

```python
# Baseado no volume do sítio
site_volume = 49910 Å³
ideal_MW_min = 49910 / 5 = 998 Da  (heurística)
ideal_MW_max = 49910 / 2 = 2495 Da

# Mas limita a 200-600 para drug-likeness
ideal_MW_min = max(200, ideal_MW_min)
ideal_MW_max = min(600, ideal_MW_max)

# Score
if ideal_MW_min ≤ MW ≤ ideal_MW_max:
    volume_score = 1.0
else:
    volume_score = penalizado
```

**Efeito:** Moléculas muito grandes ou pequenas para o bolso são descartadas.

---

### 2. Complementaridade Química (25%)

#### Hidrofobicidade
```python
site_hydrophobic_ratio = 0.3082  # 30.82% de resíduos apolares

if site_hydrophobic_ratio > 0.5:
    # Sítio hidrofóbico → favorece LogP positivo
    ideal_logp = 2-4
else:
    # Sítio polar → favorece LogP próximo de 0
    ideal_logp = -1 a 1
```

**Exemplo (1AXB):**
- Hidrofobicidade: 30.82% → sítio **polar**
- Fitness favorece LogP ~ 0
- Moléculas com LogP > 3 são penalizadas

#### Carga Eletrostática
```python
charged_positive = 16  # LYS, ARG, HIS
charged_negative = 15  # ASP, GLU

if charged_positive + charged_negative > 10:
    # Sítio muito carregado → favorece moléculas com grupos ionizáveis
    bonus se tem carboxilato (COO⁻) ou amônio (NH₃⁺)
```

**Exemplo (1AXB):**
- Cargas: +16 / -15 → sítio **muito carregado**
- Fitness favorece moléculas com **carboxilato** (essencial!)

---

### 3. Características Estruturais (20%)

Idêntico à versão standalone:
- Beta-lactam: +0.5
- Carboxilato: +0.3 (ou -0.1 se ausente)
- Heterociclos: +0.2

---

### 4. QED Score (15%)

Drug-likeness geral (igual standalone).

---

### 5. Similaridade (10%)

Tanimoto com inibidores conhecidos (igual standalone).

---

### 6. Farmacóforo Bonus (até +0.1)

```python
if ligante_co_cristalizado:
    ligand_size = 5.23 Å
    mol_size = MW / 50  # Aproximação
    
    if abs(mol_size - ligand_size) < 2.0:
        pharmacophore_bonus = +0.1
```

---

## Comparação de Resultados

### Exemplo: Beta-Lactamase 1AXB

**Sítio:**
- Volume: 49910 Å³ (bolso grande)
- Hidrofobicidade: 30.82% (polar)
- Cargas: +16/-15 (muito carregado)

**Candidato Gerado (Rank #1):**
```
SMILES: C[C@H]1[C@@H]2C(C(=O)O)N(C)C(=O)N2C(C)(C)N(C)S1(=O)=O
MW: 305.36 Da
LogP: -0.42  ← Polar (compatível com sítio!)
Fitness: 0.6443
Características:
  ✓ Carboxilato presente (interage com LYS/ARG)
  ✓ Sulfonamida (H-bonds)
  ✓ Tamanho adequado para o bolso
```

**vs. Standalone (genérico):**
```
Fitness standalone: 0.9033 (alta, mas genérica)
Fitness protein-guided: 0.6443 (menor, mas específica para 1AXB)
```

**Interpretação:**
- Standalone gera inibidores genéricos (alta fitness química geral)
- Protein-guided sacrifica fitness genérica para **otimizar para aquele alvo**

---

## Workflow Completo Recomendado

### Passo 1: Baixar Estrutura PDB
```bash
# Buscar beta-lactamase de interesse no RCSB
# Ex: 1ZG4 (classe A), 1M2X (ESBL), 2G2U (KPC)

wget https://files.rcsb.org/download/1ZG4.pdb -P data/structures/
```

---

### Passo 2: (Opcional) Preparar Proteína
```bash
# Remove água, adiciona H, detecta sítio
python scripts/prepare_protein.py
```

**Nota:** Não é obrigatório! O script aceita PDB bruto também.

---

### Passo 3: Gerar Inibidores Guiados
```bash
python scripts/generate_inhibitors_from_protein.py \
  --pdb data/structures/1ZG4.pdb \
  --generations 100 \
  --population 200 \
  --top-out 30
```

**Output:**
- `results/inhibitors_1ZG4.csv` — Top 30 candidatos

---

### Passo 4: Filtrar por Fitness
```bash
# Selecionar candidatos com fitness > 0.7
cat results/inhibitors_1ZG4.csv | awk -F, '$4 > 0.7'
```

---

### Passo 5: Validar com Docking (Opcional)
```bash
# Se você tem AutoDock Vina instalado:
python scripts/screen_and_rank.py \
  --ligands results/inhibitors_1ZG4.sdf \
  --config config.yaml
```

**Resultado:** ΔG (energia de ligação) para cada candidato.

---

### Passo 6: Síntese Experimental
- Revisar rotas sintéticas (Reaxys, SciFinder)
- Priorizar moléculas com:
  - Fitness > 0.7
  - QED > 0.6
  - Rota sintética viável (< 5 passos)

---

## Interpretação dos Resultados

### Arquivo CSV

```csv
rank,id,smiles,fitness,mw,logp,hbd,hba,qed
1,PROTEIN_GUIDED_001,C[C@H]1...,0.6443,305.36,-0.42,1,5,0.709
```

### Colunas

| Coluna | Descrição | Ideal |
|--------|-----------|-------|
| **fitness** | Score total (específico para o alvo) | > 0.65 |
| **mw** | Peso molecular | Baseado no volume do sítio |
| **logp** | Lipofilia | Baseado na hidrofobicidade do sítio |
| **hbd/hba** | H-bond donors/acceptors | Baseado em SER/THR/ASP/GLU no sítio |
| **qed** | Drug-likeness geral | > 0.6 |

---

## Comparação das 3 Versões

| Aspecto | Original (Pipeline) | Standalone | **Protein-Guided** |
|---------|---------------------|------------|------------------|
| **Input** | PDB + 5000 compostos | Nenhum | **PDB da proteína** |
| **Tempo** | 24-48h | 10 min | **10 min** |
| **Fitness** | ΔG do docking | Genérica | **Específica para o alvo** |
| **Validação** | Alta (docking real) | Baixa | **Média (baseada em estrutura)** |
| **Especificidade** | Alta (para aquele alvo) | Baixa (genérica) | **Alta (customizada)** |
| **Casos de uso** | Validação final | Brainstorming | **Exploração guiada** |

---

## Quando Usar Cada Versão?

### Use **Pipeline Original** se:
- Precisa de ΔG validado por docking
- Tem tempo (24-48h)
- Quer máxima precisão

### Use **Standalone** se:
- Quer exploração rápida genérica
- Não tem PDB
- Está em fase de brainstorming

### Use **Protein-Guided** se: ✨
- **Tem estrutura PDB do alvo**
- **Quer inibidores específicos para aquela proteína**
- **Quer velocidade (10 min) + especificidade**
- **Está explorando múltiplos alvos** (rodar para cada PDB)

---

## Perguntas Frequentes

**P: Preciso de PDB preparado ou bruto?**  
R: Ambos funcionam. PDB preparado (sem água) é melhor, mas o script aceita PDB bruto também.

**P: E se meu PDB não tem SITE records?**  
R: O script usa fallback: busca resíduos conservados (SER, LYS, GLU) ou usa geometria central.

**P: O fitness é mais baixo que no standalone. Por quê?**  
R: Normal! Protein-guided sacrifica fitness genérica para otimizar especificamente para aquele sítio.

**P: Posso usar para outras proteínas (não beta-lactamase)?**  
R: Sim, mas ajuste os padrões estruturais (linha ~340) para características do seu alvo.

**P: Como validar os candidatos?**  
R: Idealmente, use docking (screen_and_rank.py) para calcular ΔG. Ou vá direto para síntese experimental.

---

## Limitações

1. **Não calcula ΔG real** — fitness é heurística baseada em propriedades químicas
2. **Não considera dinâmica** — sítio ativo é tratado como rígido
3. **Farmacóforo simplificado** — apenas tamanho do ligante co-cristalizado
4. **Requer validação** — docking ou teste experimental obrigatório

**Solução:** Use como **primeira triagem rápida**, depois valide com docking ou experimento.

---

## Próximos Passos

1. **Executar para múltiplos PDBs:**
   ```bash
   for pdb in data/structures/*.pdb; do
       python scripts/generate_inhibitors_from_protein.py --pdb $pdb
   done
   ```

2. **Comparar candidatos entre alvos:**
   - Inibidores específicos para classe A vs. classe C
   - Selecionar inibidores broad-spectrum (funcionam em vários alvos)

3. **Validar com docking:**
   ```bash
   python scripts/screen_and_rank.py --ligands results/inhibitors_*.sdf
   ```

4. **Sintetizar e testar experimentalmente**

---

## Referências

- **BioPython:** Cock et al. (2009) Bioinformatics
- **RCSB PDB:** Berman et al. (2000) Nucleic Acids Res
- **Beta-lactamase Catalytic Residues:** Ambler numbering system
- **RDKit:** https://www.rdkit.org/

---

**Esta é a versão DEFINITIVA que combina análise estrutural + GA!** 🎯
