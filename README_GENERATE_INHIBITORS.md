# Geração de Inibidores de Beta-Lactamase com Algoritmo Genético

Algoritmo genético que analisa a estrutura 3D da beta-lactamase e gera inibidores específicos para aquele alvo.

---

## Instalação

```bash
pip install rdkit biopython numpy
```

---

## Uso Rápido

### 1. Obter estrutura PDB da beta-lactamase
```bash
# Exemplo: KPC-2
wget https://files.rcsb.org/download/1ZG4.pdb -P data/structures/
```

### 2. Gerar inibidores específicos
```bash
python scripts/generate_inhibitors_from_protein.py --pdb data/structures/1ZG4.pdb
```

**Saída (após ~10 minutos):**
- `results/inhibitors_1ZG4.sdf` — Estruturas 3D dos candidatos
- `results/inhibitors_1ZG4.csv` — Propriedades e fitness

---

## Como Funciona

### 1. Análise Automática do Sítio Ativo
```
PDB da Beta-Lactamase
        ↓
[Identifica automaticamente]
   - Resíduos catalíticos (SER70, LYS73, GLU166)
   - Hidrofobicidade do bolso
   - Cargas eletrostáticas
   - Volume do sítio
   - Ligante co-cristalizado (se houver)
        ↓
[Fitness Customizada]
   - Tamanho ideal para este bolso
   - LogP baseado na hidrofobicidade
   - Cargas complementares
        ↓
[Algoritmo Genético]
   - População inicial: inibidores conhecidos
   - Operadores: mutação (SMARTS) + crossover (BRICS)
   - Evolui por 50-100 gerações
        ↓
Inibidores Específicos para Aquele Alvo
```

### 2. Função de Fitness

```python
fitness = (
    volume_score * 0.30 +              # MW ideal para o bolso
    chemical_complementarity * 0.25 +   # Hidrofob/carga compatível
    structural_features * 0.20 +        # Beta-lactam, carboxilato
    qed_score * 0.15 +                  # Drug-likeness
    similarity * 0.10 +                 # Similar a inibidores conhecidos
    pharmacophore_bonus                 # Bônus se tem ligante
)
```

---

## Parâmetros

```bash
python scripts/generate_inhibitors_from_protein.py \
  --pdb data/structures/1ZG4.pdb \    # Arquivo PDB (obrigatório)
  --generations 100 \                 # Número de gerações (padrão: 50)
  --population 200 \                  # Tamanho da população (padrão: 100)
  --mutation-rate 0.6 \               # Taxa de mutação (padrão: 0.5)
  --top-out 30 \                      # Quantos salvar (padrão: 20)
  --seed 42                           # Semente aleatória (padrão: 42)
```

---

## Exemplo de Saída

```
══════════════════════════════════════════════════════════════════════
  GA Guiado por Estrutura de Beta-Lactamase
══════════════════════════════════════════════════════════════════════

Analisando estrutura: 1ZG4.pdb

Características do sítio ativo:
  Centro: (7.93, 10.46, 41.51)
  Volume: 49910.9 Å³
  Resíduos no sítio: 159
  Hidrofobicidade: 30.82%  ← Sítio polar
  Resíduos carregados: +16 / -15  ← Muito carregado
  Doadores H-bond: 47
  Aceptores H-bond: 46

Resíduos catalíticos identificados:
    SER70  — Serina nucleofílica
    LYS73  — Lisina catalítica
    GLU166 — Glutamato

✓ Ligante co-cristalizado encontrado: FOS
  Tamanho: 5.23 Å

[GA executa por 50 gerações]

Concluído!
  20 candidatos salvos em results/inhibitors_1ZG4.csv

Top 10 candidatos:
  #    Fitness    MW       LogP    QED     
  1    0.6443     305.36   -0.42   0.709   ← Polar (compatível!)
  2    0.6097     208.18   -0.72   0.595
  ...
```

---

## Interpretação dos Resultados

### CSV Gerado
```csv
rank,id,smiles,fitness,mw,logp,hbd,hba,qed,origin
1,PROTEIN_GUIDED_001,C[C@H]1...,0.6443,305.36,-0.42,1,5,0.709,seed:Tazobactam
```

### Critérios de Qualidade

| Propriedade | Ideal | Descrição |
|-------------|-------|-----------|
| **fitness** | > 0.65 | Score total (específico para o alvo) |
| **qed** | > 0.6 | Drug-likeness geral |
| **mw** | 250-500 | Peso molecular adequado |
| **logp** | Baseado no sítio | Hidrofilia/hidrofobicidade |
| **hbd/hba** | ≤5 / ≤10 | Lipinski compliance |

**Priorize candidatos com:**
- ✅ Fitness > 0.65
- ✅ QED > 0.6
- ✅ Presença de carboxilato (essencial para beta-lactamases)

---

## Estruturas PDB Sugeridas

### Beta-lactamases Classe A (Serina)
- **1ZG4** — KPC-2 (resistência a carbapenens)
- **1M2X** — ESBL CTX-M-9 (resistência a cefalosporinas)
- **1AXB** — TEM-1 (beta-lactamase clássica)
- **2G2U** — KPC-2 (alta resolução)

### Beta-lactamases Classe B (Metalo)
- **4S2I** — NDM-1 (requer inibidores não-beta-lactâmicos)
- **5N5I** — VIM-2

### Beta-lactamases Classe C (AmpC)
- **1KVL** — AmpC de E. coli

Busque mais em: https://www.rcsb.org/

---

## Próximos Passos

### 1. Validação Computacional (Opcional)
```bash
# Se você tem AutoDock Vina instalado:
python scripts/screen_and_rank.py --ligands results/inhibitors_1ZG4.sdf
```

### 2. Análise Química
- Verificar viabilidade sintética (Reaxys, SciFinder)
- Validar propriedades ADMET no SwissADME: http://www.swissadme.ch/
- Buscar substructures similares em bancos de dados

### 3. Síntese e Teste Experimental
**Priorize candidatos com:**
- Fitness > 0.70
- Estruturas sinteticamente acessíveis
- Scaffolds novos (potencial para patente)

**Testes in vitro:**
- Ensaios de inibição enzimática (IC₅₀)
- Testes de sinergia com antibióticos beta-lactâmicos
- Ensaios de citotoxicidade

---

## Comparação com Pipeline Original

| Aspecto | generate_inhibitors_from_protein.py | Pipeline Original (7 etapas) |
|---------|-------------------------------------|------------------------------|
| **Tempo** | 10 minutos | 24-48 horas |
| **Input** | PDB da proteína | PDB + 5000 compostos + docking |
| **Dependências** | RDKit + BioPython | Vina + OpenBabel + PyMOL |
| **Fitness** | Específica para o alvo | Baseada em ΔG do docking |
| **Validação** | Estrutural (heurística) | Alta (ΔG calculado) |
| **Uso** | Exploração guiada | Refinamento final |

**Recomendação:** Use esta versão para gerar candidatos rapidamente, depois valide os melhores com docking (Etapa 4 do pipeline original).

---

## Workflow Híbrido (Recomendado)

```bash
# 1. Gerar candidatos para vários alvos (1h)
python scripts/generate_inhibitors_from_protein.py --pdb data/structures/1ZG4.pdb
python scripts/generate_inhibitors_from_protein.py --pdb data/structures/4S2I.pdb

# 2. Filtrar manualmente (30 min)
# Selecionar top 20 por fitness + inspeção visual

# 3. Validar com docking (4h)
python scripts/screen_and_rank.py --ligands results/top20_selected.sdf

# 4. Propor para síntese
# Top 5 candidatos com melhor ΔG + fitness
```

**Tempo total:** ~6h (vs. 48h do pipeline completo)

---

## Troubleshooting

### Erro: ModuleNotFoundError
```bash
pip install rdkit biopython numpy
```

### Erro: PDB não encontrado
```bash
# Verificar caminho
ls -lh data/structures/

# Re-baixar
wget https://files.rcsb.org/download/1ZG4.pdb -P data/structures/
```

### Fitness muito baixo (< 0.5)
Normal! A fitness é específica para aquele alvo. Valores > 0.55 já são bons.

### Poucos candidatos gerados
```bash
# Aumentar população e gerações
python scripts/generate_inhibitors_from_protein.py \
  --pdb data/structures/1ZG4.pdb \
  --generations 100 \
  --population 200
```

---

## Detalhes Técnicos

### Identificação do Sítio Ativo

1. **SITE records do PDB** (prioridade)
2. **Resíduos catalíticos conservados** (SER70, LYS73, GLU166)
3. **Proximidade ao ligante co-cristalizado** (se houver)
4. **Fallback:** Centro geométrico da proteína

### Operadores Genéticos

**Mutação (SMARTS):**
- Halogênios: F ↔ Cl ↔ Br
- Grupos funcionais: CH₃ → OH → NH₂
- Ácidos/Amidas: COOH ↔ CONH₂
- Sulfonamidas: SO₂OH → SO₂NH₂

**Crossover (BRICS):**
- Fragmentação em ligações sinteticamente viáveis
- Recombinação de fragmentos de dois pais
- Reconstrói moléculas válidas

### População Inicial

6 inibidores clínicos conhecidos:
- Clavulanato (FDA approved)
- Sulbactam (FDA approved)
- Tazobactam (FDA approved)
- Avibactam (FDA approved)
- Relebactam (FDA approved)
- Vaborbactam (FDA approved)

---

## Limitações

1. **Não calcula ΔG real** — fitness é baseada em heurísticas químicas
2. **Sítio ativo estático** — não considera flexibilidade da proteína
3. **Requer validação** — docking ou teste experimental obrigatório
4. **Específico para beta-lactamases** — características estruturais hardcoded

---

## Referências

- **BioPython:** Cock et al. (2009) Bioinformatics
- **RDKit:** https://www.rdkit.org/
- **RCSB PDB:** Berman et al. (2000) Nucleic Acids Res
- **Beta-lactamase Catalytic Residues:** Ambler numbering
- **QED Score:** Bickerton et al. (2012) Nature Chemistry

---

## Documentação Completa

Para detalhes avançados, consulte `README_GA.md`

---

## Licença

Este código é fornecido para fins educacionais e de pesquisa. Candidatos gerados devem ser validados experimentalmente antes de qualquer aplicação clínica.
