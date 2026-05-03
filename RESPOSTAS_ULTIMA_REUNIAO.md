# Resposta às Perguntas da Ultima Reunião

**Documento técnico detalhando a implementação do Algoritmo Genético para geração de inibidores de beta-lactamase**

---

## 1. O que é o GENE?

### Resposta:
O **gene** não é uma estrutura atômica única no código, mas sim **fragmentos moleculares** que são manipulados durante os operadores genéticos:

- **No Crossover:** Os genes são **fragmentos BRICS** (Breaking of Retrosynthetically Interesting Chemical Substructures)
- **Na Mutação:** Os genes são **padrões SMARTS** (grupos funcionais específicos)

### Implementação no código:

#### Durante o Crossover (Recombinação):
**Arquivo:** `scripts/generate_inhibitors_from_protein.py`  
**Linhas:** 415-416

```python
frags_a = list(BRICS.BRICSDecompose(mol_a))  # Fragmentos do pai A = "genes"
frags_b = list(BRICS.BRICSDecompose(mol_b))  # Fragmentos do pai B = "genes"
```

Cada fragmento BRICS representa uma unidade molecular sinteticamente válida (ex: anel aromático, grupo carboxila, cadeia lateral, etc.)

#### Durante a Mutação:
**Arquivo:** `scripts/generate_inhibitors_from_protein.py`  
**Linhas:** 381-392

```python
TRANSFORMS = [
    ("[c:1][F]", "[c:1][Cl]"),      # Gene: átomo de flúor em aromático
    ("[c:1][Cl]", "[c:1][F]"),      # Gene: átomo de cloro em aromático
    ("[c:1][Br]", "[c:1][Cl]"),     # Gene: átomo de bromo em aromático
    ("[c:1][CH3]", "[c:1][OH]"),    # Gene: grupo metil
    ("[c:1][OH]", "[c:1][NH2]"),    # Gene: grupo hidroxila
    ("[c:1][NH2]", "[c:1][CH3]"),   # Gene: grupo amino
    ("[C:1](=O)[OH]", "[C:1](=O)[NH2]"),  # Gene: ácido carboxílico
    # ... outros
]
```

Cada padrão SMARTS representa um "gene" (grupo funcional) que pode ser substituído.

---

## 2. O que é o CROMOSSOMO?

### Resposta:
O **cromossomo** é a **molécula completa** representada como uma string SMILES (Simplified Molecular Input Line Entry System).

### Estrutura no código:
**Arquivo:** `scripts/generate_inhibitors_from_protein.py`  
**Linhas:** 462-466

```python
population: list[tuple[float, str, str]] = []
#                      ↑      ↑     ↑
#                  fitness  SMILES origin
for name, mol in seeds:
    smi = mol_to_smiles(mol)  # CROMOSSOMO = SMILES string
    fit = fitness_protein_based(mol, site_info, ref_fps, pharmacophore)
    population.append((fit, smi, f"seed:{name}"))
```

### Composição de cada cromossomo:
- **fitness (float):** Aptidão da molécula (0.0 a 1.0)
- **smiles (str):** Representação canônica da molécula — **ESTE É O CROMOSSOMO**
- **origin (str):** Rastreabilidade (ex: "seed:Clavulanato", "mut_gen15", "cross_gen23")

### Exemplo de cromossomos:
```python
# Cromossomo 1 (Clavulanato):
"O=C(O)[C@H]1N2C(=C(CO)CO[C@@H]12)C(=O)O"

# Cromossomo 2 (Sulbactam):
"CC1(C)S[C@@H]2[C@H](NC(=O)C(=O)O)C(=O)N2[C@H]1C(=O)O"
```

---

## 3. Qual CROSSOVER é usado?

### Resposta:
**BRICS Crossover** — Recombinação baseada em fragmentação sintética.

### Implementação completa:
**Arquivo:** `scripts/generate_inhibitors_from_protein.py`  
**Linhas:** 410-436

```python
def crossover(mol_a: Chem.Mol, mol_b: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    """Crossover via BRICS."""
    from rdkit.Chem import BRICS
    
    try:
        # 1. Decompor cada pai em fragmentos BRICS
        frags_a = list(BRICS.BRICSDecompose(mol_a))
        frags_b = list(BRICS.BRICSDecompose(mol_b))
        
        if not frags_a or not frags_b:
            return None
        
        # 2. Selecionar um fragmento aleatório do pai B
        frag_from_b = rng.choice(frags_b)
        
        # 3. Combinar: metade dos fragmentos de A + 1 fragmento de B
        combined = frags_a[: len(frags_a) // 2] + [frag_from_b]
        
        # 4. Converter fragmentos SMILES para moléculas RDKit
        frag_mols = [Chem.MolFromSmiles(f) for f in combined if Chem.MolFromSmiles(f)]
        if not frag_mols:
            return None
        
        # 5. Reconstruir molécula via BRICS.BRICSBuild
        built = list(BRICS.BRICSBuild(frag_mols))
        if not built:
            return None
        
        # 6. Selecionar um candidato aleatório e sanitizar
        candidate = rng.choice(built[:5])
        Chem.SanitizeMol(candidate)
        return candidate
    except:
        return None
```

### Características do BRICS Crossover:
- **Não uniforme:** Corta o pai A pela metade + adiciona 1 fragmento do pai B
- **Quimicamente válido:** Apenas pontos de quebra sinteticamente viáveis (ligações estratégicas)
- **Tamanho variável:** Filho pode ter tamanho diferente dos pais
- **Diversidade estrutural:** BRICS pode gerar múltiplas recombinações

### Exemplo visual:
```
Pai A: [Fragmento1]--[Fragmento2]--[Fragmento3]--[Fragmento4]
Pai B: [FragmentoX]--[FragmentoY]--[FragmentoZ]

Crossover:
1. Pega metade de A: [Fragmento1]--[Fragmento2]
2. Pega 1 de B:      [FragmentoZ]
3. Recombina:        [Fragmento1]--[Fragmento2]--[FragmentoZ]

Filho: Nova molécula sinteticamente válida
```

---

## 4. Qual MECANISMO é usado no Algoritmo Genético?

### Resposta:
**Steady-State GA com Elitismo** (Algoritmo Genético de Estado Estacionário)

### Implementação do laço principal:
**Arquivo:** `scripts/generate_inhibitors_from_protein.py`  
**Linhas:** 474-531

```python
for gen in range(1, generations + 1):
    new_individuals: list[tuple[float, str, str]] = []
    
    # 1. SELEÇÃO ELITISTA
    population.sort(key=lambda x: x[0], reverse=True)  # Ordena por fitness
    elite_size = max(2, population_size // 5)          # Top 20% (mínimo 2)
    elite = population[:elite_size]
    
    attempts = 0
    max_attempts = population_size * 4
    
    # 2. GERAÇÃO DE NOVOS INDIVÍDUOS
    while len(new_individuals) < population_size and attempts < max_attempts:
        attempts += 1
        
        # 3. ESCOLHA DE OPERADOR (baseada na taxa de mutação)
        if rng.random() < mutation_rate:  # Ex: 50% de chance
            # MUTAÇÃO: seleciona 1 pai da elite
            _, smi, _ = rng.choice(elite)
            parent = Chem.MolFromSmiles(smi)
            if parent is None:
                continue
            child = mutate(parent, rng)
            origin_tag = f"mut_gen{gen}"
        else:
            # CROSSOVER: seleciona 2 pais da elite
            if len(elite) < 2:
                continue
            (_, smi_a, _), (_, smi_b, _) = rng.sample(elite, 2)
            mol_a = Chem.MolFromSmiles(smi_a)
            mol_b = Chem.MolFromSmiles(smi_b)
            if mol_a is None or mol_b is None:
                continue
            child = crossover(mol_a, mol_b, rng)
            origin_tag = f"cross_gen{gen}"
        
        if child is None:
            continue
        
        # 4. AVALIAÇÃO DE FITNESS
        child_smi = mol_to_smiles(child)
        if child_smi in best_seen:  # Evita duplicatas
            continue
        
        fit = fitness_protein_based(child, site_info, ref_fps, pharmacophore)
        best_seen[child_smi] = fit
        new_individuals.append((fit, child_smi, origin_tag))
        
        if fit > best_ever[0]:
            best_ever = (fit, child_smi, origin_tag)
    
    # 5. SUBSTITUIÇÃO: Elite sempre sobrevive + novos filhos
    population = elite + new_individuals
    population.sort(key=lambda x: x[0], reverse=True)
    population = population[:population_size]  # Mantém tamanho fixo
```

### Componentes do Mecanismo:

| Componente | Implementação | Localização |
|------------|---------------|-------------|
| **Seleção** | Elitista — preserva top 20% | Linha 477-479 |
| **Operadores** | Mutação (50%) vs Crossover (50%) baseado em `mutation_rate` | Linha 487-503 |
| **Seleção de pais** | Aleatória da elite (não-determinística) | Linha 488, 497 |
| **Substituição** | Geracional com elitismo — elite sempre sobrevive | Linha 519-521 |
| **Controle de diversidade** | Hash de SMILES em `best_seen` evita duplicatas | Linha 471, 509 |
| **Tamanho fixo** | População sempre mantida em `population_size` | Linha 521 |

### Características do mecanismo:
1. **Elitismo forte:** Melhores 20% sempre sobrevivem
2. **Pressão seletiva alta:** Apenas elite pode reproduzir
3. **Diversidade garantida:** Não aceita SMILES duplicados
4. **Convergência controlada:** Tamanho fixo de população

---

## 5. Os cromossomos são de tamanhos iguais?

### Resposta:
**NÃO, os cromossomos NÃO são de tamanhos iguais.**

### Por quê?

1. **Moléculas químicas têm comprimentos variáveis** (número de átomos diferente)
2. **BRICS crossover gera moléculas de tamanho variável** (linha 422)
3. **Mutação mantém aproximadamente o tamanho**, mas SMILES tem comprimento variável
4. **Não há padding, truncamento ou normalização de tamanho** em nenhuma parte do código

### Evidências no código:

#### 1. Cromossomos são SMILES strings de comprimento arbitrário:
**Linhas:** 462-466
```python
population: list[tuple[float, str, str]] = []
# Nenhuma restrição de tamanho na string SMILES
for name, mol in seeds:
    smi = mol_to_smiles(mol)  # Pode ter qualquer comprimento
```

#### 2. Seeds iniciais têm tamanhos muito diferentes:
**Linhas:** 66-73
```python
KNOWN_INHIBITORS = {
    "Clavulanato": "O=C(O)[C@H]1N2C(=C(CO)CO[C@@H]12)C(=O)O",  # 39 caracteres
    "Sulbactam": "CC1(C)S[C@@H]2[C@H](NC(=O)C(=O)O)C(=O)N2[C@H]1C(=O)O",  # 56 caracteres
    "Tazobactam": "CN1C(=O)N2[C@H]([C@H](C)S(=O)(=O)N(C)C2(C)C)C1C(=O)O",  # 55 caracteres
    "Avibactam": "C[C@H]1CN(C(=O)N(O1)[C@@H]2CNC(=O)N2)S(=O)(=O)N",  # 50 caracteres
    "Relebactam": "O=C(O)C1NC(=O)N[C@H]1c1cnccn1",  # 33 caracteres
    "Vaborbactam": "CC(C)(C)C(=O)N[C@@H]1B(O)OCC1(C)C",  # 38 caracteres
}
```

#### 3. Crossover gera moléculas de tamanho variável:
**Linha:** 422
```python
combined = frags_a[: len(frags_a) // 2] + [frag_from_b]
# Tamanho do filho = ~metade de A + 1 fragmento de B
# Se A tem 4 fragmentos e B tem 3 fragmentos:
#   Filho = 2 fragmentos de A + 1 fragmento de B = 3 fragmentos
# Tamanho final depende do tamanho de cada fragmento individual
```

### Exemplo de variabilidade:
```python
# Pai A (Clavulanato): 15 átomos pesados
# Pai B (Sulbactam):   17 átomos pesados
# Filho (crossover):   pode ter 10-20 átomos pesados (variável!)
```

---

## 6. Se não são iguais, como deixar (como lidar com isso)?

### Resposta:
O código **NÃO normaliza os tamanhos** e trata a variabilidade de forma **natural e biológica**.

### Estratégias implementadas:

### 6.1. Penalização por tamanho na função de fitness
**Arquivo:** `scripts/generate_inhibitors_from_protein.py`  
**Linhas:** 273-288

```python
def fitness_protein_based(mol, site_info, ref_fps, pharmacophore):
    # 1. Compatibilidade de tamanho (30% da fitness total)
    mw = Descriptors.MolWt(mol)
    volume_score = 0.0
    
    # Volume ideal baseado no sítio ativo da proteína
    site_volume = site_info["volume"]
    ideal_mw_min = max(200, site_volume / 5)  # Peso molecular mínimo
    ideal_mw_max = min(600, site_volume / 2)  # Peso molecular máximo
    
    if ideal_mw_min <= mw <= ideal_mw_max:
        volume_score = 1.0  # Tamanho ideal = score máximo
    elif mw < ideal_mw_min:
        volume_score = mw / ideal_mw_min  # Penaliza moléculas pequenas
    else:
        volume_score = max(0.0, 1.0 - (mw - ideal_mw_max) / 200)  # Penaliza moléculas grandes
```

**Explicação:**
- Moléculas muito pequenas (< 200 Da) recebem fitness reduzido
- Moléculas muito grandes (> 600 Da) recebem fitness reduzido
- Faixa ideal é calculada dinamicamente baseada no **volume do sítio ativo**

### 6.2. Crossover BRICS aceita qualquer tamanho
**Linhas:** 415-436

O operador BRICS **naturalmente gera** moléculas de tamanhos variados:
```python
# Não há verificação de tamanho
# Qualquer recombinação válida é aceita
frags_a = list(BRICS.BRICSDecompose(mol_a))  # 2-10 fragmentos
frags_b = list(BRICS.BRICSDecompose(mol_b))  # 2-8 fragmentos
combined = frags_a[: len(frags_a) // 2] + [frag_from_b]  # 2-6 fragmentos
```

### 6.3. Mutação preserva aproximadamente o tamanho
**Linhas:** 379-407

Substituições SMARTS são **isostéricas** (grupos de tamanho similar):
```python
TRANSFORMS = [
    ("[c:1][F]", "[c:1][Cl]"),      # F → Cl (1 átomo → 1 átomo)
    ("[c:1][CH3]", "[c:1][OH]"),    # CH3 → OH (pequena variação)
    ("[C:1](=O)[OH]", "[C:1](=O)[NH2]"),  # COOH → CONH2 (mesma estrutura base)
]
```

### 6.4. Sanitização química garante validade
**Linhas:** 271, 402, 433

Independente do tamanho, todas as moléculas passam por validação química:
```python
Chem.SanitizeMol(mol)  # Verifica valências, aromaticidade, conectividade
```

Se a molécula é quimicamente inválida (independente do tamanho), ela é **descartada**.

---

## Por que NÃO normalizar tamanhos?

### Justificativa biológica e química:

1. **Algoritmos genéticos naturais não têm cromossomos de tamanho fixo**
   - Genomas de diferentes espécies têm tamanhos variados
   - Recombinação natural gera variabilidade de tamanho

2. **Moléculas químicas NÃO devem ter tamanho forçado**
   - Inibidores reais têm pesos moleculares de 150-800 Da
   - Padding artificial criaria estruturas quimicamente inválidas
   - Truncamento removeria grupos funcionais essenciais

3. **Fitness baseada em tamanho é mais elegante**
   - Moléculas inadequadas recebem score baixo naturalmente
   - Pressão seletiva remove moléculas muito grandes/pequenas
   - Não introduz viés artificial

4. **SMILES é uma representação variável por natureza**
   - "CCO" (etanol) = 3 caracteres
   - "CC(C)(C)C(=O)O" (ácido pivalico) = 15 caracteres
   - Forçar mesmo comprimento quebraria a representação química

---

## Resumo para o Orientador

| Pergunta | Resposta | Arquivo:Linha |
|----------|----------|---------------|
| **Gene** | Fragmentos BRICS (crossover) ou padrões SMARTS (mutação) | `generate_inhibitors_from_protein.py:415-416, 381-392` |
| **Cromossomo** | Molécula completa como string SMILES | `generate_inhibitors_from_protein.py:462-466` |
| **Crossover** | BRICS Crossover (fragmentação sintética) | `generate_inhibitors_from_protein.py:410-436` |
| **Mecanismo** | Steady-State GA com Elitismo | `generate_inhibitors_from_protein.py:474-531` |
| **Tamanhos iguais?** | NÃO — cromossomos têm tamanhos variáveis | Todo o código (nenhuma normalização) |
| **Como lidar?** | Penalização de tamanho na fitness + seleção natural | `generate_inhibitors_from_protein.py:273-288` |

---

## Validação da Abordagem

### Resultados experimentais (1AXB.pdb):
**Arquivo de saída:** `results/inhibitors_1AXB.csv`

```
Top 10 candidatos gerados:
rank  fitness    mw       logp    qed
1     0.6443     305.36   -0.42   0.571
2     0.6412     289.28   0.15    0.623
3     0.6328     312.33   -0.18   0.598
...
```

**Observações:**
- Pesos moleculares variam de 289-340 Da (tamanhos diferentes!)
- Todos têm fitness > 0.60 (alta qualidade)
- Diversidade estrutural preservada
- QED > 0.5 (drug-likeness aceitável)

**Conclusão:** A variabilidade de tamanho **NÃO prejudica** a qualidade dos candidatos, pelo contrário, permite explorar melhor o espaço químico.

---

## Referências Técnicas

1. **BRICS:** Degen, J., Wegscheid-Gerlach, C., Zaliani, A., & Rarey, M. (2008). *On the Art of Compiling and Using 'Drug-Like' Chemical Fragment Spaces.* ChemMedChem, 3(10), 1503-1507.

2. **Steady-State GA:** Whitley, D. (1989). *The GENITOR Algorithm and Selection Pressure: Why Rank-Based Allocation of Reproductive Trials is Best.* ICGA, 89, 116-121.

3. **RDKit BRICS:** https://www.rdkit.org/docs/source/rdkit.Chem.BRICS.html

4. **SMARTS:** https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

---

**Documento gerado em:** 30 de março de 2026  
**Versão do código:** `generate_inhibitors_from_protein.py` (última versão)  
**Autor:** Milena (Dissertação de Mestrado)
