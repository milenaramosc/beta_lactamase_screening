Análise Completa do Algoritmo Genético e Tratamento de PDBs
1. ALGORITMO GENÉTICO - Estrutura e Componentes
1.1 Gene (Unidade Básica de Informação)
Definição: O "gene" não é explicitamente uma estrutura atômica no código, mas sim:
- Fragmentos BRICS durante o crossover
- Padrões SMARTS durante a mutação
Localização no código:
# scripts/generate_candidate.py:206-238
def crossover(mol_a: Chem.Mol, mol_b: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    from rdkit.Chem import BRICS
    
    frags_a = list(BRICS.BRICSDecompose(mol_a))  # Fragmentos = "genes"
    frags_b = list(BRICS.BRICSDecompose(mol_b))
Por que BRICS? BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures) quebra moléculas em fragmentos quimicamente válidos que podem ser recombinados mantendo a viabilidade sintética.
---
1.2 Cromossomo (Solução Completa)
Definição: Cada cromossomo é uma molécula inteira representada como string SMILES.
Localização no código:
# scripts/generate_candidate.py:262-266
population: list[tuple[float, str, str]] = []
#                      ↑      ↑     ↑
#                  fitness  SMILES origin
for name, mol, dg in seeds:
    smi = mol_to_smiles(mol)  # Cromossomo = SMILES
    fit = fitness(mol, ref_fps)
    population.append((fit, smi, f"seed:{name}"))
Estrutura: (fitness: float, smiles: str, origin: str)
- fitness: valor de aptidão da molécula
- smiles: representação canônica da molécula (o cromossomo propriamente dito)
- origin: rastreabilidade (seed, mutate, crossover)
---
1.3 Função de Fitness
Localização: scripts/generate_candidate.py:131-144
def fitness(mol: Chem.Mol, ref_fps: list) -> float:
    """
    Função de fitness combinada:
      - Similaridade com compostos de referência (alto = estrutura promissora)
      - Penalidade por violações de Lipinski
      - Penalidade por tamanho excessivo
    """
    lip   = lipinski_score(mol)        # Regra de Lipinski (druglikeness)
    sim   = tanimoto_to_reference(mol, ref_fps)  # Similaridade estrutural
    mw    = Descriptors.MolWt(mol)
    size_pen = max(0.0, (mw - 400) / 200)  # Penaliza MW > 400
    
    return sim * lip - size_pen * 0.1
Componentes:
1. Lipinski Score (generate_candidate.py:103-116): penaliza violações das regras de Lipinski (MW>500, LogP>5, HBD>5, HBA>10) - 0.25 por violação
2. Tanimoto Similarity (generate_candidate.py:119-128): calcula similaridade máxima com os compostos seed usando Morgan fingerprints (raio 2, 2048 bits)
3. Size Penalty: penaliza moléculas muito grandes
---
1.4 Operador de Crossover
Tipo: BRICS Crossover (recombinação de fragmentos sintéticos)
Localização: scripts/generate_candidate.py:206-238
def crossover(mol_a: Chem.Mol, mol_b: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    from rdkit.Chem import BRICS
    
    try:
        # 1. Decompor cada pai em fragmentos BRICS
        frags_a = list(BRICS.BRICSDecompose(mol_a))
        frags_b = list(BRICS.BRICSDecompose(mol_b))
        
        if not frags_a or not frags_b:
            return None
        
        # 2. Selecionar fragmento aleatório do pai B
        frag_from_b = rng.choice(frags_b)
        
        # 3. Combinar: metade dos fragmentos de A + 1 fragmento de B
        combined = frags_a[: len(frags_a) // 2] + [frag_from_b]
        
        # 4. Converter fragmentos SMILES para moléculas
        frag_mols = [Chem.MolFromSmiles(f) for f in combined if Chem.MolFromSmiles(f)]
        
        # 5. Reconstruir molécula via BRICS.BRICSBuild
        built = list(BRICS.BRICSBuild(frag_mols))
        if not built:
            return None
        
        # 6. Selecionar um candidato e sanitizar
        candidate = rng.choice(built[:5])
        Chem.SanitizeMol(candidate)
        return candidate
    except Exception:
        return None
Características:
- Não é uniforme: corta o pai A pela metade e adiciona 1 fragmento de B
- Respeita química: fragmentos BRICS são sinteticamente válidos
- Tamanho variável: filho pode ter tamanho diferente dos pais
---
1.5 Operador de Mutação
Tipo: SMARTS Pattern Substitution (substituição de grupos funcionais)
Localização: scripts/generate_candidate.py:167-203
def mutate(mol: Chem.Mol, rng: random.Random) -> Chem.Mol | None:
    from rdkit.Chem import AllChem
    
    # Lista de transformações químicas válidas
    TRANSFORMS = [
        ("[c:1][F]", "[c:1][Cl]"),      # F → Cl em aromático
        ("[c:1][Cl]", "[c:1][F]"),      # Cl → F em aromático
        ("[c:1][Br]", "[c:1][Cl]"),     # Br → Cl em aromático
        ("[c:1][CH3]", "[c:1][OH]"),    # Metil → Hidroxil
        ("[c:1][OH]", "[c:1][NH2]"),    # Hidroxil → Amina
        ("[c:1][NH2]", "[c:1][CH3]"),   # Amina → Metil
        ("[C:1](=O)[OH]", "[C:1](=O)[NH2]"),  # Ácido → Amida
    ]
    
    rng.shuffle(TRANSFORMS)
    
    for reactant_smarts, product_smarts in TRANSFORMS:
        try:
            rxn = AllChem.ReactionFromSmarts(f"{reactant_smarts}>>{product_smarts}")
            products = rxn.RunReactants((mol,))
            if products:
                product = products[0][0]
                Chem.SanitizeMol(product)
                return product
        except Exception:
            continue
    
    return None
Características:
- Local: altera apenas um grupo funcional
- Química válida: preserva valências e aromaticidade
- Tamanho quase constante: substituições isostéricas ou pequenas variações
---
1.6 Mecanismo do Algoritmo Genético
Tipo: Steady-State GA com Elitismo
Laço principal: scripts/generate_candidate.py:272-327
for gen in range(1, generations + 1):
    new_individuals: list[tuple[float, str, str]] = []
    
    # 1. SELEÇÃO: Ordenar por fitness e selecionar elite
    population.sort(key=lambda x: x[0], reverse=True)
    elite = population[: max(2, population_size // 5)]  # Top 20% (mínimo 2)
    
    attempts = 0
    while len(new_individuals) < population_size and attempts < population_size * 4:
        attempts += 1
        
        # 2. ESCOLHA DE OPERADOR: 50% mutação, 50% crossover
        op = rng.choice(["mutate", "crossover"])
        
        try:
            if op == "mutate" and elite:
                # Mutação: seleciona 1 pai da elite aleatoriamente
                _, smi, origin = rng.choice(elite)
                parent = Chem.MolFromSmiles(smi)
                if parent is None:
                    continue
                child = mutate(parent, rng)
                origin_tag = f"mutate(gen{gen})"
            
            else:  # crossover
                # Crossover: seleciona 2 pais da elite sem reposição
                if len(elite) < 2:
                    continue
                (_, smi_a, _), (_, smi_b, _) = rng.sample(elite, 2)
                mol_a = Chem.MolFromSmiles(smi_a)
                mol_b = Chem.MolFromSmiles(smi_b)
                if mol_a is None or mol_b is None:
                    continue
                child = crossover(mol_a, mol_b, rng)
                origin_tag = f"crossover(gen{gen})"
            
            if child is None:
                continue
            
            # 3. AVALIAÇÃO DE FITNESS
            child_smi = mol_to_smiles(child)
            if child_smi in best_seen:
                continue  # Evita duplicatas
            
            fit = fitness(child, ref_fps)
            best_seen[child_smi] = fit
            new_individuals.append((fit, child_smi, origin_tag))
        
        except Exception:
            continue
    
    # 4. SUBSTITUIÇÃO: Elite + novos filhos, corta para população fixa
    population = elite + new_individuals
    population.sort(key=lambda x: x[0], reverse=True)
    population = population[:population_size]  # Mantém tamanho fixo
Componentes do mecanismo:
1. Seleção: Elitista - preserva top 20% (função generate_candidate.py:276-277)
2. Operadores: Probabilidade implícita 50/50 mutação/crossover (generate_candidate.py:282)
3. Seleção de pais: Aleatória da elite (não-determinística)
4. Substituição: Geracional com elitismo - elite sempre sobrevive
5. Controle de diversidade: Hash de SMILES em best_seen evita duplicatas (generate_candidate.py:270, 308)
6. Convergência: Tamanho de população fixo (generate_candidate.py:320)
---
1.7 Tamanho dos Cromossomos
Resposta: NÃO SÃO DE TAMANHOS IGUAIS
Por quê?
- Moléculas químicas têm comprimentos variáveis (número de átomos/fragmentos)
- BRICS crossover pode gerar moléculas maiores ou menores que os pais
- Mutação mantém aproximadamente o tamanho, mas SMILES tem comprimento variável
Evidência no código:
# Nenhum padding, truncamento ou normalização de tamanho
# Cromossomos são SMILES strings de comprimento arbitrário
# Exemplo de variabilidade:
# Pai A: "CCO" (3 átomos, etanol)
# Pai B: "c1ccccc1C(=O)O" (benzoato)
# Filho: pode ter 5-15 átomos dependendo dos fragmentos
Como o GA lida com tamanhos variáveis:
1. Penalização por tamanho no fitness:
# generate_candidate.py:141
size_pen = max(0.0, (mw - 400) / 200)  # Penaliza MW > 400
return sim * lip - size_pen * 0.1
2. Crossover BRICS aceita qualquer tamanho:
# generate_candidate.py:222
combined = frags_a[: len(frags_a) // 2] + [frag_from_b]
# Tamanho do filho = ~metade de A + 1 fragmento de B
3. Não há normalização de comprimento: O código aceita moléculas de qualquer tamanho válido, desde que passem pela sanitização química.
---
2. TRATAMENTO DE ARQUIVOS PDB
2.1 Por que os PDBs são tratados?
Resposta: Para preparar a proteína para docking molecular com AutoDock Vina, que requer:
1. Estrutura limpa (sem água/ligantes co-cristalizados)
2. Hidrogênios adicionados (não presentes em cristalografia de raios-X)
3. Formato PDBQT com cargas parciais e tipos atômicos
2.2 Pipeline de Tratamento
Etapa 1: Preparação da Proteína
Script: scripts/prepare_protein.py
Passo 1.1: Remoção de heteroátomos (prepare_protein.py:71-88)
def strip_heteroatoms(pdb_path: Path, out_path: Path) -> int:
    """
    Remove HETATM records (água + ligantes co-cristalizados).
    Mantém ATOM (cadeia proteica) e registros de cabeçalho/SITE.
    """
    kept = []
    removed = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):  # Água, ligantes, íons
                removed += 1
            else:
                kept.append(line)  # ATOM, SITE, HEADER, etc.
    with open(out_path, "w") as f:
        f.writelines(kept)
    return removed
Por quê? 
- Água cristalográfica (HOH) interfere no docking
- Ligantes co-cristalizados não fazem parte do receptor para novo docking
- Mantém SITE records para detectar o sítio ativo
---
Passo 1.2: Adição de hidrogênios (prepare_protein.py:92-110)
def add_hydrogens(pdb_path: Path, out_path: Path) -> bool:
    """
    Adiciona hidrogênios via OpenBabel e salva o PDB resultante.
    """
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdb")
    
    mol = ob.OBMol()
    if not conv.ReadFile(mol, str(pdb_path)):
        raise RuntimeError(f"OpenBabel não conseguiu ler {pdb_path}")
    
    mol.AddHydrogens()  # Adiciona H polares e apolares
    conv.WriteFile(mol, str(out_path))
    return True
Por quê?
- Cristalografia de raios-X não detecta hidrogênios (muito leves)
- AutoDock Vina precisa de H para calcular ligações de hidrogênio
- OpenBabel posiciona H em geometrias quimicamente corretas
---
Passo 1.3: Detecção do sítio ativo (prepare_protein.py:140-192)
def detect_binding_site(pdb_path: Path) -> dict:
    """
    Detecta sítio ativo e retorna dict com center_x/y/z e size_x/y/z.
    Prioridade:
      1. Registros SITE do PDB
      2. Fallback: centroide de todos os CA + caixa 30 Å
    """
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("X", str(pdb_path))
    model  = struct[0]
    
    # 1. Tentar SITE records
    site_residues = parse_site_records(pdb_path)
    coords = []
    
    if site_residues:
        for resname, chain, resnum in site_residues:
            try:
                res = model[chain][resnum]
                ca  = res["CA"]
                coords.append(ca.get_vector().get_array())
            except (KeyError, Exception):
                pass
    
    if coords:
        coords  = np.array(coords)
        center  = coords.mean(axis=0)  # Centróide 3D
        extents = coords.max(axis=0) - coords.min(axis=0)
        box     = np.maximum(extents + 10, 20.0)  # +10 Å padding, mínimo 20 Å
        method  = "SITE records"
    else:
        # 2. Fallback: centróide de todos os CA
        for chain in model:
            for res in chain:
                if "CA" in res:
                    coords.append(res["CA"].get_vector().get_array())
        coords = np.array(coords)
        center = coords.mean(axis=0)
        box    = np.array([30.0, 30.0, 30.0])
        method = "centroide CA (fallback)"
    
    return {
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x":   round(float(box[0]), 3),
        "size_y":   round(float(box[1]), 3),
        "size_z":   round(float(box[2]), 3),
        "method":   method,
    }
Por quê?
- AutoDock Vina precisa de uma search box (caixa de busca) para docking
- SITE records do PDB indicam resíduos funcionalmente importantes
- Fallback garante que sempre há um sítio definido
Saída: config.yaml atualizado com coordenadas (prepare_protein.py:329-337)
---
Etapa 2: Conversão PDB → PDBQT
Script: scripts/screen_and_rank.py
Passo 2.1: Receptor PDBQT (screen_and_rank.py:170-185)
def prepare_receptor_pdbqt(pdb_path: Path, out_dir: Path) -> Path:
    """Convert protein PDB → PDBQT (rigid receptor) using OpenBabel."""
    out_path = out_dir / f"{pdb_path.stem}_receptor.pdbqt"
    if out_path.exists():
        return out_path
    
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdbqt")
    conv.AddOption("r", ob.OBConversion.OUTOPTIONS)  # rigid (receptor rígido)
    mol = ob.OBMol()
    ok = conv.ReadFile(mol, str(pdb_path))
    if not ok or mol.NumAtoms() == 0:
        raise RuntimeError(f"OpenBabel could not read {pdb_path}")
    conv.WriteFile(mol, str(out_path))
    return out_path
Por quê?
- AutoDock Vina requer formato PDBQT (PDB + Qt = cargas + tipos atômicos)
- OpenBabel adiciona automaticamente cargas Gasteiger e tipos AD4
- Flag "r" marca o receptor como rígido (sem rotâmeros internos)
---
Passo 2.2: Ligante PDBQT (screen_and_rank.py:202-225)
def prepare_ligand_pdbqt(mol, name: str, out_dir: Path) -> Path | None:
    """Convert RDKit Mol (with 3D coords) → PDBQT using Meeko."""
    out_path = out_dir / f"{name}.pdbqt"
    if out_path.exists():
        return out_path
    if mol_has_unsupported_elements(mol):  # Verifica elementos suportados
        return None
    try:
        prep = MoleculePreparation()
        molsetups = prep.prepare(mol)  # Detecta ligações rotáveis
        if not molsetups:
            return None
        pdbqt_str, ok, msg = PDBQTWriterLegacy.write_string(molsetups[0])
        if not ok:
            return None
        out_path.write_text(pdbqt_str)
        return out_path
    except Exception:
        return None
Por quê?
- Meeko detecta ligações rotáveis (graus de liberdade) do ligante
- AutoDock Vina trata ligante como flexível (amostra conformações)
- Filtra elementos não suportados (B, Si, metais exóticos)
---
Etapa 3: Docking e Saída de Poses
Execução do Vina: (screen_and_rank.py:262-310)
def run_vina(
    vina_bin: str,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    out_pdbqt: Path,
    center_x: float, center_y: float, center_z: float,
    size_x: float, size_y: float, size_z: float,
) -> float | None:
    """Run AutoDock Vina and parse binding affinity from output."""
    cmd = [
        vina_bin,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--out", str(out_pdbqt),  # Grava poses em PDBQT
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--cpu", "1",
        "--exhaustiveness", "8",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    
    # Parse ΔG do output
    for line in result.stdout.split("\n"):
        if "   1 " in line:  # Primeiro modo (melhor afinidade)
            match = re.search(r"(-?\d+\.\d+)", line)
            if match:
                return float(match.group(1))
    return None
Saída: 
- results/poses/<target>__<compound>.pdbqt com até 9 conformações ranqueadas
- results/ranking.csv com ΔG (kcal/mol)
---
Etapa 4: Conversão PDBQT → PDB para Visualização
Script: scripts/save_and_visualize.py
Conversão PDBQT→PDB (save_and_visualize.py:47-72)
def pdbqt_to_pdb_lines(pdbqt_path: Path) -> list[str]:
    """
    Extract the first MODEL from a PDBQT file and convert to basic PDB lines.
    Strips AutoDock-specific columns (partial charges, atom type at end).
    """
    in_model = False
    pdb_lines = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                break  # only first model (melhor pose)
            if not in_model:
                continue
            if line.startswith(("ATOM", "HETATM")):
                # Trunca colunas extras do PDBQT (cargas, tipo atômico)
                pdb_line = line[:66].rstrip()
                # Renomeia resíduo para LIG
                pdb_line = pdb_line[:17] + "LIG" + pdb_line[20:]
                pdb_lines.append(pdb_line)
    return pdb_lines
Criação de complexo (save_and_visualize.py:77-111)
def save_complex(
    protein_pdb_path: Path,
    pose_pdbqt_path: Path,
    out_path: Path,
    compound_id: str,
) -> Path:
    """Combine protein PDB + ligand pose into a single complex PDB file."""
    # Lê proteína (apenas ATOM/TER/SEQRES)
    with open(protein_pdb_path) as f:
        protein_lines = [
            line.rstrip()
            for line in f
            if line.startswith(("ATOM", "TER", "SEQRES"))
        ]
    
    # Converte ligante PDBQT para PDB
    ligand_lines = pdbqt_to_pdb_lines(pose_pdbqt_path)
    
    # Combina em um único PDB
    with open(out_path, "w") as f:
        f.write(f"REMARK Complex: {protein_pdb_path.stem} + {compound_id}\n")
        for line in protein_lines:
            f.write(line + "\n")
        f.write("TER\n")
        f.write(f"REMARK Ligand: {compound_id}\n")
        for line in ligand_lines:
            # Marca ligante como HETATM na chain L
            if line.startswith("ATOM"):
                line = "HETATM" + line[6:21] + "L" + line[22:]
            f.write(line + "\n")
        f.write("END\n")
    
    return out_path
Saída: results/complexes/<target>_<compound>.pdb para visualização 3D
---
2.3 Resumo: Por que toda essa tratativa?
Etapa	Entrada	Saída	Por quê?
prepare_protein.py	data/structures/<ID>.pdb (bruto)	data/prepared/<ID>_prepared.pdb	Remove água/ligantes, adiciona H, detecta sítio ativo
screen_and_rank.py	PDB preparado	receptor.pdbqt + ligand.pdbqt	Formato exigido pelo AutoDock Vina (cargas + tipos atômicos)
AutoDock Vina	PDBQTs	results/poses/<ID>.pdbqt	Gera conformações docadas com ΔG
save_and_visualize.py	PDBQT poses	results/complexes/<ID>.pdb	Formato padrão para visualização (PyMOL, py3Dmol)
Fluxo completo:
PDB bruto → [limpa] → PDB preparado → [converte] → PDBQT → [docka] → PDBQT poses → [converte] → PDB complexo
Objetivo: Simular a ligação de moléculas candidatas à proteína beta-lactamase, ranqueando por afinidade de ligação (ΔG), para priorizar as melhores candidatas para síntese experimental.
---
3. Conexão GA ↔ PDB
O algoritmo genético (generate_candidate.py) NÃO manipula PDBs diretamente.
Fluxo integrado:
1. Docking inicial (screen_and_rank.py): PDBs de proteína + biblioteca química → ranking por ΔG
2. GA (generate_candidate.py): Top compostos do ranking → SMILES seed → evolução molecular → novos SMILES
3. Validação (manual ou próximo ciclo): Novos SMILES → docking contra PDBs → ranquear → realimentar GA
O GA opera em "espaço químico" (SMILES/fingerprints), não em "espaço estrutural" (coordenadas 3D/PDB).
▣  Build · claude-sonnet-4.5 · 2m 15s
tenho uma duvida no arquivo README.md as etapas de 3 a 6 são realmente necessaria para executar a etapa 7?
Thinking: O usuário está perguntando se as etapas 3 a 6 do README.md são realmente necessárias para executar a etapa 7. Preciso ler o README.md para ver quais são essas etapas e entender a dependência entre elas.
A etapa 7 é o generate_candidate.py que já analisamos. Vamos ver o que ele precisa como entrada.
→ Read README.md 
→ Read scripts/generate_candidate.py [limit=100, offset=59]
Resposta: NÃO, as etapas 3 a 6 NÃO são estritamente necessárias para executar a Etapa 7.
Análise das Dependências
O que a Etapa 7 (generate_candidate.py) realmente precisa?
Olhando o código em generate_candidate.py:59-98:
def load_seed_compounds(top_n: int) -> list[tuple[str, Chem.Mol, float]]:
    ranking_csv = PROJECT_DIR / "results" / "ranking.csv"      # ← Precisa disso
    sdf_path    = PROJECT_DIR / "data" / "compounds" / "compounds.sdf"  # ← E disso
Arquivos OBRIGATÓRIOS:
1. results/ranking.csv (gerado pela Etapa 4)
2. data/compounds/compounds.sdf (gerado pela Etapa 3)
---
Mapeamento Etapa → Arquivo
Etapa	Script	Saída
3	fetch_compounds.py	data/compounds/compounds.sdf
4	screen_and_rank.py	results/ranking.csv
5	admet_analysis.py	results/admet_report.csv, results/final_candidates.csv
6	save_and_visualize.py	results/complexes/, results/visualizations/
---
Conclusão
Etapas NECESSÁRIAS para rodar a Etapa 7:
- ✅ Etapa 1: fetch_structures.py (indiretamente, pois Etapa 4 precisa)
- ✅ Etapa 2: select_targets.py (indiretamente, pois Etapa 4 precisa)
- ✅ Etapa 2.5: prepare_protein.py (indiretamente, pois Etapa 4 precisa)
- ✅ Etapa 3: fetch_compounds.py → compounds.sdf
- ✅ Etapa 4: screen_and_rank.py → ranking.csv
- ❌ Etapa 5: admet_analysis.py → NÃO é necessária
- ❌ Etapa 6: save_and_visualize.py → NÃO é necessária
Sequência mínima:
# Etapas 1-2.5 (preparação de alvos)
python scripts/fetch_structures.py --query "beta-lactamase" --max 20
python scripts/select_targets.py
python scripts/prepare_protein.py
# Etapa 3 (compostos)
python scripts/fetch_compounds.py --source chembl --max 5000
# Etapa 4 (docking - gera o ranking)
python scripts/screen_and_rank.py
# Etapa 7 (GA - usa ranking.csv + compounds.sdf)
python scripts/generate_candidate.py
Etapas 5 e 6 são opcionais — servem apenas para análise ADMET e visualização dos resultados do docking, mas não são inputs do algoritmo genético.