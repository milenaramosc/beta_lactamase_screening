from rdkit import Chem
from rdkit.Chem import Draw

# Carregar o arquivo PDB
# mol = Chem.MolFromPDBFile('/home/milena/beta_lactamase_screening/results/complexes/6C78_ZINC000000033518.pdb')

# # Adicionar hidrogênios
# mol = Chem.AddHs(mol)

# # Gerar SMILES
# smiles = Chem.MolToSmiles(mol)

# Exemplo: SMILES de cafeína
smiles = "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C.C1C2N(C1=O)C(C(=CCO)O2)C(=O)[O-].[K+]"

# Converta SMILES em molécula (objeto RDKit)
mol = Chem.MolFromSmiles(smiles)

# Gere uma imagem 2D da molécula
img = Draw.MolToImage(mol, size=(800, 600))

# Salve ou exiba a imagem
# img.save("molecula_2d.png")
img.show()
# print(smiles)   