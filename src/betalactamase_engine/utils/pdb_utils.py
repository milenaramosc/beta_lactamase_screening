import warnings
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from pathlib import Path


# Silence PDB construction warnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

def get_pdb_parser():
    """Returns a Bio.PDB parser instance."""
    return PDBParser(QUIET=True)

def load_structure(pdb_path: Path):
    """Loads a structure from a PDB file using Bio.PDB."""
    parser = get_pdb_parser()
    structure_id = pdb_path.stem
    return parser.get_structure(structure_id, str(pdb_path))

def calculate_centroid(atoms):
    """Calculates the geometric center (centroid) of a list of atoms."""
    if not atoms:
        return None
    
    x_coords = [atom.get_coord()[0] for atom in atoms]
    y_coords = [atom.get_coord()[1] for atom in atoms]
    z_coords = [atom.get_coord()[2] for atom in atoms]
    
    centroid = [
        float(sum(x_coords) / len(atoms)),
        float(sum(y_coords) / len(atoms)),
        float(sum(z_coords) / len(atoms))
    ]
    return centroid
