from rdkit import Chem
from rdkit.Chem import AllChem

def parse_molecule(smiles: str):
    """
    Parses a SMILES string into an RDKit molecule object with 3D coordinates.
    Returns None if parsing fails.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Add hydrogens for proper geometry
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates using ETKDGv3 algorithm
        params = AllChem.ETKDGv3()
        result = AllChem.EmbedMolecule(mol, params)
        if result != 0:
            return None

        # Optimize geometry (optional but improves bond lengths)
        AllChem.UFFOptimizeMolecule(mol)

        return mol

    except Exception as e:
        print(f"Error parsing molecule: {e}")
        return None
