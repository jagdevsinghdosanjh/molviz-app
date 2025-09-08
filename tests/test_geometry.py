import pytest # noqa
from rdkit import Chem
from app.geometry import compute_geometry

def parse_molecule(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    # Continue with embedding and optimization...


def test_geometry_output_structure():
    mol = Chem.MolFromSmiles("CCO")
    if mol is None:
        pytest.fail("Failed to parse SMILES string.")
    mol = Chem.AddHs(mol)
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)

    geometry = compute_geometry(mol)
    assert isinstance(geometry, dict)
    assert "geometry" in geometry
    assert "angles" in geometry
    assert "lengths" in geometry
    assert all(isinstance(angle, float) for angle in geometry["angles"])
    assert all(isinstance(length, float) for length in geometry["lengths"])
