import pytest # noqa
from app.parser import parse_molecule

def test_valid_smiles():
    mol = parse_molecule("CCO")  # Ethanol
    assert mol is not None
    assert mol.GetNumAtoms() > 0

def test_invalid_smiles():
    mol = parse_molecule("XYZ123")
    assert mol is None
