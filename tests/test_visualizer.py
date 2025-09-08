import pytest # noqa
from rdkit import Chem
from app.visualizer import render_3d

def test_render_3d_runs_without_error():
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None, "Chem.MolFromSmiles returned None"
    mol = Chem.AddHs(mol)
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)

    # This test checks that render_3d executes without crashing
    try:
        render_3d(mol)
    except Exception as e:
        pytest.fail(f"render_3d failed with error: {e}")
