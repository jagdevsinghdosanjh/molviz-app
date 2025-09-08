import streamlit as st
from rdkit import Chem  # noqa
from rdkit.Chem import AllChem  # noqa
from rdkit.Chem import rdDepictor

from app.ui import render_ui
from app.parser import parse_molecule
from app.geometry import compute_geometry
from app.visualizer import render_3d
from app.exporter import export_geometry_pdf

# Sample molecules
sample_smiles = {
    "Ethanol": "CCO",
    "Methane": "C",
    "Water": "O",
    "Ammonia": "N",
    "Benzene": "c1ccccc1",
    "BHC":"Cl[C@@H]1Cl[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H]1Cl",
    "Aniline":"c1ccccc1N",
    "Lactic Acid":"CC(C(=O)O)O",
    "Glycerine":"OCC(O)CO"
}

def main():
    st.set_page_config(page_title="Molecular Geometry Visualizer", layout="wide")
    st.title("🔬 Molecular Geometry Explorer")

    col1, col2 = st.columns([2, 1])
    with col1:
        selected = st.selectbox("Choose a sample molecule:", list(sample_smiles.keys()))
        smiles = sample_smiles[selected]
        st.text_input("SMILES (editable):", value=smiles, key="smiles_input")

    smiles = st.session_state.smiles_input

    view_mode = st.radio("View Mode:", ["3D", "2D"])

    mol = parse_molecule(smiles)
    if mol is None:
        st.error("Invalid SMILES string. Please try again.")
        return

    geometry_data = compute_geometry(mol)
    render_ui(geometry_data)

    if view_mode == "3D":
        render_3d(mol)
    else:
        rdDepictor.Compute2DCoords(mol)
        img = Chem.MolToImage(mol, size=(300, 300))
        st.image(img, caption="2D Structure")

    if st.button("📤 Export Geometry as PDF"):
        success = export_geometry_pdf("geometry_report.pdf", geometry_data, smiles)
        if success:
            st.success("PDF exported successfully!")
        else:
            st.error("PDF export failed.")

if __name__ == "__main__":
    main()
