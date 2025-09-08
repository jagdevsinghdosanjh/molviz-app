import streamlit as st
from app.ui import render_ui
from app.parser import parse_molecule
from app.geometry import compute_geometry
from app.visualizer import render_3d
from app.exporter import export_geometry_pdf

# Export button and logic will be placed inside main() after geometry_data and smiles are defined

def main():
    st.set_page_config(page_title="Molecular Geometry Visualizer", layout="wide")
    st.title("🔬 Molecular Geometry Explorer")

    smiles = st.text_input("Enter SMILES string:", "CCO")  # Example: ethanol
    if smiles:
        mol = parse_molecule(smiles)
        if mol:
            geometry_data = compute_geometry(mol)
            render_ui(geometry_data)
            render_3d(mol)
        else:
            st.error("Invalid SMILES string. Please try again.")
            geometry_data = compute_geometry(mol)
            render_ui(geometry_data)
            render_3d(mol)
            if st.button("📤 Export Geometry as PDF"):
                success = export_geometry_pdf("geometry_report.pdf", geometry_data, smiles)
                if success:
                    st.success("PDF exported successfully!")
                else:
                    st.error("PDF export failed.")
        else:
            st.error("Invalid SMILES string. Please try again.")
