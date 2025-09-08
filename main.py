# import streamlit as st
# from app.ui import render_ui
# from app.parser import parse_molecule
# from app.geometry import compute_geometry
# from app.visualizer import render_3d
# from app.exporter import export_geometry_pdf
import streamlit as st
from app.ui import render_ui
from app.parser import parse_molecule
from app.geometry import compute_geometry
from app.visualizer import render_3d
from app.exporter import export_geometry_pdf

def main():
    st.set_page_config(page_title="Molecular Geometry Visualizer", layout="wide")
    st.title("🔬 Molecular Geometry Explorer")

    smiles = st.text_input("Enter SMILES string:", "CCO")  # Example: ethanol

    if not smiles:
        st.warning("Please enter a SMILES string to begin.")
        return

    mol = parse_molecule(smiles)
    if mol is None:
        st.error("Invalid SMILES string. Please try again.")
        return

    geometry_data = compute_geometry(mol)
    render_ui(geometry_data)
    render_3d(mol)

    if st.button("📤 Export Geometry as PDF"):
        success = export_geometry_pdf("geometry_report.pdf", geometry_data, smiles)
        if success:
            st.success("PDF exported successfully!")
        else:
            st.error("PDF export failed.")

if __name__ == "__main__":
    main()
