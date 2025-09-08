import streamlit.components.v1 as components
import py3Dmol
from rdkit import Chem

def render_3d(mol):
    """
    Renders a 3D molecular viewer using py3Dmol inside Streamlit.
    """
    try:
        mb = Chem.MolToMolBlock(mol)

        viewer = py3Dmol.view(width=600, height=400)
        viewer.addModel(mb, "mol")
        viewer.setStyle({'stick': {}})
        viewer.zoomTo()

        # Embed in Streamlit
        components.html(viewer._make_html(), height=400)

    except Exception as e:
        import streamlit as st
        st.error(f"3D rendering failed: {e}")
