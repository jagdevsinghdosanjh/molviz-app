import streamlit as st
from rdkit import Chem # noqa
from rdkit.Chem import rdDepictor # noqa
from rdkit.Chem import Draw

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
    "Lactic Acid": "CC(C(=O)O)O",
    "Phenol": "c1ccccc1O",
    "Aniline": "c1ccccc1N"
}

# Page setup
st.set_page_config(page_title="Molecular Geometry Explorer", layout="wide")
st.title("🔬 Molecular Geometry Explorer")

# Molecule selector
col1, col2 = st.columns([2, 1])
with col1:
    selected = st.selectbox("Choose a sample molecule:", list(sample_smiles.keys()))
    smiles = sample_smiles[selected]
    st.text_input("SMILES (editable):", value=smiles, key="smiles_input")

smiles = st.session_state.smiles_input
view_mode = st.radio("View Mode:", ["3D", "2D"])

# Parse and validate molecule
mol = parse_molecule(smiles)
if mol is None:
    st.error("Invalid SMILES string. Please try again.")
else:
    geometry_data = compute_geometry(mol)
    render_ui(geometry_data)

    if view_mode == "3D":
        render_3d(mol)
    else:
        img = Draw.MolToImage(mol, size=(300, 300))
        st.image(img, caption="2D Structure")
        st.image(img, caption="2D Structure")

    if st.button("📤 Export Geometry as PDF"):
        filename = f"{selected.replace(' ', '_')}.pdf"
        success = export_geometry_pdf(filename, geometry_data, smiles)
        if success:
            st.success("PDF exported successfully!")
        else:
            st.error("PDF export failed.")

# SMILES documentation panel
with st.expander("📘 SMILES Documentation"):
    st.markdown("""
### 🔬 What is SMILES?

**SMILES** stands for **Simplified Molecular Input Line Entry System**. It’s a compact way to represent chemical structures using plain text. Instead of drawing molecules, we encode atoms, bonds, rings, and stereochemistry in a linear string format.

---

### 🧪 SMILES Syntax Rules

| Feature            | Rule / Symbol                  | Example                    |
|--------------------|--------------------------------|----------------------------|
| **Atoms**          | Use atomic symbols (`C`, `O`, `N`) | `CO` for methanol         |
| **Bonds**          | Single (implied), double `=`, triple `#` | `C=C` for ethene          |
| **Branches**       | Use parentheses for side chains | `CC(C)C` for isobutane     |
| **Rings**          | Use numbers to indicate closures | `C1CCCCC1` for cyclohexane |
| **Aromatic atoms** | Lowercase letters (`c`, `n`)    | `c1ccccc1` for benzene     |
| **Charges**        | Use brackets with `+` or `-`     | `[NH4+]`, `[O-]`           |
| **Stereochemistry**| Use `@` and `@@` for chiral centers | `C[C@H](O)C(=O)O` for lactic acid |
| **Isotopes**       | Prefix with isotope number       | `[13C]`, `[2H]`            |

---

### 🧠 Why Use SMILES?

- ✅ Compact and readable
- ✅ Easy to store, search, and share
- ✅ Compatible with cheminformatics tools like RDKit
- ✅ Ideal for batch processing and molecular visualization

---

### 📌 Notes for Students

- SMILES does **not** show 3D geometry—just connectivity and stereochemistry
- Different SMILES can represent the same molecule (canonical SMILES solves this)
- Always validate SMILES before using them in modeling or visualization

---

### 📚 Example SMILES Library

Name,SMILES Methane,C Ethanol,CCO Benzene,c1ccccc1 Acetic Acid,CC(=O)O Lactic Acid,CC(C(=O)O)O
""")
