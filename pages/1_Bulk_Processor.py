import streamlit as st
import pandas as pd
from bulk_processor import process_molecule_row

st.title("📁 Bulk Molecular Processor")

uploaded_file = st.file_uploader("Upload a CSV of molecules", type="csv")
if uploaded_file:
    df = pd.read_csv(uploaded_file)
    results = []
    for _, row in df.iterrows():
        name = row.get("Name", "Unnamed")
        smiles = row.get("SMILES", "")
        result = process_molecule_row(name, smiles, export=True)
        results.append(result)

    for res in results:
        st.subheader(res["name"])
        st.write(f"SMILES: `{res['smiles']}`")
        if "geometry" in res:
            st.json(res["geometry"])
        else:
            st.error(res["error"])
