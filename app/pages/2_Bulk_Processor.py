import streamlit as st
import pandas as pd # noqa
from bulk_processor import process_csv

st.header("📁 Bulk Molecular Processor")

uploaded_file = st.file_uploader("Upload a CSV of molecules", type="csv")
if uploaded_file:
    import tempfile
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
        tmp.write(uploaded_file.read())
        tmp_path = tmp.name
    results = process_csv(tmp_path, export=True)
    for res in results:
        st.subheader(res["name"])
        st.write(f"SMILES: `{res['smiles']}`")
        if "geometry" in res:
            st.json(res["geometry"])
        else:
            st.error(res["error"])
