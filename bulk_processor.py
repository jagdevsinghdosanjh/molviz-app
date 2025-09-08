import pandas as pd
from app.parser import parse_molecule
from app.geometry import compute_geometry
from app.exporter import export_geometry_pdf
import zipfile
import os

def bundle_pdfs_to_zip(pdf_folder: str, zip_name: str = "molecule_reports.zip") -> str:
    with zipfile.ZipFile(zip_name, "w") as zipf:
        for filename in os.listdir(pdf_folder):
            if filename.endswith(".pdf"):
                zipf.write(os.path.join(pdf_folder, filename), arcname=filename)
    return zip_name


def process_molecule_row(name: str, smiles: str, export: bool = True) -> dict:
    mol = parse_molecule(smiles)
    if mol is None:
        return {"name": name, "smiles": smiles, "error": "Invalid SMILES"}

    geometry = compute_geometry(mol)
    if export:
        filename = f"{name.replace(' ', '_')}.pdf"
        export_geometry_pdf(filename, geometry, smiles)

    return {
        "name": name,
        "smiles": smiles,
        "geometry": geometry,
        "status": "Success"
    }

def process_csv(path: str, export: bool = True) -> list:
    df = pd.read_csv(path)
    results = []
    for _, row in df.iterrows():
        name = row.get("Name", "Unnamed")
        smiles = row.get("SMILES", "")
        result = process_molecule_row(name, smiles, export)
        results.append(result)
    return results
