# 🔬 Molecular Geometry Explorer

**molviz-app** is an interactive Streamlit-based tool for visualizing molecular geometry from SMILES input. Designed for educators, coordinators, and chemistry learners, it computes bond angles, lengths, and renders 3D structures with export-ready reports.

---

## 🚀 Features

- 🧪 **SMILES Parsing** — Input or select molecules like ethanol, methane, benzene, and more
- 📐 **Geometry Computation** — Calculates bond angles and bond lengths using RDKit
- 🧬 **3D Visualization** — Interactive molecular models rendered via py3Dmol
- 🖼️ **2D Viewer Toggle** — Switch between 3D and annotated 2D diagrams
- 📤 **PDF Export** — Generate printable geometry reports for handouts or documentation
- 📁 **Batch Mode (Coming Soon)** — Upload `.csv` of SMILES for bulk processing

---

## 🧰 Tech Stack

| Component      | Purpose                          |
|----------------|----------------------------------|
| Streamlit      | UI framework                     |
| RDKit          | Molecular parsing & geometry     |
| py3Dmol        | 3D visualization                 |
| reportlab      | PDF generation                   |
| pandas, numpy  | Data handling                    |
| pytest         | Modular testing                  |

---

## 📦 Installation

```bash
git clone https://github.com/jagdevsinghdosanjh/molviz-app
cd molviz-app
python -m venv .venv
.\.venv\Scripts\activate  # On Windows
pip install -r requirements.txt
streamlit run main.py

🧪 Testing
bash
pytest tests/
Includes tests for:

SMILES parsing

Geometry computation

3D rendering

📁 Project Structure
Code
molviz-app/
├── app/
│   ├── parser.py        # SMILES to RDKit molecule
│   ├── geometry.py      # Bond angles and lengths
│   ├── visualizer.py    # 3D viewer with py3Dmol
│   ├── exporter.py      # PDF export logic
│   ├── ui.py            # Streamlit layout components
│   ├── config.py        # Centralized settings
├── main.py              # App entry point
├── tests/               # Pytest suite
├── requirements.txt     # Dependencies
├── README.md            # This file
📚 Educational Use
This tool is ideal for:

Chemistry classrooms

ICT coordinators

Curriculum designers

Self-learners exploring molecular structure

🧠 Credits
Developed by Jagdev Singh Dosanjh, educator and modular app designer. Built to empower clarity, reproducibility, and student engagement in molecular chemistry.

📜 License
MIT License. See LICENSE file for details.
MIT License

Copyright (c) 2025 Jagdev Singh Dosanjh : https://dosanjhpubsasr.org

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


