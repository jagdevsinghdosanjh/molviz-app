from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import HybridizationType
import numpy as np

def compute_geometry(mol):
    """
    Computes bond lengths, bond angles, and dominant geometry type.
    """
    try:
        AllChem.EmbedMolecule(mol, maxAttempts=1000)
        AllChem.UFFOptimizeMolecule(mol)
        conformer = mol.GetConformer()

        bond_lengths = []
        bond_angles = []
        hybridizations = []

        # Bond lengths
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            pos1 = conformer.GetAtomPosition(idx1)
            pos2 = conformer.GetAtomPosition(idx2)
            distance = np.linalg.norm(np.array(pos1) - np.array(pos2))
            bond_lengths.append(distance)

        # Bond angles
        for atom in mol.GetAtoms():
            neighbors = atom.GetNeighbors()
            if len(neighbors) >= 2:
                for i in range(len(neighbors)):
                    for j in range(i + 1, len(neighbors)):
                        idx1 = neighbors[i].GetIdx()
                        idx2 = atom.GetIdx()
                        idx3 = neighbors[j].GetIdx()

                        vec1 = np.array(conformer.GetAtomPosition(idx1)) - np.array(conformer.GetAtomPosition(idx2))
                        vec2 = np.array(conformer.GetAtomPosition(idx3)) - np.array(conformer.GetAtomPosition(idx2))

                        cosine_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
                        angle_rad = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
                        angle_deg = np.degrees(angle_rad)
                        bond_angles.append(angle_deg)

            # Hybridization-based geometry classification
            hyb = atom.GetHybridization()
            if hyb == HybridizationType.SP3:
                hybridizations.append("tetrahedral")
            elif hyb == HybridizationType.SP2:
                hybridizations.append("trigonal planar")
            elif hyb == HybridizationType.SP:
                hybridizations.append("linear")

        # Determine dominant geometry
        if hybridizations:
            geometry_type = max(set(hybridizations), key=hybridizations.count)
        else:
            geometry_type = "unknown"

        return {
            "geometry": geometry_type,
            "angles": [round(angle, 2) for angle in bond_angles],
            "lengths": [round(length, 2) for length in bond_lengths]
        }

    except Exception as error:
        return {
            "geometry": "error",
            "angles": [],
            "lengths": [],
            "error": str(error)
        }
