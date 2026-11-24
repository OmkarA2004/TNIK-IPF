from pathlib import Path
import numpy as np
import pandas as pd

RECEPTOR = Path("06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt")
LIGAND = Path("06_docking/virtual_screen/top10_redocking/poses/drug_1593_out.pdbqt")

def read_atoms(path, first_model=False):
    atoms = []
    seen_model = False

    for line in open(path):
        if line.startswith("MODEL"):
            if first_model and seen_model:
                break
            seen_model = True
            continue
        if line.startswith("ENDMDL"):
            if first_model:
                break
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue

        atoms.append({
            "name": line[12:16].strip(),
            "residue": line[17:20].strip(),
            "resid": line[22:26].strip(),
            "chain": line[21:22].strip(),
            "type": line[77:79].strip().upper(),
            "xyz": np.array([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])
        })
    return atoms

def angle(a, b, c):
    v1 = a - b
    v2 = c - b
    cosine = np.dot(v1, v2) / (
        np.linalg.norm(v1) * np.linalg.norm(v2)
    )
    return np.degrees(np.arccos(np.clip(cosine, -1, 1)))

def find_hbonds(donor_atoms, all_donor_atoms, acceptors, direction):
    results = []

    for donor in donor_atoms:
        hydrogens = [
            h for h in all_donor_atoms
            if h["type"] == "HD"
            and np.linalg.norm(h["xyz"] - donor["xyz"]) < 1.25
        ]

        for h in hydrogens:
            for acceptor in acceptors:
                da = np.linalg.norm(donor["xyz"] - acceptor["xyz"])
                if da > 3.5:
                    continue

                ha = np.linalg.norm(h["xyz"] - acceptor["xyz"])
                dha = angle(donor["xyz"], h["xyz"], acceptor["xyz"])

                if ha <= 2.7 and dha >= 120:
                    results.append({
                        "Direction": direction,
                        "Donor": donor["name"],
                        "Donor_residue": donor["residue"],
                        "Donor_resid": donor["resid"],
                        "Acceptor": acceptor["name"],
                        "Acceptor_residue": acceptor["residue"],
                        "Acceptor_resid": acceptor["resid"],
                        "H_A_distance_A": round(ha, 3),
                        "D_A_distance_A": round(da, 3),
                        "D_H_A_angle_deg": round(dha, 1)
                    })

    return results

protein = read_atoms(RECEPTOR)
ligand = read_atoms(LIGAND, first_model=True)

# AutoDock atom types: N is a donor; OA and NA are acceptors.
protein_donors = [a for a in protein if a["type"] == "N"]
ligand_donors = [a for a in ligand if a["type"] == "N"]

protein_acceptors = [a for a in protein if a["type"] in ("OA", "NA")]
ligand_acceptors = [a for a in ligand if a["type"] in ("OA", "NA")]

results = []
results += find_hbonds(
    ligand_donors, ligand, protein_acceptors, "Ligand to protein"
)
results += find_hbonds(
    protein_donors, protein, ligand_acceptors, "Protein to ligand"
)

df = pd.DataFrame(results)
out = Path("09_results/tables/TNIK_ponatinib_hbond_geometry.csv")
df.to_csv(out, index=False)

print("Ligand atoms:", len(ligand))
print("Ligand polar hydrogens:", sum(a["type"] == "HD" for a in ligand))
print("Receptor polar hydrogens:", sum(a["type"] == "HD" for a in protein))
print("Hydrogen-bond geometry candidates:", len(df))
print(df.to_string(index=False) if len(df) else "No candidates passed the geometry filters.")
print("Saved:", out)
