
from pathlib import Path
import numpy as np
import pandas as pd

RECEPTOR = Path("06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt")

DRUGS = {
    "Ponatinib": 1593,
    "Tepotinib": 2015,
    "Entrectinib": 1969,
    "Nilotinib": 1310,
}

POSE_DIR = Path("06_docking/virtual_screen/top10_redocking/poses")
OUTPUT_DIR = Path("09_results/tables")


def read_atoms(path, first_model=False):
    atoms = []
    seen_model = False

    with open(path) as handle:
        for line in handle:
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

    denominator = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denominator == 0:
        return np.nan

    cosine = np.dot(v1, v2) / denominator
    return np.degrees(np.arccos(np.clip(cosine, -1, 1)))


def find_hbonds(donor_atoms, all_donor_atoms, acceptors, direction):
    results = []

    for donor in donor_atoms:
        hydrogens = [
            h for h in all_donor_atoms
            if h["type"] == "HD"
            and np.linalg.norm(h["xyz"] - donor["xyz"]) < 1.25
        ]

        for hydrogen in hydrogens:
            for acceptor in acceptors:
                da = np.linalg.norm(donor["xyz"] - acceptor["xyz"])

                if da > 3.5:
                    continue

                ha = np.linalg.norm(
                    hydrogen["xyz"] - acceptor["xyz"]
                )

                dha = angle(
                    donor["xyz"],
                    hydrogen["xyz"],
                    acceptor["xyz"]
                )

                if ha <= 2.7 and dha >= 120:
                    results.append({
                        "Direction": direction,
                        "Donor_atom": donor["name"],
                        "Donor_residue": donor["residue"],
                        "Donor_resid": donor["resid"],
                        "Hydrogen_atom": hydrogen["name"],
                        "Acceptor_atom": acceptor["name"],
                        "Acceptor_residue": acceptor["residue"],
                        "Acceptor_resid": acceptor["resid"],
                        "H_A_distance_A": round(ha, 3),
                        "D_A_distance_A": round(da, 3),
                        "D_H_A_angle_deg": round(dha, 1)
                    })

    return results


protein = read_atoms(RECEPTOR)

protein_donors = [
    a for a in protein if a["type"] == "N"
]

protein_acceptors = [
    a for a in protein if a["type"] in ("OA", "NA")
]

all_results = []
summary = []

for drug_name, drug_id in DRUGS.items():
    ligand_path = POSE_DIR / f"drug_{drug_id}_out.pdbqt"

    if not ligand_path.exists():
        print(f"Missing pose: {ligand_path}")
        continue

    ligand = read_atoms(ligand_path, first_model=True)

    ligand_donors = [
        a for a in ligand if a["type"] == "N"
    ]

    ligand_acceptors = [
        a for a in ligand if a["type"] in ("OA", "NA")
    ]

    results = []

    results += find_hbonds(
        ligand_donors,
        ligand,
        protein_acceptors,
        "Ligand to protein"
    )

    results += find_hbonds(
        protein_donors,
        protein,
        ligand_acceptors,
        "Protein to ligand"
    )

    for result in results:
        result["Drug"] = drug_name
        result["Drug_ID"] = drug_id

    all_results.extend(results)

    summary.append({
        "Drug": drug_name,
        "Drug_ID": drug_id,
        "Ligand_atoms": len(ligand),
        "Ligand_polar_hydrogens": sum(
            a["type"] == "HD" for a in ligand
        ),
        "Geometry_candidates": len(results)
    })

    print(
        f"{drug_name}: {len(results)} "
        "geometrically plausible candidates"
    )


OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

results_df = pd.DataFrame(all_results)
summary_df = pd.DataFrame(summary)

results_file = OUTPUT_DIR / "TNIK_four_drug_hbond_geometry.csv"
summary_file = OUTPUT_DIR / "TNIK_four_drug_hbond_geometry_summary.csv"

results_df.to_csv(results_file, index=False)
summary_df.to_csv(summary_file, index=False)

print("\n===== SUMMARY =====")
print(summary_df.to_string(index=False))

print("\n===== GEOMETRY CANDIDATES =====")
if len(results_df):
    print(results_df.to_string(index=False))
else:
    print("No candidates passed the geometry filters.")

print("\nSaved:", results_file)
print("Saved:", summary_file)
