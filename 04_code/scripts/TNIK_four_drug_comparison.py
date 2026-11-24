from pathlib import Path
import pandas as pd
from scipy.spatial.distance import cdist

BASE = Path.cwd()
RECEPTOR = BASE / "06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt"
OUT = BASE / "09_results/tables"
OUT.mkdir(parents=True, exist_ok=True)

DRUGS = {
    "Ponatinib": 1593,
    "Tepotinib": 2015,
    "Entrectinib": 1969,
    "Nilotinib": 1310
}

def read_atoms(path, first_model_only=False):
    atoms = []
    inside = not first_model_only
    model_seen = False

    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if first_model_only:
                    if model_seen:
                        break
                    model_seen = True
                    inside = True
                continue

            if line.startswith("ENDMDL"):
                if first_model_only:
                    break
                continue

            if not inside:
                continue

            if line.startswith(("ATOM", "HETATM")):
                typ = line[77:79].strip().upper()
                if typ in ("H", "HD"):
                    continue

                atoms.append({
                    "atom": line[12:16].strip(),
                    "residue": line[17:20].strip(),
                    "residue_number": line[22:26].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "type": typ
                })

    return pd.DataFrame(atoms)

receptor = read_atoms(RECEPTOR)
rec_xyz = receptor[["x", "y", "z"]].to_numpy()

summary = []

for drug, drug_id in DRUGS.items():
    print(f"\nAnalysing {drug} (ID {drug_id})...")

    ligand_file = (
        BASE / f"06_docking/virtual_screen/top10_redocking/poses/drug_{drug_id}_out.pdbqt"
    )

    if not ligand_file.exists():
        print("Missing pose file:", ligand_file)
        continue

    ligand = read_atoms(ligand_file, first_model_only=True)
    lig_xyz = ligand[["x", "y", "z"]].to_numpy()
    distances = cdist(lig_xyz, rec_xyz)

    # Minimum distance for each nearby residue
    nearby = []
    for j, atom in receptor.iterrows():
        d = distances[:, j].min()
        if d <= 4.0:
            nearby.append({
                "residue": atom["residue"],
                "residue_number": atom["residue_number"],
                "minimum_distance_A": d
            })

    nearby_df = pd.DataFrame(nearby)
    nearby_df = (
        nearby_df.groupby(["residue", "residue_number"], as_index=False)
        ["minimum_distance_A"].min()
        .sort_values("minimum_distance_A")
    )
    nearby_df["minimum_distance_A"] = nearby_df["minimum_distance_A"].round(3)

    # N/O contacts within 3.5 A
    contacts = []
    for i, la in ligand.iterrows():
        if la["type"] not in ("N", "NA", "O", "OA"):
            continue

        for j, ra in receptor.iterrows():
            if ra["type"] not in ("N", "NA", "O", "OA"):
                continue

            d = distances[i, j]
            if d <= 3.5:
                contacts.append({
                    "ligand_atom": la["atom"],
                    "ligand_type": la["type"],
                    "receptor_residue": ra["residue"],
                    "receptor_residue_number": ra["residue_number"],
                    "receptor_atom": ra["atom"],
                    "receptor_type": ra["type"],
                    "distance_A": round(d, 3)
                })

    contacts_df = pd.DataFrame(contacts)

    # All heavy-atom contacts <= 3.2 A
    close = []
    for i, la in ligand.iterrows():
        for j, ra in receptor.iterrows():
            d = distances[i, j]
            if d <= 3.2:
                close.append({
                    "ligand_atom": la["atom"],
                    "ligand_type": la["type"],
                    "receptor_residue": ra["residue"],
                    "receptor_residue_number": ra["residue_number"],
                    "receptor_atom": ra["atom"],
                    "receptor_type": ra["type"],
                    "distance_A": round(d, 3)
                })

    close_df = pd.DataFrame(close)

    prefix = f"TNIK_{drug.lower()}"

    nearby_df.to_csv(OUT / f"{prefix}_nearby_residues.csv", index=False)
    contacts_df.to_csv(OUT / f"{prefix}_polar_contacts.csv", index=False)
    close_df.to_csv(OUT / f"{prefix}_close_atom_contacts.csv", index=False)

    # Docking score from first model
    score = None
    with open(ligand_file) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                score = float(line.split()[3])
                break

    summary.append({
        "Drug": drug,
        "Drug_ID": drug_id,
        "Docking_score_kcal_mol": score,
        "Ligand_heavy_atoms": len(ligand),
        "Nearby_residues_within_4A": len(nearby_df),
        "Potential_N_O_contacts": len(contacts_df),
        "Closest_heavy_atom_contact_A": round(distances.min(), 3)
    })

    print("Nearby residues:", len(nearby_df))
    print("Potential N/O contacts:", len(contacts_df))
    print("Close atom pairs:", len(close_df))

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUT / "TNIK_four_drug_comparison.csv", index=False)

print("\n========== FOUR-DRUG COMPARISON ==========")
print(summary_df.to_string(index=False))
print("\nAll individual and combined files saved in:", OUT)
