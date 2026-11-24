
from pathlib import Path
import pandas as pd
from scipy.spatial.distance import cdist

PROJECT = Path.cwd()

RECEPTOR_FILE = (
    PROJECT / "06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt"
)

POSE_DIR = (
    PROJECT / "06_docking/virtual_screen/top10_redocking/poses"
)

OUTPUT_DIR = PROJECT / "09_results/tables"

DRUGS = {
    "Ponatinib": 1593,
    "Tepotinib": 2015,
    "Entrectinib": 1969,
    "Nilotinib": 1310,
}

CUTOFF = 4.0


def read_atoms(path, first_model_only=False):
    atoms = []
    inside_first_model = not first_model_only
    model_seen = False

    with open(path) as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if first_model_only:
                    if model_seen:
                        break
                    model_seen = True
                    inside_first_model = True
                continue

            if line.startswith("ENDMDL"):
                if first_model_only:
                    break
                continue

            if not inside_first_model:
                continue

            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_type = line[77:79].strip().upper()

            if atom_type in ("H", "HD"):
                continue

            atoms.append({
                "atom": line[12:16].strip(),
                "residue": line[17:20].strip(),
                "residue_number": line[22:26].strip(),
                "type": atom_type,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })

    return pd.DataFrame(atoms)


receptor = read_atoms(RECEPTOR_FILE)

if receptor.empty:
    raise ValueError("No receptor atoms were read.")

# Carbon atom types in AutoDock Vina:
# C = aliphatic carbon; A = aromatic carbon.
receptor_nonpolar = receptor[
    receptor["type"].isin(["C", "A"])
].copy()

if receptor_nonpolar.empty:
    raise ValueError("No receptor carbon atoms were found.")

receptor_xyz = receptor_nonpolar[
    ["x", "y", "z"]
].to_numpy()

all_contacts = []
summary = []

for drug_name, drug_id in DRUGS.items():

    ligand_file = POSE_DIR / f"drug_{drug_id}_out.pdbqt"

    if not ligand_file.exists():
        print("Missing ligand pose:", ligand_file)
        continue

    ligand = read_atoms(ligand_file, first_model_only=True)

    ligand_nonpolar = ligand[
        ligand["type"].isin(["C", "A"])
    ].copy()

    if ligand_nonpolar.empty:
        print("No ligand carbon atoms:", drug_name)
        continue

    ligand_xyz = ligand_nonpolar[
        ["x", "y", "z"]
    ].to_numpy()

    distances = cdist(ligand_xyz, receptor_xyz)

    contacts = []

    for i, ligand_atom in enumerate(
        ligand_nonpolar.to_dict("records")
    ):
        for j, receptor_atom in enumerate(
            receptor_nonpolar.to_dict("records")
        ):
            distance = distances[i, j]

            if distance <= CUTOFF:
                contacts.append({
                    "Drug": drug_name,
                    "Drug_ID": drug_id,
                    "Ligand_atom": ligand_atom["atom"],
                    "Ligand_type": ligand_atom["type"],
                    "Receptor_residue": receptor_atom["residue"],
                    "Receptor_residue_number":
                        receptor_atom["residue_number"],
                    "Receptor_atom": receptor_atom["atom"],
                    "Receptor_type": receptor_atom["type"],
                    "Distance_A": round(distance, 3),
                })

    contact_df = pd.DataFrame(contacts)

    if not contact_df.empty:
        all_contacts.extend(contacts)

        unique_residues = contact_df[
            ["Receptor_residue", "Receptor_residue_number"]
        ].drop_duplicates()

        closest_distance = contact_df["Distance_A"].min()
        residue_count = len(unique_residues)
        pair_count = len(contact_df)
    else:
        closest_distance = None
        residue_count = 0
        pair_count = 0

    summary.append({
        "Drug": drug_name,
        "Drug_ID": drug_id,
        "Carbon_atom_pairs_within_4A": pair_count,
        "Unique_receptor_residues": residue_count,
        "Closest_carbon_contact_A": closest_distance,
    })

    print(
        f"{drug_name}: {pair_count} carbon-carbon pairs; "
        f"{residue_count} receptor residues"
    )


OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

contacts_output = (
    OUTPUT_DIR / "TNIK_four_drug_hydrophobic_contacts.csv"
)

summary_output = (
    OUTPUT_DIR / "TNIK_four_drug_hydrophobic_summary.csv"
)

pd.DataFrame(all_contacts).to_csv(
    contacts_output, index=False
)

summary_df = pd.DataFrame(summary)
summary_df.to_csv(summary_output, index=False)

print("\n===== HYDROPHOBIC CONTACT SUMMARY =====")
print(summary_df.to_string(index=False))

print("\nSaved:", contacts_output)
print("Saved:", summary_output)
