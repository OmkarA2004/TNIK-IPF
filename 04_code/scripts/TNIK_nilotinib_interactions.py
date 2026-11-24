from pathlib import Path
import pandas as pd
from scipy.spatial.distance import cdist

project = Path.cwd()

ligand_file = project / "06_docking/virtual_screen/top10_redocking/poses/drug_1310_out.pdbqt"
receptor_file = project / "06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt"

def read_atoms(path, first_model_only=False):
    atoms = []
    inside_first_model = not first_model_only
    model_seen = False

    with open(path) as f:
        for line in f:
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

            if line.startswith(("ATOM", "HETATM")):
                atom_type = line[77:79].strip().upper()
                if atom_type in ("H", "HD"):
                    continue

                atoms.append({
                    "atom": line[12:16].strip(),
                    "residue": line[17:20].strip(),
                    "residue_number": line[22:26].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "type": atom_type
                })

    return pd.DataFrame(atoms)

ligand = read_atoms(ligand_file, first_model_only=True)
receptor = read_atoms(receptor_file)

if ligand.empty or receptor.empty:
    raise ValueError("Could not read ligand or receptor atoms.")

ligand_xyz = ligand[["x", "y", "z"]].to_numpy()
receptor_xyz = receptor[["x", "y", "z"]].to_numpy()

distances = cdist(ligand_xyz, receptor_xyz)

# Nearby receptor residues: minimum heavy-atom distance <= 4 Å
nearby = []
for j, atom in receptor.iterrows():
    minimum = distances[:, j].min()
    if minimum <= 4.0:
        nearby.append({
            "residue": atom["residue"],
            "residue_number": atom["residue_number"],
            "minimum_distance_A": round(minimum, 3)
        })

nearby_df = pd.DataFrame(nearby)
nearby_df = (
    nearby_df.groupby(["residue", "residue_number"], as_index=False)
    ["minimum_distance_A"].min()
    .sort_values("minimum_distance_A")
)

# Potential polar contacts: N/O atom pairs within 3.5 Å
contacts = []
for i, lig_atom in ligand.iterrows():
    if lig_atom["type"] not in ("N", "NA", "O", "OA"):
        continue

    for j, rec_atom in receptor.iterrows():
        if rec_atom["type"] not in ("N", "NA", "O", "OA"):
            continue

        distance = distances[i, j]
        if distance <= 3.5:
            contacts.append({
                "ligand_atom": lig_atom["atom"],
                "ligand_type": lig_atom["type"],
                "receptor_residue": rec_atom["residue"],
                "receptor_residue_number": rec_atom["residue_number"],
                "receptor_atom": rec_atom["atom"],
                "receptor_type": rec_atom["type"],
                "distance_A": round(distance, 3)
            })

contacts_df = pd.DataFrame(contacts)

out_dir = project / "09_results/tables"
out_dir.mkdir(parents=True, exist_ok=True)

nearby_df.to_csv(out_dir / "TNIK_nilotinib_nearby_residues.csv", index=False)
contacts_df.to_csv(out_dir / "TNIK_nilotinib_polar_contacts_model1.csv", index=False)

print("Ligand heavy atoms in Model 1:", len(ligand))
print("Receptor heavy atoms:", len(receptor))
print("Nearby residues within 4 Å:", len(nearby_df))
print("Potential N/O atom contacts within 3.5 Å:", len(contacts_df))
print("\nPotential polar contacts:")
print(contacts_df.to_string(index=False))
print("\nSaved corrected results in:", out_dir)
