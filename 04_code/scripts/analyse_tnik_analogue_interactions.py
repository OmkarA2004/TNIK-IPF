from pathlib import Path
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

ROOT = Path.cwd()

RECEPTOR = ROOT / "06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt"
POSE_DIR = ROOT / "06_docking/virtual_screen/analogue_candidates"
OUT = ROOT / "06_docking/interaction_analysis/analogue_candidates"

OUT.mkdir(parents=True, exist_ok=True)

CANDIDATES = [
    "TNIK_BRICS_0001",
    "TNIK_BRICS_0110",
    "TNIK_BRICS_0112",
    "TNIK_BRICS_0179",
    "TNIK_BRICS_0189",
]

def read_heavy_atoms(path, first_model_only=False):
    atoms = []
    inside_first_model = False
    model_seen = False

    with open(path) as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if model_seen:
                    break
                model_seen = True
                inside_first_model = True
                continue

            if line.startswith("ENDMDL"):
                if first_model_only and inside_first_model:
                    break
                continue

            if first_model_only and model_seen and not inside_first_model:
                continue

            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            atom_type = line[77:79].strip().upper()
            if not atom_type:
                atom_type = line.split()[-1].upper()

            if atom_type in ("H", "HD"):
                continue

            try:
                atoms.append({
                    "atom": line[12:16].strip(),
                    "residue": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "residue_number": line[22:26].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "type": atom_type,
                })
            except (ValueError, IndexError):
                continue

    return pd.DataFrame(atoms)

if not RECEPTOR.exists():
    raise FileNotFoundError(f"Receptor missing: {RECEPTOR}")

receptor = read_heavy_atoms(RECEPTOR)

if receptor.empty:
    raise ValueError("No receptor heavy atoms could be read.")

receptor_xyz = receptor[["x", "y", "z"]].to_numpy()

all_nearby = []
all_polar = []
all_hydrophobic = []
all_clashes = []
summary = []

for candidate in CANDIDATES:
    pose_file = POSE_DIR / f"{candidate}_docked.pdbqt"

    if not pose_file.exists():
        print(f"Skipping missing pose: {candidate}")
        continue

    ligand = read_heavy_atoms(pose_file, first_model_only=True)

    if ligand.empty:
        print(f"Skipping empty ligand pose: {candidate}")
        continue

    ligand_xyz = ligand[["x", "y", "z"]].to_numpy()
    distances = cdist(ligand_xyz, receptor_xyz)

    # Minimum distance from each receptor atom to any ligand atom.
    min_distances = distances.min(axis=0)

    nearby = receptor.copy()
    nearby["minimum_distance_A"] = min_distances
    nearby = nearby[nearby["minimum_distance_A"] <= 4.0].copy()

    if not nearby.empty:
        nearby["candidate_id"] = candidate
        nearby = (
            nearby.groupby(
                ["candidate_id", "chain", "residue", "residue_number"],
                as_index=False
            )["minimum_distance_A"].min()
            .sort_values("minimum_distance_A")
        )
        all_nearby.append(nearby)

    # Potential polar contacts: N/O atom pairs within 3.5 A.
    polar_pairs = []
    for i, lig_atom in ligand.iterrows():
        if lig_atom["type"] not in ("N", "NA", "O", "OA"):
            continue

        for j, rec_atom in receptor.iterrows():
            if rec_atom["type"] not in ("N", "NA", "O", "OA"):
                continue

            d = distances[i, j]
            if d <= 3.5:
                polar_pairs.append({
                    "candidate_id": candidate,
                    "ligand_atom": lig_atom["atom"],
                    "ligand_type": lig_atom["type"],
                    "receptor_residue": rec_atom["residue"],
                    "receptor_chain": rec_atom["chain"],
                    "receptor_residue_number": rec_atom["residue_number"],
                    "receptor_atom": rec_atom["atom"],
                    "receptor_type": rec_atom["type"],
                    "distance_A": round(float(d), 3),
                })

    all_polar.extend(polar_pairs)

    # Potential hydrophobic contacts: carbon-type atom pairs within 4 A.
    hydrophobic_pairs = []
    for i, lig_atom in ligand.iterrows():
        if lig_atom["type"] not in ("C", "A"):
            continue

        for j, rec_atom in receptor.iterrows():
            if rec_atom["type"] not in ("C", "A"):
                continue

            d = distances[i, j]
            if d <= 4.0:
                hydrophobic_pairs.append({
                    "candidate_id": candidate,
                    "ligand_atom": lig_atom["atom"],
                    "receptor_residue": rec_atom["residue"],
                    "receptor_chain": rec_atom["chain"],
                    "receptor_residue_number": rec_atom["residue_number"],
                    "receptor_atom": rec_atom["atom"],
                    "distance_A": round(float(d), 3),
                })

    all_hydrophobic.extend(hydrophobic_pairs)

    # Very close heavy-atom contacts: flag for visual inspection.
    clash_i, clash_j = np.where(distances < 2.0)
    clashes = []

    for i, j in zip(clash_i, clash_j):
        clashes.append({
            "candidate_id": candidate,
            "ligand_atom": ligand.iloc[i]["atom"],
            "receptor_residue": receptor.iloc[j]["residue"],
            "receptor_chain": receptor.iloc[j]["chain"],
            "receptor_residue_number": receptor.iloc[j]["residue_number"],
            "receptor_atom": receptor.iloc[j]["atom"],
            "distance_A": round(float(distances[i, j]), 3),
        })

    all_clashes.extend(clashes)

    summary.append({
        "candidate_id": candidate,
        "ligand_heavy_atoms": len(ligand),
        "receptor_heavy_atoms": len(receptor),
        "nearby_residues_within_4A": (
            nearby[["residue", "residue_number"]].drop_duplicates().shape[0]
        ),
        "potential_polar_atom_pairs_within_3_5A": len(polar_pairs),
        "potential_hydrophobic_atom_pairs_within_4A": len(hydrophobic_pairs),
        "very_close_contacts_below_2A": len(clashes),
        "minimum_interatomic_distance_A": round(float(distances.min()), 3),
    })

    print(
        f"{candidate}: polar={len(polar_pairs)}, "
        f"hydrophobic={len(hydrophobic_pairs)}, "
        f"very-close={len(clashes)}"
    )

def save_table(rows, filename):
    path = OUT / filename

    if isinstance(rows, pd.DataFrame):
        df = rows.copy()
    elif rows and all(isinstance(item, pd.DataFrame) for item in rows):
        df = pd.concat(rows, ignore_index=True)
    else:
        df = pd.DataFrame.from_records(rows)

    df.to_csv(path, index=False)
    print(f"Saved: {path} ({len(df)} rows)")
    return df

summary_df = save_table(
    summary, "TNIK_analogue_interaction_summary.csv"
)
save_table(all_nearby, "TNIK_analogue_nearby_residues.csv")
save_table(all_polar, "TNIK_analogue_potential_polar_contacts.csv")
save_table(all_hydrophobic, "TNIK_analogue_hydrophobic_contacts.csv")
save_table(all_clashes, "TNIK_analogue_very_close_contacts.csv")

print("\n===== INTERACTION SUMMARY =====")
print(summary_df.to_string(index=False))

print("\nSaved interaction reports in:")
print(OUT)
