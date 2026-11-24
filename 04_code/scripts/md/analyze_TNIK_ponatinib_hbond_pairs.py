from pathlib import Path
from collections import Counter
import re
import pandas as pd

base = Path("06_docking/molecular_dynamics")

systems = {
    "Ponatinib_N32_H": {
        "gro": base / "energy_minimization/ponatinib_N32_H/npt_1ns/npt_1ns.gro",
        "ndx": base / "analysis/npt_stability/hbond_pairs/ponatinib_N32_H_framewise_hbonds.ndx",
    },
    "Ponatinib_N37_H": {
        "gro": base / "energy_minimization/ponatinib_N37_H/npt_1ns/npt_1ns.gro",
        "ndx": base / "analysis/npt_stability/hbond_pairs/ponatinib_N37_H_framewise_hbonds.ndx",
    },
}

outdir = base / "analysis/npt_stability/hbond_pairs"
outdir.mkdir(parents=True, exist_ok=True)

all_rows = []

for system, files in systems.items():
    gro = files["gro"]
    ndx = files["ndx"]

    if not gro.exists() or not ndx.exists():
        raise FileNotFoundError(f"Missing input for {system}: {gro} or {ndx}")

    # Map GROMACS atom indices to residue and atom identities.
    atom_map = {}
    lines = gro.read_text().splitlines()

    for line in lines[2:-1]:
        try:
            index = int(line[15:20])
            residue_number = int(line[0:5])
            residue_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_map[index] = (residue_number, residue_name, atom_name)
        except ValueError:
            continue

    # Parse frame-specific hydrogen-bond sections.
    sections = re.split(r"(?m)^\[\s*(.*?)\s*\]\s*$", ndx.read_text())
    frame_data = {}

    for i in range(1, len(sections), 2):
        name = sections[i]
        data = sections[i + 1]
        match = re.search(r"frame_(\d+)", name)

        if not match:
            continue

        frame = int(match.group(1))
        triples = []

        for line in data.splitlines():
            values = line.split()
            if len(values) == 3 and all(v.isdigit() for v in values):
                triples.append(tuple(map(int, values)))

        frame_data[frame] = triples

    total_frames = len(frame_data)
    counts = Counter(
        triple
        for triples in frame_data.values()
        for triple in set(triples)
    )

    print(f"\n===== {system} =====")
    print("Frames:", total_frames)
    print("Distinct atom triples:", len(counts))

    for triple, count in counts.most_common():
        identities = [atom_map.get(index, ("NOT FOUND", "", "")) for index in triple]
        protein_atoms = [x for x in identities if x[1] != "MOL"]
        ligand_atoms = [x for x in identities if x[1] == "MOL"]

        row = {
            "system": system,
            "atom_indices": "-".join(map(str, triple)),
            "frames_present": count,
            "total_frames": total_frames,
            "occupancy_percent": 100 * count / total_frames,
            "atom_1": str(identities[0]),
            "atom_2": str(identities[1]),
            "atom_3": str(identities[2]),
            "protein_atoms": "; ".join(map(str, protein_atoms)),
            "ligand_atoms": "; ".join(map(str, ligand_atoms)),
        }
        all_rows.append(row)

        print(
            f"{triple} | {count}/{total_frames} frames "
            f"({100 * count / total_frames:.2f}%)"
        )
        for index, identity in zip(triple, identities):
            print(f"  Atom {index}: {identity}")

    if total_frames != 101:
        print("WARNING: Expected 101 frames; found", total_frames)

df = pd.DataFrame(all_rows)
outfile = outdir / "ponatinib_hbond_pair_occupancy.csv"
df.to_csv(outfile, index=False)

print("\nSaved:", outfile)
