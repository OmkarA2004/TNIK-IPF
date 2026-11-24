import csv
import itertools
import math
from pathlib import Path
from collections import defaultdict

root = Path.home() / "Desktop/TNIK_IPF_Project"
folder = root / "06_docking/reproducibility"
drugs = ["Ponatinib", "Tepotinib", "Entrectinib", "Nilotinib"]
seeds = [11, 22, 33, 44, 55]

def read_first_pose(path):
    atoms = {}
    counts = defaultdict(int)
    in_first_model = True
    saw_model = False

    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            if saw_model:
                break
            saw_model = True
            in_first_model = True
            continue
        if line.startswith("ENDMDL") and saw_model:
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue

        atom_type = line[77:79].strip().upper()
        if not atom_type:
            atom_type = line.split()[-1].upper()
        if atom_type.startswith("H"):
            continue

        atom_name = line[12:16].strip()
        counts[atom_name] += 1
        key = (atom_name, counts[atom_name])
        xyz = tuple(float(line[i:i+8]) for i in (30, 38, 46))
        atoms[key] = xyz

    if not atoms:
        raise ValueError(f"No heavy atoms found in {path}")
    return atoms

def rmsd(a, b):
    if set(a) != set(b):
        raise ValueError("Atom names differ between poses")
    keys = sorted(a)
    squared = [
        sum((a[k][i] - b[k][i])**2 for i in range(3))
        for k in keys
    ]
    return math.sqrt(sum(squared) / len(squared))

pair_rows = []
summary_rows = []

for drug in drugs:
    poses = {}
    for seed in seeds:
        path = folder / f"{drug.lower()}_seed_{seed}_out.pdbqt"
        if not path.exists():
            raise FileNotFoundError(path)
        poses[seed] = read_first_pose(path)

    values = []
    for s1, s2 in itertools.combinations(seeds, 2):
        value = rmsd(poses[s1], poses[s2])
        values.append(value)
        pair_rows.append([drug, s1, s2, round(value, 3)])

    summary_rows.append([
        drug, len(values), round(sum(values)/len(values), 3),
        round(min(values), 3), round(max(values), 3)
    ])

with open(folder / "TNIK_pose_pairwise_RMSD.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Drug", "Seed_1", "Seed_2", "Heavy_atom_RMSD_A"])
    w.writerows(pair_rows)

with open(folder / "TNIK_pose_consistency_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Drug", "Pairwise_comparisons", "Mean_RMSD_A",
                "Min_RMSD_A", "Max_RMSD_A"])
    w.writerows(summary_rows)

print("Pose comparison completed.")
print("Drug          Mean RMSD (Å)   Min (Å)   Max (Å)")
for row in summary_rows:
    print(f"{row[0]:13s} {row[2]:10.3f} {row[3]:10.3f} {row[4]:10.3f}")
print("Results saved in:", folder)
