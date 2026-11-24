import subprocess
import csv
import re
from pathlib import Path

root = Path.home() / "Desktop/TNIK_IPF_Project"
outdir = root / "06_docking/reproducibility"
outdir.mkdir(parents=True, exist_ok=True)

receptor = root / "06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt"

drugs = {
    "Ponatinib": 1593,
    "Tepotinib": 2015,
    "Entrectinib": 1969,
    "Nilotinib": 1310
}

seeds = [11, 22, 33, 44, 55]
results = []

for drug, drug_id in drugs.items():
    ligand = root / f"06_docking/ligand_prep/batch_pdbqt/drug_{drug_id}.pdbqt"

    for seed in seeds:
        name = f"{drug.lower()}_seed_{seed}"
        pose = outdir / f"{name}_out.pdbqt"
        log = outdir / f"{name}.log"

        cmd = [
            "vina", "--receptor", str(receptor),
            "--ligand", str(ligand),
            "--center_x", "10.794",
            "--center_y", "10.872",
            "--center_z", "6.540",
            "--size_x", "22", "--size_y", "22", "--size_z", "22",
            "--exhaustiveness", "8",
            "--seed", str(seed),
            "--out", str(pose)
        ]

        print(f"Running {drug}, seed {seed}...", flush=True)
        run = subprocess.run(cmd, capture_output=True, text=True)
        output = run.stdout + "\n" + run.stderr
        log.write_text(output)

        match = re.search(r"^\s*1\s+(-?\d+\.\d+)", output, re.MULTILINE)

        if run.returncode == 0 and match:
            score = float(match.group(1))
            print(f"Score: {score:.3f} kcal/mol")
            results.append([drug, seed, score])
        else:
            print(f"Check log: {log.name}")
            results.append([drug, seed, "FAILED"])

with open(outdir / "TNIK_reproducibility_scores.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Drug", "Seed", "Docking_score_kcal_mol"])
    writer.writerows(results)

print("Finished. Results saved in:", outdir)
