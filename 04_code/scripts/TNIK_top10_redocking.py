import subprocess
import pandas as pd
from pathlib import Path

receptor = Path("06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt")
ligand_dir = Path("06_docking/ligand_prep/batch_pdbqt")

top_file = Path("09_results/tables/TNIK_top10_named_candidates.csv")

output_dir = Path("06_docking/virtual_screen/top10_redocking")
pose_dir = output_dir / "poses"
log_dir = output_dir / "logs"

pose_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

top = pd.read_csv(top_file, dtype={"ID": str})
results = []
failures = []

for _, row in top.iterrows():
    drug_id = str(row["ID"]).strip()
    drug_name = row["name_sdf"]

    ligand = ligand_dir / f"drug_{drug_id}.pdbqt"
    pose = pose_dir / f"drug_{drug_id}_out.pdbqt"
    log = log_dir / f"drug_{drug_id}.log"

    print(f"Docking {drug_name} (ID {drug_id})...", flush=True)

    command = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", "10.794",
        "--center_y", "10.872",
        "--center_z", "6.540",
        "--size_x", "22",
        "--size_y", "22",
        "--size_z", "22",
        "--exhaustiveness", "32",
        "--num_modes", "9",
        "--seed", "42",
        "--out", str(pose)
    ]

    try:
        run = subprocess.run(
            command, capture_output=True, text=True, timeout=1800
        )

        log_text = run.stdout + "\n" + run.stderr
        log.write_text(log_text)

        affinity = None
        for line in log_text.splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "1":
                try:
                    affinity = float(parts[1])
                    break
                except ValueError:
                    pass

        if run.returncode == 0 and pose.exists() and affinity is not None:
            results.append({
                "ID": drug_id,
                "Drug_Name": drug_name,
                "Best_Affinity_kcal_mol": affinity
            })
            print(f"  Score: {affinity} kcal/mol", flush=True)
        else:
            failures.append({
                "ID": drug_id,
                "Drug_Name": drug_name,
                "Reason": "Vina error or affinity not found"
            })
            print("  Docking failed. Check log.", flush=True)

    except Exception as e:
        failures.append({
            "ID": drug_id,
            "Drug_Name": drug_name,
            "Reason": str(e)
        })
        print(f"  Error: {e}", flush=True)

pd.DataFrame(results).to_csv(
    output_dir / "TNIK_top10_redocking_results.csv", index=False
)

pd.DataFrame(failures).to_csv(
    output_dir / "TNIK_top10_redocking_failures.csv", index=False
)

print("\nRedocking finished.")
print("Successful:", len(results))
print("Failed:", len(failures))
