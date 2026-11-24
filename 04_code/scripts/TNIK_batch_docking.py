import subprocess
import pandas as pd
from pathlib import Path

receptor = "06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt"
ligand_dir = Path("06_docking/ligand_prep/batch_pdbqt")
output_dir = Path("06_docking/virtual_screen/poses")
log_dir = Path("06_docking/virtual_screen/logs")

output_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

results = []
failures = []

ligands = sorted(ligand_dir.glob("drug_*.pdbqt"))

for i, ligand in enumerate(ligands, 1):
    ligand_id = ligand.stem.replace("drug_", "")
    output = output_dir / f"{ligand.stem}_out.pdbqt"
    log_file = log_dir / f"{ligand.stem}.log"

    command = [
        "vina",
        "--receptor", receptor,
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
        "--out", str(output)
    ]

    try:
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=1800
        )

        log_text = result.stdout + "\n" + result.stderr
        log_file.write_text(log_text)

        if result.returncode == 0 and output.exists():
            affinity = None

            for line in log_text.splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "1":
                    try:
                        affinity = float(parts[1])
                        break
                    except ValueError:
                        pass

            if affinity is not None:
                results.append((ligand_id, affinity))
            else:
                failures.append((ligand_id, "Affinity not found in log"))
        else:
            failures.append((ligand_id, "Vina error or missing output"))

    except subprocess.TimeoutExpired:
        failures.append((ligand_id, "Timeout after 30 minutes"))
    except Exception as e:
        failures.append((ligand_id, str(e)))

    if i % 25 == 0:
        pd.DataFrame(
            results, columns=["ID", "Best_Affinity_kcal_mol"]
        ).to_csv(
            "06_docking/virtual_screen/TNIK_docking_results_partial.csv",
            index=False
        )
        print(f"Processed {i}/{len(ligands)} | "
              f"Success: {len(results)} | Failed: {len(failures)}",
              flush=True)

pd.DataFrame(
    results, columns=["ID", "Best_Affinity_kcal_mol"]
).to_csv(
    "06_docking/virtual_screen/TNIK_docking_results.csv",
    index=False
)

pd.DataFrame(
    failures, columns=["ID", "Failure_Reason"]
).to_csv(
    "06_docking/virtual_screen/TNIK_docking_failures.csv",
    index=False
)

print("\nBatch docking finished.")
print("Successful:", len(results))
print("Failed:", len(failures))
