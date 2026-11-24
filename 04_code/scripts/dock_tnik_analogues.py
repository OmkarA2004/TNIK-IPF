from pathlib import Path
import subprocess
import shutil
import re
import pandas as pd

receptor = Path("06_docking/receptor_prep/TNIK_6RA7_rigid.pdbqt")
ligand_dir = Path("06_docking/ligand_prep/analogue_candidates")
outdir = Path("06_docking/virtual_screen/analogue_candidates")
outdir.mkdir(parents=True, exist_ok=True)

vina = shutil.which("vina")

if vina is None:
    raise SystemExit("ERROR: AutoDock Vina not found.")

if not receptor.exists():
    raise SystemExit(f"ERROR: Receptor not found: {receptor}")

candidates = [
    "TNIK_BRICS_0001",
    "TNIK_BRICS_0110",
    "TNIK_BRICS_0112",
    "TNIK_BRICS_0179",
    "TNIK_BRICS_0189"
]

results = []

for candidate in candidates:
    ligand = ligand_dir / f"{candidate}.pdbqt"
    output = outdir / f"{candidate}_docked.pdbqt"
    log_file = outdir / f"{candidate}_vina.log"

    if not ligand.exists():
        print(f"Skipping {candidate}: ligand file missing.")
        continue

    print(f"\n===== Docking {candidate} =====")

    command = [
        vina,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", "10.794",
        "--center_y", "10.872",
        "--center_z", "6.540",
        "--size_x", "22",
        "--size_y", "22",
        "--size_z", "22",
        "--exhaustiveness", "32",
        "--seed", "42",
        "--num_modes", "20",
        "--energy_range", "5",
        "--out", str(output)
    ]

    run = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    log_file.write_text(run.stdout + "\n" + run.stderr)
    print(run.stdout)

    if run.returncode != 0 or not output.exists():
        print(f"FAILED: {candidate}")
        print(run.stderr[-1500:])
        results.append({
            "candidate_id": candidate,
            "status": "Failed",
            "best_affinity_kcal_mol": None
        })
        continue

    scores = re.findall(
        r"^\s*\d+\s+(-?\d+\.\d+)",
        run.stdout,
        flags=re.MULTILINE
    )

    best_score = float(scores[0]) if scores else None

    results.append({
        "candidate_id": candidate,
        "status": "Docked",
        "best_affinity_kcal_mol": best_score,
        "docked_pose_file": str(output),
        "log_file": str(log_file)
    })

summary = pd.DataFrame(results)
summary_file = outdir / "TNIK_analogue_docking_summary.csv"
summary.to_csv(summary_file, index=False)

print("\n===== ANALOGUE DOCKING SUMMARY =====")
print(summary.to_string(index=False))
print("\nSaved:", summary_file)
