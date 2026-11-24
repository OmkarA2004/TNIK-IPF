import numpy as np
import pandas as pd
from pathlib import Path

analysis = Path("06_docking/molecular_dynamics/analysis")
results = Path("09_results/molecular_dynamics")
results.mkdir(parents=True, exist_ok=True)

file = analysis / "TNIK_apo_backbone_RMSD.xvg"

data = []
with open(file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        values = line.split()
        data.append([float(values[0]), float(values[1])])

df = pd.DataFrame(data, columns=["Time_ns", "RMSD_nm"])

# Analyse three consecutive 200 ps blocks
blocks = [(0, 0.2), (0.2, 0.4), (0.4, 0.6)]
rows = []

for start, end in blocks:
    subset = df[(df.Time_ns >= start) & (df.Time_ns < end)]
    if subset.empty:
        continue

    rows.append({
        "Block_ns": f"{start:.1f}-{end:.1f}",
        "Frames": len(subset),
        "Mean_RMSD_nm": subset.RMSD_nm.mean(),
        "SD_RMSD_nm": subset.RMSD_nm.std(),
        "Min_RMSD_nm": subset.RMSD_nm.min(),
        "Max_RMSD_nm": subset.RMSD_nm.max()
    })

summary = pd.DataFrame(rows)
summary.to_csv(results / "TNIK_apo_RMSD_block_summary.csv", index=False)

with open(results / "TNIK_apo_RMSD_block_report.txt", "w") as f:
    f.write("APO-TNIK BACKBONE RMSD BLOCK ANALYSIS\n")
    f.write("====================================\n")
    f.write(f"Total frames: {len(df)}\n")
    f.write(f"Trajectory duration: {df.Time_ns.max():.3f} ns\n\n")
    f.write(summary.round(4).to_string(index=False))
    f.write("\n")

print("RMSD block analysis completed.")
print(summary.round(4).to_string(index=False))
print("\nReport saved in:", results)
