import numpy as np
import pandas as pd
from pathlib import Path

source = Path("06_docking/molecular_dynamics/analysis/TNIK_apo_hbond_number.xvg")
output = Path("09_results/molecular_dynamics")
output.mkdir(parents=True, exist_ok=True)

data = []

with open(source) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        values = line.split()
        data.append([float(values[0]), float(values[1])])

df = pd.DataFrame(data, columns=["Time_ns", "Hydrogen_bonds"])

if df.empty:
    raise ValueError("No hydrogen-bond data found.")

blocks = [(0, 0.2), (0.2, 0.4), (0.4, 0.6)]
rows = []

for start, end in blocks:
    subset = df[(df.Time_ns >= start) & (df.Time_ns < end)]
    if not subset.empty:
        rows.append({
            "Time_block_ns": f"{start:.1f}-{end:.1f}",
            "Frames": len(subset),
            "Mean_Hbonds": subset.Hydrogen_bonds.mean(),
            "SD": subset.Hydrogen_bonds.std(),
            "Minimum": subset.Hydrogen_bonds.min(),
            "Maximum": subset.Hydrogen_bonds.max()
        })

summary = pd.DataFrame(rows)

df.to_csv(output / "TNIK_apo_hbond_timeseries.csv", index=False)
summary.to_csv(output / "TNIK_apo_hbond_block_summary.csv", index=False)

with open(output / "TNIK_apo_hbond_report.txt", "w") as f:
    f.write("APO-TNIK INTRAPROTEIN HYDROGEN-BOND ANALYSIS\n")
    f.write("============================================\n")
    f.write(f"Frames: {len(df)}\n")
    f.write(f"Duration: {df.Time_ns.max():.3f} ns\n")
    f.write(f"Overall mean: {df.Hydrogen_bonds.mean():.2f}\n")
    f.write(f"Overall SD: {df.Hydrogen_bonds.std():.2f}\n")
    f.write(f"Minimum: {df.Hydrogen_bonds.min():.0f}\n")
    f.write(f"Maximum: {df.Hydrogen_bonds.max():.0f}\n\n")
    f.write(summary.round(2).to_string(index=False))

print("Hydrogen-bond analysis completed.")
print(f"Frames: {len(df)}")
print(f"Overall mean: {df.Hydrogen_bonds.mean():.2f}")
print(f"Overall SD: {df.Hydrogen_bonds.std():.2f}")
print("\nBlock-wise results:")
print(summary.round(2).to_string(index=False))
print("\nReports saved in:", output)
