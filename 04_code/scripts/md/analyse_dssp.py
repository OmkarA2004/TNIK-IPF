import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

base = Path("06_docking/molecular_dynamics/analysis")
xvg = base / "TNIK_apo_secondary_structure.xvg"
out = Path("09_results/molecular_dynamics")
out.mkdir(parents=True, exist_ok=True)

# Read category names and numerical data
labels = []
data = []

with open(xvg) as f:
    for line in f:
        line = line.strip()

        match = re.match(r'@\s+s(\d+)\s+legend\s+"(.*)"', line)
        if match:
            labels.append((int(match.group(1)), match.group(2)))

        if line and not line.startswith(("#", "@")):
            try:
                values = [float(x) for x in line.split()]
                data.append(values)
            except ValueError:
                pass

if not data:
    raise ValueError("No numerical data found in DSSP XVG file.")

arr = np.array(data)
labels = [name for _, name in sorted(labels)]

if arr.shape[1] != len(labels) + 1:
    raise ValueError(
        f"Column mismatch: {arr.shape[1]-1} data columns, "
        f"{len(labels)} labels. Please check the XVG file."
    )

df = pd.DataFrame(arr, columns=["Time_ns"] + labels)

# Calculate category percentages at every frame
category_cols = labels
totals = df[category_cols].sum(axis=1)

if (totals <= 0).any():
    raise ValueError("Some frames have zero total secondary-structure counts.")

percent = df[category_cols].div(totals, axis=0) * 100
percent.insert(0, "Time_ns", df["Time_ns"])

# Save numerical results
percent.to_csv(out / "TNIK_apo_DSSP_percentages.csv", index=False)

summary = percent[category_cols].agg(["mean", "std", "min", "max"]).T
summary.to_csv(out / "TNIK_apo_DSSP_summary.csv")

# Plot percentage composition over time
plt.figure(figsize=(12, 6))

for col in category_cols:
    plt.plot(percent["Time_ns"], percent[col], label=col)

plt.xlabel("Time (ns)")
plt.ylabel("Secondary-structure composition (%)")
plt.title("Apo-TNIK Secondary-Structure Composition")
plt.xlim(percent["Time_ns"].min(), percent["Time_ns"].max())
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
plt.tight_layout()

plt.savefig(
    out / "TNIK_apo_DSSP_composition.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Save a compact report
with open(out / "TNIK_apo_DSSP_report.txt", "w") as f:
    f.write("APO-TNIK DSSP ANALYSIS\n")
    f.write("=====================\n")
    f.write(f"Frames analysed: {len(df)}\n")
    f.write(f"Trajectory duration: {df['Time_ns'].iloc[-1]:.3f} ns\n\n")
    f.write("Mean secondary-structure composition (%):\n")
    f.write(summary["mean"].round(2).to_string())
    f.write("\n\nStandard deviation (%):\n")
    f.write(summary["std"].round(2).to_string())

print("DSSP analysis completed successfully.")
print(f"Frames analysed: {len(df)}")
print(f"Duration: {df['Time_ns'].iloc[-1]:.3f} ns")
print("\nMean composition (%):")
print(summary["mean"].round(2).to_string())
print("\nFiles saved in:", out)
