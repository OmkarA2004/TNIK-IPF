import numpy as np
from pathlib import Path

path = Path(
    "06_docking/molecular_dynamics/analysis/"
    "TNIK_npt_properties.xvg"
)

data = []
for line in path.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith(("#", "@")):
        continue
    values = [float(x) for x in line.split()]
    data.append(values)

a = np.asarray(data)
if a.ndim != 2 or a.shape[1] < 6:
    raise SystemExit(f"Unexpected XVG dimensions: {a.shape}")

# Columns: time, Potential, Temperature, Pressure, Volume, Density
names = ["Potential", "Temperature", "Pressure", "Volume", "Density"]
report = []

report.append(f"Frames: {len(a)}")
report.append(f"Time range: {a[0,0]:.2f}–{a[-1,0]:.2f} ps")

for i, name in enumerate(names, start=1):
    x = a[:, i]
    first = x[a[:, 0] < 50]
    last = x[a[:, 0] >= 50]

    report.append(f"\n{name}")
    report.append(f"  Overall mean: {np.mean(x):.4f}")
    report.append(f"  Overall SD:   {np.std(x, ddof=1):.4f}")
    report.append(f"  First 50 ps:  {np.mean(first):.4f}")
    report.append(f"  Last 50 ps:   {np.mean(last):.4f}")
    report.append(f"  Last 20 ps:   {np.mean(x[a[:,0] >= 80]):.4f}")

out = Path(
    "06_docking/molecular_dynamics/logs/"
    "TNIK_npt_stability_report.txt"
)
out.write_text("\n".join(report) + "\n")
print("NPT stability report created.")
print("\n".join(report))
