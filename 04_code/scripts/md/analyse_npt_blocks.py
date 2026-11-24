import numpy as np
from pathlib import Path

path = Path(
    "06_docking/molecular_dynamics/analysis/"
    "TNIK_npt_600ps_properties.xvg"
)

rows = []
for line in path.read_text().splitlines():
    line = line.strip()
    if line and not line.startswith(("#", "@")):
        rows.append([float(x) for x in line.split()])

a = np.asarray(rows)

names = [
    "Potential (kJ/mol)",
    "Temperature (K)",
    "Pressure (bar)",
    "Volume (nm3)",
    "Density (kg/m3)"
]

print("NPT FINAL BLOCK ANALYSIS")
print(f"Frames: {len(a)}")
print(f"Time: {a[0,0]:.1f}–{a[-1,0]:.1f} ps")

for start, end in [(300,400), (400,500), (500,600)]:
    block = a[(a[:,0] >= start) & (a[:,0] <= end)]
    print(f"\n{start}–{end} ps")
    for col, name in enumerate(names, start=1):
        print(
            f"{name}: mean={np.mean(block[:,col]):.3f}, "
            f"SD={np.std(block[:,col], ddof=1):.3f}"
        )

report = Path(
    "06_docking/molecular_dynamics/logs/"
    "TNIK_npt_final_block_analysis.txt"
)
report.write_text("Final 300 ps block analysis completed.\n")
print("\nAnalysis completed.")
