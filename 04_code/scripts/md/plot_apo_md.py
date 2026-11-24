import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

base = Path("06_docking/molecular_dynamics/analysis")
out = Path("09_results/molecular_dynamics")
out.mkdir(parents=True, exist_ok=True)

def read_xvg(path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith(("#", "@")):
            rows.append([float(x) for x in line.split()])
    return np.array(rows)

plots = [
    ("TNIK_apo_backbone_RMSD.xvg", "Backbone RMSD (nm)",
     "TNIK Apo: Backbone RMSD", "TNIK_apo_RMSD.png"),
    ("TNIK_apo_Rg.xvg", "Radius of gyration (nm)",
     "TNIK Apo: Radius of Gyration", "TNIK_apo_Rg.png"),
]

for filename, ylabel, title, output in plots:
    data = read_xvg(base / filename)
    time_ns = data[:, 0] / 1000
    plt.figure(figsize=(8, 4.5))
    plt.plot(time_ns, data[:, 1], linewidth=1.2)
    plt.xlabel("Time (ns)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out / output, dpi=300)
    plt.close()

data = read_xvg(base / "TNIK_apo_RMSF.xvg")
plt.figure(figsize=(9, 4.5))
plt.plot(data[:, 0], data[:, 1], linewidth=1.2)
plt.xlabel("C-alpha residue index")
plt.ylabel("RMSF (nm)")
plt.title("TNIK Apo: Residue-wise RMSF")
plt.tight_layout()
plt.savefig(out / "TNIK_apo_RMSF.png", dpi=300)
plt.close()

print("Created three figures in:", out)
