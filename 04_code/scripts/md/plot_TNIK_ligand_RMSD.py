from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

folder = Path("06_docking/molecular_dynamics/analysis/npt_stability")
plot_dir = folder / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

systems = {
    "Ponatinib N32-H": "ponatinib_N32_H_ligand_RMSD.xvg",
    "Ponatinib N37-H": "ponatinib_N37_H_ligand_RMSD.xvg",
    "Tepotinib +1": "tepotinib_pH74_plus1_ligand_RMSD.xvg",
}

plt.figure(figsize=(10, 6))

for label, filename in systems.items():
    data = np.loadtxt(folder / filename, comments=("@", "#"))
    plt.plot(data[:, 0], data[:, 1], linewidth=1.8, label=label)

plt.xlabel("Simulation time (ns)")
plt.ylabel("Ligand heavy-atom RMSD (nm)")
plt.title("TNIK–Ligand Complexes: Ligand RMSD during 1 ns NPT")
plt.xlim(0, 1.0)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

output = plot_dir / "TNIK_three_complexes_ligand_RMSD.png"
plt.savefig(output, dpi=300)
plt.close()

print(f"Plot saved: {output}")
