from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

folder = Path("06_docking/molecular_dynamics/analysis/npt_stability")
plot_dir = folder / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

systems = {
    "Ponatinib N32-H": "ponatinib_N32_H",
    "Ponatinib N37-H": "ponatinib_N37_H",
    "Tepotinib +1": "tepotinib_pH74_plus1",
}

def read_xvg(path):
    data = np.loadtxt(path, comments=("@", "#"))
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected XVG format: {path}")
    return data[:, :2]

def summarise(label, direction, data):
    time = data[:, 0]
    hbonds = data[:, 1]

    return {
        "System": label,
        "Direction": direction,
        "Frames": len(hbonds),
        "Mean_H_bonds": np.mean(hbonds),
        "SD_H_bonds": np.std(hbonds, ddof=1),
        "Maximum_H_bonds": np.max(hbonds),
        "Frames_with_H_bonds_percent": 100 * np.mean(hbonds > 0),
        "First_200ps_mean": np.mean(hbonds[time <= 200]),
        "Last_200ps_mean": np.mean(hbonds[time >= 800]),
    }

all_data = {}
summary_rows = []

for label, prefix in systems.items():
    protein_to_ligand = read_xvg(
        folder / f"{prefix}_protein_to_ligand_hbonds.xvg"
    )
    ligand_to_protein = read_xvg(
        folder / f"{prefix}_ligand_to_protein_hbonds.xvg"
    )

    if len(protein_to_ligand) != len(ligand_to_protein):
        raise ValueError(f"Frame-count mismatch for {label}")

    if not np.allclose(protein_to_ligand[:, 0],
                       ligand_to_protein[:, 0]):
        raise ValueError(f"Time-point mismatch for {label}")

    all_data[label] = {
        "Protein to ligand": protein_to_ligand,
        "Ligand to protein": ligand_to_protein,
    }

    summary_rows.append(
        summarise(label, "Protein to ligand", protein_to_ligand)
    )
    summary_rows.append(
        summarise(label, "Ligand to protein", ligand_to_protein)
    )

    combined = protein_to_ligand.copy()
    combined[:, 1] += ligand_to_protein[:, 1]

    summary_rows.append(
        summarise(label, "Combined directions", combined)
    )

summary = pd.DataFrame(summary_rows)
output = folder / "TNIK_hydrogen_bond_summary.csv"
summary.to_csv(output, index=False)

print("\nHYDROGEN-BOND SUMMARY")
print(summary.to_string(index=False,
                        float_format=lambda x: f"{x:.3f}"))
print(f"\nSaved: {output}")

# Plot both directions separately
fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True)

for ax, (label, directions) in zip(axes, all_data.items()):
    for direction, data in directions.items():
        ax.plot(data[:, 0] / 1000, data[:, 1],
                linewidth=1.5, label=direction)

    ax.set_title(label)
    ax.set_ylabel("H-bond count")
    ax.grid(alpha=0.3)
    ax.legend()

axes[-1].set_xlabel("Simulation time (ns)")
axes[-1].set_xlim(0, 1)

plt.suptitle("TNIK–Ligand Hydrogen Bonds During 1 ns NPT",
             fontsize=14)
plt.tight_layout()

plot1 = plot_dir / "TNIK_hydrogen_bonds_by_direction.png"
plt.savefig(plot1, dpi=300)
plt.close()

# Plot combined directional counts
plt.figure(figsize=(10, 6))

for label, directions in all_data.items():
    combined = directions["Protein to ligand"].copy()
    combined[:, 1] += directions["Ligand to protein"][:, 1]

    plt.plot(combined[:, 0] / 1000, combined[:, 1],
             linewidth=1.7, label=label)

plt.xlabel("Simulation time (ns)")
plt.ylabel("Combined directional H-bond count")
plt.title("TNIK Complexes: Combined Hydrogen-Bond Counts")
plt.xlim(0, 1)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()

plot2 = plot_dir / "TNIK_combined_hydrogen_bonds.png"
plt.savefig(plot2, dpi=300)
plt.close()

print("\nPlots saved:")
print(plot1)
print(plot2)
