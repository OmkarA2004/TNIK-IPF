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
    return np.loadtxt(path, comments=("@", "#"))

summary_rows = []
contact_data = {}
distance_data = {}

for label, prefix in systems.items():
    contact_file = folder / f"{prefix}_contact_count.xvg"
    distance_file = folder / f"{prefix}_min_distance.xvg"

    if not contact_file.exists() or not distance_file.exists():
        raise FileNotFoundError(
            f"Missing output for {label}: {contact_file} or {distance_file}"
        )

    contacts = read_xvg(contact_file)
    distances = read_xvg(distance_file)

    contact_data[label] = contacts
    distance_data[label] = distances

    summary_rows.append({
        "System": label,
        "Frames": len(contacts),
        "Mean_contacts": np.mean(contacts[:, 1]),
        "SD_contacts": np.std(contacts[:, 1], ddof=1),
        "Minimum_contacts": np.min(contacts[:, 1]),
        "Maximum_contacts": np.max(contacts[:, 1]),
        "Final_contacts": contacts[-1, 1],
        "Mean_min_distance_nm": np.mean(distances[:, 1]),
        "Final_min_distance_nm": distances[-1, 1],
    })

summary = pd.DataFrame(summary_rows)
summary_file = folder / "TNIK_ligand_contact_summary.csv"
summary.to_csv(summary_file, index=False)

print("\nPROTEIN-LIGAND CONTACT SUMMARY")
print(summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print(f"\nSaved: {summary_file}")

# Contact count over time
plt.figure(figsize=(10, 6))
for label, data in contact_data.items():
    plt.plot(data[:, 0] / 1000, data[:, 1], label=label, linewidth=1.7)

plt.xlabel("Simulation time (ns)")
plt.ylabel("Protein-ligand atom-pair contacts")
plt.title("TNIK Complexes: Contacts within 0.35 nm")
plt.xlim(0, 1)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(plot_dir / "TNIK_ligand_contact_counts.png", dpi=300)
plt.close()

# Minimum distance over time
plt.figure(figsize=(10, 6))
for label, data in distance_data.items():
    plt.plot(data[:, 0] / 1000, data[:, 1], label=label, linewidth=1.7)

plt.xlabel("Simulation time (ns)")
plt.ylabel("Minimum protein-ligand distance (nm)")
plt.title("TNIK Complexes: Minimum Protein-Ligand Distance")
plt.xlim(0, 1)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(plot_dir / "TNIK_ligand_minimum_distances.png", dpi=300)
plt.close()

print("\nPlots saved:")
print(plot_dir / "TNIK_ligand_contact_counts.png")
print(plot_dir / "TNIK_ligand_minimum_distances.png")
