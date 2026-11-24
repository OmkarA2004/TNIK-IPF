from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT = Path.cwd()
DATA_DIR = (
    PROJECT
    / "06_docking/molecular_dynamics/analysis/npt_stability"
)
PLOT_DIR = DATA_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

systems = {
    "Ponatinib N32-H": "ponatinib_N32_H_backbone_RMSD.xvg",
    "Ponatinib N37-H": "ponatinib_N37_H_backbone_RMSD.xvg",
    "Tepotinib +1": "tepotinib_pH74_plus1_backbone_RMSD.xvg",
}

trajectories = {}
summary = []

for name, filename in systems.items():
    filepath = DATA_DIR / filename

    if not filepath.exists():
        raise FileNotFoundError(f"Missing RMSD file: {filepath}")

    data = np.loadtxt(
        filepath,
        comments=("@", "#"),
        usecols=(0, 1)
    )

    time_ns = data[:, 0]
    rmsd_nm = data[:, 1]

    trajectories[name] = (time_ns, rmsd_nm)

    first_200 = rmsd_nm[time_ns <= 0.2]
    last_200 = rmsd_nm[time_ns >= 0.8]

    summary.append({
        "System": name,
        "Frames": len(rmsd_nm),
        "Mean_RMSD_nm": np.mean(rmsd_nm),
        "SD_RMSD_nm": np.std(rmsd_nm, ddof=1),
        "Maximum_RMSD_nm": np.max(rmsd_nm),
        "Final_RMSD_nm": rmsd_nm[-1],
        "First_200ps_mean_nm": np.mean(first_200),
        "Last_200ps_mean_nm": np.mean(last_200),
    })

summary_df = pd.DataFrame(summary)

csv_path = DATA_DIR / "TNIK_backbone_RMSD_summary.csv"
summary_df.to_csv(csv_path, index=False)

plt.figure(figsize=(10, 6))

for name, (time_ns, rmsd_nm) in trajectories.items():
    plt.plot(time_ns, rmsd_nm, linewidth=1.8, label=name)

plt.xlabel("Simulation time (ns)")
plt.ylabel("Backbone RMSD (nm)")
plt.title("TNIK–Ligand Complexes: Backbone RMSD during 1 ns NPT")
plt.xlim(0, 1.0)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

plot_path = PLOT_DIR / "TNIK_three_complexes_backbone_RMSD.png"
plt.savefig(plot_path, dpi=300)
plt.close()

print("\nBACKBONE RMSD SUMMARY")
print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

print(f"\nCSV saved: {csv_path}")
print(f"Plot saved: {plot_path}")
