
import os
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# TNIK-IPF PROJECT: SECTION 10 MD FIGURE GENERATION
# All figures saved at 300 DPI
# ============================================================

OUTPUT_DIR = "09_results/figures/MD_Figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

systems = [
    "Ponatinib N32-H",
    "Ponatinib N37-H",
    "Tepotinib +1"
]

# Consistent figure settings
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 14,
    "figure.dpi": 150,
    "savefig.dpi": 300
})

def save_figure(filename):
    path = os.path.join(OUTPUT_DIR, filename)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", path)


# ============================================================
# FIGURE 1: PROTEIN BACKBONE RMSD
# ============================================================

backbone_mean = [0.021963781, 0.022270339, 0.021968940]
backbone_sd = [0.002235927, 0.002264980, 0.002233077]

plt.figure(figsize=(8, 5))
plt.bar(
    systems,
    backbone_mean,
    yerr=backbone_sd,
    capsize=6
)
plt.ylabel("Backbone RMSD (nm)")
plt.xlabel("TNIK–ligand system")
plt.title("TNIK Protein Backbone RMSD Comparison")
plt.grid(axis="y", linestyle="--", alpha=0.4)
save_figure("MD_Figure_1_TNIK_Backbone_RMSD.png")


# ============================================================
# FIGURE 2: LIGAND RMSD
# ============================================================

ligand_mean = [0.031710042, 0.029538217, 0.029146428]
ligand_sd = [0.005891407, 0.004869704, 0.005645948]

plt.figure(figsize=(8, 5))
plt.bar(
    systems,
    ligand_mean,
    yerr=ligand_sd,
    capsize=6
)
plt.ylabel("Ligand RMSD (nm)")
plt.xlabel("TNIK–ligand system")
plt.title("Ligand RMSD Comparison")
plt.grid(axis="y", linestyle="--", alpha=0.4)
save_figure("MD_Figure_2_TNIK_Ligand_RMSD.png")


# ============================================================
# FIGURE 3: LIGAND–PROTEIN CONTACTS
# ============================================================

contact_mean = [31.198020, 32.009901, 33.376238]
contact_sd = [1.265068, 1.473059, 0.925759]

plt.figure(figsize=(8, 5))
plt.bar(
    systems,
    contact_mean,
    yerr=contact_sd,
    capsize=6
)
plt.ylabel("Mean number of contacts")
plt.xlabel("TNIK–ligand system")
plt.title("TNIK–Ligand Contact Comparison")
plt.grid(axis="y", linestyle="--", alpha=0.4)
save_figure("MD_Figure_3_TNIK_Ligand_Protein_Contacts.png")


# ============================================================
# FIGURE 4: NPT TEMPERATURE
# ============================================================

temperature_mean = [299.979539, 299.983732, 299.940157]
temperature_sd = [1.317537, 1.303057, 1.307184]

plt.figure(figsize=(8, 5))
plt.bar(
    systems,
    temperature_mean,
    yerr=temperature_sd,
    capsize=6
)
plt.axhline(
    300,
    linestyle="--",
    label="Target temperature: 300 K"
)
plt.ylabel("Temperature (K)")
plt.xlabel("TNIK–ligand system")
plt.title("Mean NPT Temperature Comparison")
plt.legend()
plt.grid(axis="y", linestyle="--", alpha=0.4)
save_figure("MD_Figure_4_TNIK_NPT_Temperature.png")


# ============================================================
# FIGURE 5: NPT DENSITY
# ============================================================

density_mean = [1001.568405, 1001.178543, 1001.754674]
density_sd = [26.093785, 27.431012, 25.752516]

plt.figure(figsize=(8, 5))
plt.bar(
    systems,
    density_mean,
    yerr=density_sd,
    capsize=6
)
plt.ylabel("Density (kg/m³)")
plt.xlabel("TNIK–ligand system")
plt.title("Mean NPT Density Comparison")
plt.grid(axis="y", linestyle="--", alpha=0.4)
save_figure("MD_Figure_5_TNIK_NPT_Density.png")


# ============================================================
# FIGURE 6: TEPOTINIB HYDROGEN-BOND PAIR OCCUPANCY
# ============================================================

hb_pairs = [
    "GLU17–Tepotinib",
    "CYS96–Tepotinib"
]

occupancy = [79.207921, 83.168317]

plt.figure(figsize=(8, 5))
bars = plt.bar(hb_pairs, occupancy)

plt.ylabel("Occupancy (%)")
plt.xlabel("Hydrogen-bond interaction pair")
plt.title("Tepotinib Hydrogen-Bond Pair Occupancy")
plt.ylim(0, 100)
plt.grid(axis="y", linestyle="--", alpha=0.4)

for bar, value in zip(bars, occupancy):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 1,
        f"{value:.2f}%",
        ha="center"
    )

save_figure("MD_Figure_6_TNIK_Tepotinib_HBond_Occupancy.png")


# ============================================================
# COMPLETION
# ============================================================

print("\n======================================")
print("ALL SIX MD FIGURES GENERATED")
print("======================================")
print("Output folder:", OUTPUT_DIR)

