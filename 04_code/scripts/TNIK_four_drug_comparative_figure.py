import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

input_file = "09_results/tables/TNIK_four_drug_master_structural_comparison.csv"
output_file = "09_results/figures/TNIK_four_drug_comparative_interactions.png"

df = pd.read_csv(input_file)
df = df.set_index("Drug").loc[
    ["Ponatinib", "Tepotinib", "Entrectinib", "Nilotinib"]
].reset_index()

fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle(
    "Comparative Docking and Structural Interaction Analysis of Four TNIK Candidates",
    fontsize=14, fontweight="bold"
)

plots = [
    ("Docking_score_kcal_mol", "Docking score (kcal/mol)",
     "Docking score\n(more negative is more favourable)"),
    ("Geometry_Hbond_candidates", "Hydrogen-bond geometry candidates",
     "Potential hydrogen-bond\ngeometry candidates"),
    ("Nearby_residues_within_4A", "Nearby residues within 4 Å",
     "Nearby TNIK residues\n(within 4 Å)"),
    ("Carbon_atom_pairs_within_4A", "Carbon atom pairs within 4 Å",
     "Potential nonpolar carbon\ncontacts (atom pairs)")
]

for ax, (column, ylabel, title) in zip(axes.flat, plots):
    values = df[column]
    bars = ax.bar(df["Drug"], values)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=15)
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{value:g}",
            ha="center",
            va="bottom" if value >= 0 else "top",
            fontsize=9
        )

fig.text(
    0.5, 0.015,
    "Computational indicators from docking poses; contact counts are not binding-affinity measurements.",
    ha="center", fontsize=9
)

plt.tight_layout(rect=[0, 0.045, 1, 0.93])
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close()

print("Figure saved successfully:")
print(output_file)
