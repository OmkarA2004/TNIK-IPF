from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

folder = Path(
    "06_docking/molecular_dynamics/analysis/npt_stability/residue_contacts"
)

source = folder / "TNIK_contact_pattern_classification.csv"

if not source.exists():
    raise FileNotFoundError(source)

df = pd.read_csv(source)

systems = [
    "Ponatinib_N32_H",
    "Ponatinib_N37_H",
    "Tepotinib_plus1",
]

required = [
    "pdb_residue_number",
    "pdb_residue_name",
    "contact_class",
    *systems,
]

missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# Verify that the input contains the expected classification categories.
expected = {
    "Shared by all three",
    "Tepotinib-associated",
    "Ponatinib-associated",
    "Below threshold",
}

observed = set(df["contact_class"].dropna().unique())
if observed != expected:
    raise ValueError(
        f"Unexpected classification categories: {observed}"
    )

# Sort residues by category and original PDB numbering.
category_order = [
    "Shared by all three",
    "Ponatinib-associated",
    "Tepotinib-associated",
    "Below threshold",
]

df["category_order"] = df["contact_class"].map(
    {name: i for i, name in enumerate(category_order)}
)

df = df.sort_values(
    ["category_order", "pdb_residue_number"]
).reset_index(drop=True)

df["residue_label"] = (
    df["pdb_residue_name"].astype(str)
    + df["pdb_residue_number"].astype(int).astype(str)
)

# Save the exact ordered data used in the heatmap.
ordered_csv = folder / "TNIK_contact_classification_heatmap_data.csv"
df.to_csv(ordered_csv, index=False)

matrix = df[systems].astype(float).values

fig_height = max(9, len(df) * 0.38)
fig, ax = plt.subplots(figsize=(11, fig_height))

im = ax.imshow(
    matrix,
    aspect="auto",
    interpolation="nearest",
    cmap="viridis",
    norm=Normalize(vmin=0, vmax=100),
)

ax.set_xticks(range(3))
ax.set_xticklabels(
    ["Ponatinib N32-H", "Ponatinib N37-H", "Tepotinib +1"],
    rotation=15,
    ha="right",
)

ax.set_yticks(range(len(df)))
ax.set_yticklabels(df["residue_label"], fontsize=9)

# Annotate every cell with contact occupancy.
for i in range(len(df)):
    for j in range(3):
        value = matrix[i, j]
        ax.text(
            j, i, f"{value:.1f}",
            ha="center",
            va="center",
            fontsize=7,
            color="white" if value < 45 else "black",
        )

# Draw separators between classification categories.
for i in range(1, len(df)):
    if df.loc[i, "contact_class"] != df.loc[i - 1, "contact_class"]:
        ax.axhline(i - 0.5, color="red", linewidth=1.2)

# Add category names in the right margin.
for category in category_order:
    positions = df.index[df["contact_class"] == category].tolist()
    if positions:
        midpoint = (positions[0] + positions[-1]) / 2
        ax.text(
            3.05,
            midpoint,
            category,
            va="center",
            ha="left",
            fontsize=9,
            clip_on=False,
        )

ax.set_xlim(-0.5, 3.0)
ax.set_title(
    "TNIK–Ligand Contact Classification\n"
    "PDB 6RA7 Chain A | 0.35 nm cutoff | 101 frames",
    pad=18,
)
ax.set_xlabel("Ligand protonation model")
ax.set_ylabel("TNIK residue (original PDB numbering)")

cbar = fig.colorbar(im, ax=ax, pad=0.25)
cbar.set_label("Contact occupancy (%)")

fig.tight_layout()

output = folder / "TNIK_contact_classification_heatmap.png"
fig.savefig(output, dpi=300, bbox_inches="tight")
plt.close(fig)

print("CLASSIFICATION HEATMAP COMPLETED")
print("Heatmap:", output)
print("Ordered data:", ordered_csv)
print("\nResidues included:", len(df))
print("\nClassification counts:")
print(df["contact_class"].value_counts().to_string())
print("\nCategory-wise residues:")
for category in category_order:
    subset = df[df["contact_class"] == category]
    labels = ", ".join(subset["residue_label"].tolist())
    print(f"{category} ({len(subset)}): {labels}")
