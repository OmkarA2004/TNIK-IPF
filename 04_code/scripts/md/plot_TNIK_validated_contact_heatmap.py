from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

folder = Path(
    "06_docking/molecular_dynamics/analysis/npt_stability/residue_contacts"
)

source = folder / "TNIK_residue_contact_comparison_VALIDATED_PDB.csv"
out_png = folder / "TNIK_validated_residue_contact_heatmap.png"
out_csv = folder / "TNIK_top_contact_residues_for_heatmap.csv"

systems = [
    "Ponatinib_N32_H",
    "Ponatinib_N37_H",
    "Tepotinib_plus1",
]

df = pd.read_csv(source)

required = [
    "pdb_residue_number",
    "pdb_residue_name",
    *systems,
]

missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing expected columns: {missing}")

# Select residues with at least 50% occupancy in one or more systems.
selected = df[df[systems].max(axis=1) >= 50].copy()
selected = selected.sort_values("pdb_residue_number")

selected["residue_label"] = (
    selected["pdb_residue_name"].astype(str)
    + selected["pdb_residue_number"].astype(int).astype(str)
)

selected.to_csv(out_csv, index=False)

matrix = selected.set_index("residue_label")[systems]
matrix.columns = [
    "Ponatinib N32-H",
    "Ponatinib N37-H",
    "Tepotinib +1",
]

fig_height = max(6, len(matrix) * 0.34)
fig, ax = plt.subplots(figsize=(9, fig_height))

image = ax.imshow(
    matrix.values,
    aspect="auto",
    interpolation="nearest",
    vmin=0,
    vmax=100,
)

ax.set_xticks(range(len(matrix.columns)))
ax.set_xticklabels(matrix.columns, rotation=15, ha="right")

ax.set_yticks(range(len(matrix.index)))
ax.set_yticklabels(matrix.index)

ax.set_title(
    "TNIK–Ligand Residue Contact Occupancy\nPDB 6RA7 Chain A"
)
ax.set_xlabel("Ligand protonation model")
ax.set_ylabel("TNIK residue (original PDB numbering)")

# Annotate each cell with its occupancy percentage.
for row in range(matrix.shape[0]):
    for col in range(matrix.shape[1]):
        value = matrix.iloc[row, col]
        ax.text(
            col, row, f"{value:.0f}",
            ha="center", va="center", fontsize=8
        )

colorbar = fig.colorbar(image, ax=ax)
colorbar.set_label("Contact occupancy (%)")

fig.tight_layout()
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)

print("CONTACT HEATMAP COMPLETED")
print("Image:", out_png)
print("Selected residue data:", out_csv)
print("Residues included:", len(matrix))
print("\nTop contact residues:")
print(
    selected[
        ["pdb_residue_number", "pdb_residue_name", *systems]
    ].head(20).to_string(index=False)
)
