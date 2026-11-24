from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path.cwd()
BASE = ROOT / "06_docking/interaction_analysis"
DATA = BASE / "comparative_analysis"
OUT = BASE / "comparative_analysis/figures"
OUT.mkdir(parents=True, exist_ok=True)

details = pd.read_csv(DATA / "TNIK_PLIP_interaction_details.csv")
summary = pd.read_csv(DATA / "TNIK_PLIP_interaction_summary.csv")
comparison = pd.read_csv(DATA / "TNIK_PLIP_MD_comparison.csv")

ligands = ["Ponatinib", "Tepotinib"]

# Figure 1: Interaction type counts
pivot = summary.pivot(
    index="interaction_type",
    columns="ligand",
    values="interaction_records"
).fillna(0)

pivot = pivot.reindex(
    ["Hydrophobic", "Hydrogen bond", "Salt bridge", "Halogen bond",
     "Pi stacking", "Pi-cation", "Water bridge", "Metal complex"],
    fill_value=0
)

ax = pivot.plot(kind="bar", figsize=(10, 6), width=0.75)
ax.set_title("PLIP interaction records in representative TNIK docking poses")
ax.set_xlabel("Interaction type")
ax.set_ylabel("Number of interaction records")
ax.legend(title="Ligand")
ax.tick_params(axis="x", rotation=25)
ax.grid(axis="y", alpha=0.25)
plt.tight_layout()
plt.savefig(OUT / "TNIK_PLIP_interaction_counts.png", dpi=300)
plt.close()

# Figure 2: Residue-level interaction presence
residues = (
    details[["residue_name", "pdb_residue_number"]]
    .drop_duplicates()
    .sort_values("pdb_residue_number")
)
residue_labels = [
    f"{r.residue_name}{r.pdb_residue_number}"
    for r in residues.itertuples()
]

matrix = pd.DataFrame(
    0,
    index=residue_labels,
    columns=ligands,
    dtype=int
)

for row in details.itertuples():
    label = f"{row.residue_name}{row.pdb_residue_number}"
    matrix.loc[label, row.ligand] = 1

fig, ax = plt.subplots(figsize=(7, max(6, len(matrix) * 0.28)))
im = ax.imshow(matrix.values, aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(len(ligands)))
ax.set_xticklabels(ligands)
ax.set_yticks(range(len(matrix)))
ax.set_yticklabels(matrix.index)
ax.set_title("Residues with PLIP-detected interactions")
ax.set_xlabel("Ligand")
ax.set_ylabel("TNIK residue")

for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        ax.text(j, i, "✓" if matrix.iloc[i, j] else "–",
                ha="center", va="center")

fig.colorbar(im, ax=ax, ticks=[0, 1],
             label="Interaction detected (1) / not detected (0)")
plt.tight_layout()
plt.savefig(OUT / "TNIK_PLIP_residue_presence.png", dpi=300)
plt.close()

# Figure 3: MD contact occupancy at PLIP-detected residues
occupancy_columns = {
    "Ponatinib N32-H": "Ponatinib_N32_H_contact_occupancy_percent",
    "Ponatinib N37-H": "Ponatinib_N37_H_contact_occupancy_percent",
    "Tepotinib +1": "Tepotinib_plus1_contact_occupancy_percent",
}

unique_residues = comparison[
    ["residue_name", "pdb_residue_number"]
].drop_duplicates().sort_values("pdb_residue_number")

labels = [
    f"{r.residue_name}{r.pdb_residue_number}"
    for r in unique_residues.itertuples()
]

occ = pd.DataFrame(index=labels)

for display, source in occupancy_columns.items():
    values = {}
    for row in comparison.itertuples():
        label = f"{row.residue_name}{row.pdb_residue_number}"
        values[label] = getattr(row, source)
    occ[display] = pd.Series(values)

fig, ax = plt.subplots(figsize=(9, max(6, len(occ) * 0.32)))
image = ax.imshow(
    occ.values.astype(float),
    aspect="auto",
    vmin=0,
    vmax=100
)

ax.set_xticks(range(len(occ.columns)))
ax.set_xticklabels(occ.columns, rotation=20, ha="right")
ax.set_yticks(range(len(occ.index)))
ax.set_yticklabels(occ.index)
ax.set_title("MD contact occupancy at PLIP-identified TNIK residues")
ax.set_xlabel("MD simulation system")
ax.set_ylabel("TNIK residue")

for i in range(occ.shape[0]):
    for j in range(occ.shape[1]):
        value = occ.iloc[i, j]
        if pd.notna(value):
            ax.text(j, i, f"{value:.0f}",
                    ha="center", va="center", fontsize=7)

fig.colorbar(image, ax=ax, label="MD contact occupancy (%)")
plt.tight_layout()
plt.savefig(OUT / "TNIK_PLIP_MD_contact_occupancy.png", dpi=300)
plt.close()

print("Figures created:")
for path in sorted(OUT.glob("*.png")):
    print(f"{path.name}: {path.stat().st_size:,} bytes")

print("\nFigure directory:", OUT)
print("All figures use validated PDB residue numbering.")
print("MD contact occupancy is not equivalent to PLIP interaction occupancy.")
