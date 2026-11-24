from pathlib import Path
import pandas as pd

folder = Path(
    "06_docking/molecular_dynamics/analysis/npt_stability/residue_contacts"
)

input_file = folder / "TNIK_all_complexes_residue_contact_occupancy.csv"
comparison_file = folder / "TNIK_residue_contact_comparison.csv"

df = pd.read_csv(input_file)

# Confirmed sequence correspondence:
# GROMACS residue 1 = original 6RA7 residue 13.
# Both structures contain 302 protein residues.
df["original_PDB_residue_number"] = df["residue_number"] + 12

# Put original PDB numbering first in the output.
columns = [
    "system",
    "residue_number",
    "original_PDB_residue_number",
    "residue_name",
    "frames_in_contact",
    "total_frames",
    "occupancy_percent",
]
df = df[columns]

# Save mapped complete results.
mapped_file = folder / "TNIK_all_complexes_residue_contact_occupancy_PDB_numbering.csv"
df.to_csv(mapped_file, index=False)

# Create mapped comparison table.
comparison = df.pivot_table(
    index=["residue_number", "original_PDB_residue_number", "residue_name"],
    columns="system",
    values="occupancy_percent",
    fill_value=0
).reset_index()

comparison_file_out = folder / "TNIK_residue_contact_comparison_PDB_numbering.csv"
comparison.to_csv(comparison_file_out, index=False)

# Create concise report.
report = folder / "TNIK_residue_contact_PDB_numbering_report.txt"

with report.open("w") as f:
    f.write("TNIK RESIDUE CONTACT ANALYSIS — PDB NUMBERING\n")
    f.write("=" * 65 + "\n")
    f.write("Original structure: PDB 6RA7, chain A\n")
    f.write("GROMACS residue 1 corresponds to original PDB residue 13.\n")
    f.write("Mapping: PDB residue number = GROMACS residue number + 12.\n\n")

    for system in df["system"].unique():
        subset = df[df["system"] == system].nlargest(
            20, "occupancy_percent"
        )

        f.write("\n" + system + "\n")
        f.write("-" * 50 + "\n")
        f.write(
            subset[
                [
                    "residue_number",
                    "original_PDB_residue_number",
                    "residue_name",
                    "frames_in_contact",
                    "occupancy_percent",
                ]
            ].to_string(index=False)
        )
        f.write("\n")

print("RESIDUE MAPPING COMPLETED")
print("Mapped full results:", mapped_file)
print("Mapped comparison:", comparison_file_out)
print("Mapped report:", report)
print("\nFirst 15 mapped rows:")
print(df.head(15).to_string(index=False))
