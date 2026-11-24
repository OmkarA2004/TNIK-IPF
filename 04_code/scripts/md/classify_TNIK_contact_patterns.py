from pathlib import Path
import pandas as pd

folder = Path(
    "06_docking/molecular_dynamics/analysis/npt_stability/residue_contacts"
)

source = folder / "TNIK_residue_contact_comparison_VALIDATED_PDB.csv"
df = pd.read_csv(source)

n32 = "Ponatinib_N32_H"
n37 = "Ponatinib_N37_H"
tep = "Tepotinib_plus1"

# Contact threshold: 50% of trajectory frames.
threshold = 50.0

def classify(row):
    a, b, c = row[n32], row[n37], row[tep]

    if a >= threshold and b >= threshold and c >= threshold:
        return "Shared by all three"
    if a >= threshold and b >= threshold and c < threshold:
        return "Ponatinib-associated"
    if c >= threshold and a < threshold and b < threshold:
        return "Tepotinib-associated"
    if a >= threshold and b < threshold:
        return "N32-H favoured"
    if b >= threshold and a < threshold:
        return "N37-H favoured"
    if c >= threshold:
        return "Tepotinib-involved"
    if a >= threshold or b >= threshold:
        return "Ponatinib-involved"
    return "Below threshold"

df["contact_class"] = df.apply(classify, axis=1)
df["ponatinib_model_difference"] = (df[n32] - df[n37]).abs()

df["residue"] = (
    df["pdb_residue_name"].astype(str)
    + df["pdb_residue_number"].astype(int).astype(str)
)

out = folder / "TNIK_contact_pattern_classification.csv"
df.sort_values(
    ["contact_class", "pdb_residue_number"]
).to_csv(out, index=False)

summary = folder / "TNIK_contact_pattern_summary.txt"
with summary.open("w") as f:
    f.write("TNIK CONTACT PATTERN CLASSIFICATION\n")
    f.write("Threshold: 50% trajectory-frame occupancy\n")
    f.write("Residue numbering: original PDB 6RA7 chain A\n\n")

    for category, group in df.groupby("contact_class"):
        f.write(f"\n{category} ({len(group)} residues)\n")
        f.write(", ".join(group["residue"].tolist()) + "\n")

    f.write("\nLargest Ponatinib protonation-model differences:\n")
    top = df.nlargest(10, "ponatinib_model_difference")
    f.write(
        top[
            ["residue", n32, n37, "ponatinib_model_difference"]
        ].to_string(index=False)
    )

print("CONTACT CLASSIFICATION COMPLETED")
print("CSV:", out)
print("Summary:", summary)
print("\nClassification counts:")
print(df["contact_class"].value_counts().to_string())
print("\nLargest Ponatinib model differences:")
print(
    df.nlargest(10, "ponatinib_model_difference")[
        ["residue", n32, n37, "ponatinib_model_difference"]
    ].to_string(index=False)
)
