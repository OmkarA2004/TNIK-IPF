
import os
import pandas as pd

# Project directory
project_dir = os.path.expanduser("~/Desktop/TNIK_IPF_Project")

# Input file
input_file = os.path.join(
    project_dir,
    "07_explainability/shap/substructure_mapping/TNIK_SHAP_substructure_mapping.csv"
)

# Output directory
output_dir = os.path.join(
    project_dir,
    "07_explainability/shap/substructure_mapping"
)

os.makedirs(output_dir, exist_ok=True)

# Load mapping data
df = pd.read_csv(input_file)

print("Total mapping records:", len(df))

# Remove records without substructure information
df = df.dropna(subset=["Substructure_SMILES"])

# Summarize recurring substructures
summary = (
    df.groupby(
        ["Fingerprint_Bit", "Substructure_SMILES"],
        as_index=False
    )
    .agg(
        Compound_Count=("Compound_ID", "nunique"),
        Median_pIC50=("pIC50", "median"),
        Mean_pIC50=("pIC50", "mean"),
        Occurrence_Count=("Compound_ID", "size")
    )
)

# Sort within each fingerprint bit
summary = summary.sort_values(
    ["Fingerprint_Bit", "Compound_Count"],
    ascending=[True, False]
)

# Save complete summary
summary_file = os.path.join(
    output_dir,
    "TNIK_SHAP_substructure_summary.csv"
)

summary.to_csv(summary_file, index=False)

# Extract the most frequent substructures for each bit
top_summary = (
    summary.groupby("Fingerprint_Bit", group_keys=False)
    .head(5)
)

top_file = os.path.join(
    output_dir,
    "TNIK_SHAP_top_substructures.csv"
)

top_summary.to_csv(top_file, index=False)

# Display results
print("\nTop recurring substructures for each fingerprint bit:")

for bit in summary["Fingerprint_Bit"].unique():

    print("\n", bit)

    bit_data = top_summary[
        top_summary["Fingerprint_Bit"] == bit
    ]

    print(
        bit_data[
            [
                "Substructure_SMILES",
                "Compound_Count",
                "Median_pIC50"
            ]
        ].to_string(index=False)
    )

print("\nSubstructure summary completed!")

print("\nFiles saved:")
print(summary_file)
print(top_file)
