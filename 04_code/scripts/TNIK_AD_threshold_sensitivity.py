
import os
import pandas as pd
import matplotlib.pyplot as plt

# Project directory
project_dir = os.path.expanduser("~/Desktop/TNIK_IPF_Project")

# Input file
input_file = os.path.join(
    project_dir,
    "07_explainability/applicability_domain/TNIK_applicability_domain_results.csv"
)

# Output locations
table_dir = os.path.join(project_dir, "09_results/tables")
figure_dir = os.path.join(project_dir, "09_results/figures")

os.makedirs(table_dir, exist_ok=True)
os.makedirs(figure_dir, exist_ok=True)

# Load Applicability Domain results
df = pd.read_csv(input_file)

# Thresholds to evaluate
thresholds = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55]

results = []

for threshold in thresholds:

    inside = df[
        df["Maximum_Tanimoto_Similarity"] >= threshold
    ]

    outside = df[
        df["Maximum_Tanimoto_Similarity"] < threshold
    ]

    results.append({
        "Threshold": threshold,
        "Inside_Count": len(inside),
        "Outside_Count": len(outside),
        "Inside_Mean_Error": inside["Absolute_Error"].mean(),
        "Outside_Mean_Error": outside["Absolute_Error"].mean(),
        "Inside_Median_Error": inside["Absolute_Error"].median(),
        "Outside_Median_Error": outside["Absolute_Error"].median()
    })

# Create summary table
results_df = pd.DataFrame(results)

print("\nAPPLICABILITY DOMAIN THRESHOLD SENSITIVITY")
print(results_df.round(4).to_string(index=False))

# Save results
output_file = os.path.join(
    table_dir,
    "TNIK_AD_threshold_sensitivity.csv"
)

results_df.to_csv(output_file, index=False)

# Plot: number of inside and outside compounds
plt.figure(figsize=(10, 6))

plt.plot(
    results_df["Threshold"],
    results_df["Inside_Count"],
    marker="o",
    label="Inside AD"
)

plt.plot(
    results_df["Threshold"],
    results_df["Outside_Count"],
    marker="s",
    label="Outside AD"
)

plt.xlabel("Tanimoto Similarity Threshold")
plt.ylabel("Number of Test Compounds")
plt.title("Applicability Domain Threshold Sensitivity")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plot_file = os.path.join(
    figure_dir,
    "TNIK_AD_threshold_sensitivity.png"
)

plt.savefig(plot_file, dpi=300, bbox_inches="tight")
plt.show()
plt.close()

print("\nFiles saved successfully:")
print(output_file)
print(plot_file)
