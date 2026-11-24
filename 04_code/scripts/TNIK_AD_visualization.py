
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

# Output directory
output_dir = os.path.join(
    project_dir,
    "09_results/figures"
)

os.makedirs(output_dir, exist_ok=True)

# Load data
df = pd.read_csv(input_file)

threshold = 0.40

# Plot 1: Similarity of each test compound
plt.figure(figsize=(15, 7))

plt.scatter(
    range(1, len(df) + 1),
    df["Maximum_Tanimoto_Similarity"],
    s=45,
    alpha=0.8
)

for i, row in df.reset_index(drop=True).iterrows():
    plt.annotate(
        str(row["molecule_chembl_id"]),
        (i + 1, row["Maximum_Tanimoto_Similarity"]),
        xytext=(3, 5),
        textcoords="offset points",
        fontsize=5,
        rotation=45
    )

plt.axhline(
    threshold,
    linestyle="--",
    linewidth=2,
    label="AD threshold = 0.40"
)

plt.xlabel("Test compound number")
plt.ylabel("Maximum Tanimoto similarity")
plt.title("Applicability Domain: Test Compound Similarity")
plt.ylim(0, 1.05)
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()

similarity_file = os.path.join(
    output_dir,
    "TNIK_AD_similarity_plot.png"
)

plt.savefig(similarity_file, dpi=300, bbox_inches="tight")
plt.show()
plt.close()

# Plot 2: Similarity versus prediction error
plt.figure(figsize=(10, 7))

plt.scatter(
    df["Maximum_Tanimoto_Similarity"],
    df["Absolute_Error"],
    s=55,
    alpha=0.8
)

for _, row in df.iterrows():
    plt.annotate(
        str(row["molecule_chembl_id"]),
        (
            row["Maximum_Tanimoto_Similarity"],
            row["Absolute_Error"]
        ),
        xytext=(3, 4),
        textcoords="offset points",
        fontsize=6
    )

plt.axvline(
    threshold,
    linestyle="--",
    linewidth=2,
    label="AD threshold = 0.40"
)

plt.xlabel("Maximum Tanimoto similarity")
plt.ylabel("Absolute prediction error")
plt.title("Applicability Domain versus Prediction Error")
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()

error_file = os.path.join(
    output_dir,
    "TNIK_AD_similarity_vs_error.png"
)

plt.savefig(error_file, dpi=300, bbox_inches="tight")
plt.show()
plt.close()

print("\nApplicability Domain plots saved:")
print(similarity_file)
print(error_file)
