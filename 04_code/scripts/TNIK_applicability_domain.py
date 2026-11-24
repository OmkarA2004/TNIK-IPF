
import os
import joblib
import numpy as np
import pandas as pd

from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# Project directory
project_dir = os.path.expanduser("~/Desktop/TNIK_IPF_Project")

# Input files
train_file = os.path.join(
    project_dir,
    "05_models/splits/TNIK_train.csv"
)

test_file = os.path.join(
    project_dir,
    "05_models/splits/TNIK_test.csv"
)

prediction_file = os.path.join(
    project_dir,
    "05_models/predictions/TNIK_RF_test_predictions.csv"
)

# Output directory
output_dir = os.path.join(
    project_dir,
    "07_explainability/applicability_domain"
)

os.makedirs(output_dir, exist_ok=True)

# Load datasets
train_df = pd.read_csv(train_file)
test_df = pd.read_csv(test_file)
pred_df = pd.read_csv(prediction_file)

# Identify fingerprint columns
excluded = [
    "molecule_chembl_id",
    "pIC50",
    "measurement_count"
]

feature_columns = [
    col for col in train_df.columns
    if col not in excluded
]

X_train = train_df[feature_columns].values
X_test = test_df[feature_columns].values

# Convert binary fingerprint arrays to RDKit bit vectors
def to_bitvector(array):
    return DataStructs.CreateFromBitString(
        "".join(str(int(x)) for x in array)
    )

train_fps = [to_bitvector(row) for row in X_train]
test_fps = [to_bitvector(row) for row in X_test]

# Calculate nearest-neighbour Tanimoto similarity
print("Calculating nearest-neighbour similarities...")

max_similarities = []

for test_fp in test_fps:

    similarities = DataStructs.BulkTanimotoSimilarity(
        test_fp,
        train_fps
    )

    max_similarities.append(max(similarities))

# Create results table
ad_df = pd.DataFrame({
    "molecule_chembl_id": test_df["molecule_chembl_id"],
    "Actual_pIC50": test_df["pIC50"],
    "Maximum_Tanimoto_Similarity": max_similarities
})

# Add predictions
ad_df = ad_df.merge(
    pred_df[
        ["molecule_chembl_id", "Predicted_pIC50", "Absolute_Error"]
    ],
    on="molecule_chembl_id",
    how="left"
)

# Provisional threshold
threshold = 0.40

ad_df["AD_Status"] = np.where(
    ad_df["Maximum_Tanimoto_Similarity"] >= threshold,
    "Inside",
    "Outside"
)

# Save results
output_file = os.path.join(
    output_dir,
    "TNIK_applicability_domain_results.csv"
)

ad_df.to_csv(output_file, index=False)

# Summary statistics
print("\nApplicability Domain Summary")

print("Total test compounds:", len(ad_df))

print("\nAD classification:")
print(ad_df["AD_Status"].value_counts())

print("\nMaximum similarity statistics:")
print(ad_df["Maximum_Tanimoto_Similarity"].describe())

print("\nPrediction error by AD status:")
print(
    ad_df.groupby("AD_Status")["Absolute_Error"]
    .agg(["count", "mean", "median"])
)

print("\nResults saved to:")
print(output_file)
