
import os
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_val_score

# Project paths
project_dir = os.path.expanduser("~/Desktop/TNIK_IPF_Project")

input_file = os.path.join(
    project_dir,
    "03_processed_data/activity_clean/TNIK_Morgan_Fingerprints_unique.csv"
)

output_file = os.path.join(
    project_dir,
    "05_models/metrics/TNIK_Y_randomization_100_results.csv"
)

# Load dataset
df = pd.read_csv(input_file)

# Select fingerprint features
feature_columns = [
    col for col in df.columns
    if col not in ["molecule_chembl_id", "pIC50", "measurement_count"]
]

X = df[feature_columns]
y = df["pIC50"].values

# Cross-validation setup
cv = KFold(n_splits=5, shuffle=True, random_state=42)

results = []

print("Starting 100 Y-randomization tests...")
print("This may take several minutes.")

# Perform 100 random permutations
for i in range(100):

    y_random = np.random.RandomState(i).permutation(y)

    model = RandomForestRegressor(
        n_estimators=500,
        random_state=42,
        n_jobs=-1
    )

    scores = cross_val_score(
        model,
        X,
        y_random,
        cv=cv,
        scoring="r2",
        n_jobs=1
    )

    results.append({
        "Randomization": i + 1,
        "Mean_R2": scores.mean()
    })

    if (i + 1) % 10 == 0:
        print(f"Completed {i + 1}/100 permutations")

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, index=False)

print("\nY-randomization completed!")
print("Results saved to:", output_file)

print("\nSummary:")
print(results_df["Mean_R2"].describe())

print("\nOriginal model mean random CV R2: 0.5439")
print("Average randomized R2:", round(results_df["Mean_R2"].mean(), 4))
