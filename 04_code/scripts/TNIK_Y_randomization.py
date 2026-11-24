import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_validate

# File paths
input_file = "../../03_processed_data/activity_clean/TNIK_Morgan_Fingerprints_unique.csv"

output_file = "../../05_models/metrics/TNIK_Y_randomization_results.csv"

# Load dataset
df = pd.read_csv(input_file)

# Separate features and target
exclude = [
    "molecule_chembl_id",
    "pIC50",
    "measurement_count"
]

feature_columns = [c for c in df.columns if c not in exclude]

X = df[feature_columns]
y = df["pIC50"].values

# Cross-validation setup
cv = KFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

results = []

# Run 10 Y-randomization experiments
for i in range(10):

    print(f"Running randomization {i+1}/10")

    # Shuffle activity values
    y_shuffled = np.random.RandomState(i).permutation(y)

    model = RandomForestRegressor(
        n_estimators=500,
        random_state=42,
        n_jobs=-1
    )

    scores = cross_validate(
        model,
        X,
        y_shuffled,
        cv=cv,
        scoring={
            "R2": "r2",
            "RMSE": "neg_root_mean_squared_error",
            "MAE": "neg_mean_absolute_error"
        }
    )

    results.append({
        "Randomization": i+1,
        "Mean_R2": scores["test_R2"].mean(),
        "Mean_RMSE": -scores["test_RMSE"].mean(),
        "Mean_MAE": -scores["test_MAE"].mean()
    })

# Save results
results_df = pd.DataFrame(results)

results_df.to_csv(output_file, index=False)

print("\nY-RANDOMIZATION RESULTS")
print(results_df.to_string(index=False))

print("\nAVERAGE RANDOMIZED RESULTS")
print("Mean R2:", round(results_df["Mean_R2"].mean(), 4))
print("Mean RMSE:", round(results_df["Mean_RMSE"].mean(), 4))
print("Mean MAE:", round(results_df["Mean_MAE"].mean(), 4))

print("\nSaved:", output_file)


