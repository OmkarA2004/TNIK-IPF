
import os
import joblib
import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt

# Project directory
project_dir = os.path.expanduser("~/Desktop/TNIK_IPF_Project")

# Input files
model_file = os.path.join(
    project_dir,
    "05_models/rdkit_rf/TNIK_RandomForest_baseline.pkl"
)

test_file = os.path.join(
    project_dir,
    "05_models/splits/TNIK_test.csv"
)

# Output directories
shap_dir = os.path.join(
    project_dir,
    "07_explainability/shap"
)

fig_dir = os.path.join(
    project_dir,
    "09_results/figures"
)

os.makedirs(shap_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

# Load model and test data
print("Loading model and test data...")

model = joblib.load(model_file)
test_df = pd.read_csv(test_file)

excluded_columns = [
    "molecule_chembl_id",
    "pIC50",
    "measurement_count"
]

feature_columns = [
    col for col in test_df.columns
    if col not in excluded_columns
]

X_test = test_df[feature_columns]

print("Test compounds:", len(X_test))
print("Fingerprint features:", len(feature_columns))

# Calculate SHAP values
print("\nCalculating SHAP values...")

explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)

# Save SHAP values
shap_df = pd.DataFrame(
    shap_values,
    columns=feature_columns
)

shap_df.insert(
    0,
    "molecule_chembl_id",
    test_df["molecule_chembl_id"].values
)

shap_file = os.path.join(
    shap_dir,
    "TNIK_SHAP_values.csv"
)

shap_df.to_csv(shap_file, index=False)

# Calculate mean absolute SHAP importance
importance = pd.DataFrame({
    "Fingerprint_Bit": feature_columns,
    "Mean_Absolute_SHAP": np.abs(shap_values).mean(axis=0)
})

importance = importance.sort_values(
    "Mean_Absolute_SHAP",
    ascending=False
)

importance_file = os.path.join(
    shap_dir,
    "TNIK_SHAP_feature_importance.csv"
)

importance.to_csv(importance_file, index=False)

# Save SHAP summary plot
plt.figure()
shap.summary_plot(
    shap_values,
    X_test,
    max_display=20,
    show=False
)

plt.tight_layout()

summary_file = os.path.join(
    fig_dir,
    "TNIK_SHAP_summary.png"
)

plt.savefig(summary_file, dpi=300, bbox_inches="tight")
plt.close()

# Save SHAP bar plot
plt.figure()
shap.summary_plot(
    shap_values,
    X_test,
    plot_type="bar",
    max_display=20,
    show=False
)

plt.tight_layout()

bar_file = os.path.join(
    fig_dir,
    "TNIK_SHAP_bar.png"
)

plt.savefig(bar_file, dpi=300, bbox_inches="tight")
plt.close()

print("\nSHAP analysis completed!")

print("\nTop 10 important fingerprint bits:")
print(importance.head(10).to_string(index=False))

print("\nSaved files:")
print(shap_file)
print(importance_file)
print(summary_file)
print(bar_file)
