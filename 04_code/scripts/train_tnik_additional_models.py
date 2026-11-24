from pathlib import Path
import json
import joblib

# Import sklearn first to avoid the OpenMP import conflict
from sklearn.model_selection import KFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

import numpy as np
import pandas as pd

SEED = 42

train_path = Path("05_models/splits/TNIK_train.csv")
test_path = Path("05_models/splits/TNIK_test.csv")

out_dir = Path("05_models/additional_ml")
out_dir.mkdir(parents=True, exist_ok=True)

train = pd.read_csv(train_path)
test = pd.read_csv(test_path)

target = "pIC50"
id_col = "molecule_chembl_id"

# Identify the existing binary Morgan fingerprint columns
excluded = {target, id_col, "canonical_smiles"}
candidate_cols = [
    c for c in train.columns
    if c not in excluded and c in test.columns
]

feature_cols = []
for col in candidate_cols:
    a = pd.to_numeric(train[col], errors="coerce")
    b = pd.to_numeric(test[col], errors="coerce")
    values = set(a.dropna().unique()) | set(b.dropna().unique())

    if not a.isna().any() and not b.isna().any() and values.issubset({0, 1}):
        feature_cols.append(col)

if len(feature_cols) != 2048:
    raise ValueError(
        f"Expected 2048 Morgan fingerprint columns, found {len(feature_cols)}. "
        "Stopping rather than using an incorrect feature set."
    )

X_train = train[feature_cols].astype(np.float32)
X_test = test[feature_cols].astype(np.float32)
y_train = train[target].astype(float)
y_test = test[target].astype(float)

print("Training compounds:", len(train))
print("Test compounds:", len(test))
print("Morgan fingerprint features:", len(feature_cols))

models = {
    "SVR": Pipeline([
        ("scaler", StandardScaler()),
        ("model", SVR(kernel="rbf", C=10, epsilon=0.2, gamma="scale"))
    ]),
    "Extra Trees": ExtraTreesRegressor(
        n_estimators=500,
        min_samples_leaf=2,
        max_features=0.8,
        random_state=SEED,
        n_jobs=-1
    ),
    "Gradient Boosting": GradientBoostingRegressor(
        n_estimators=300,
        learning_rate=0.03,
        max_depth=2,
        loss="huber",
        random_state=SEED
    )
}

cv = KFold(n_splits=5, shuffle=True, random_state=SEED)
metrics = []
predictions = []

for name, model in models.items():
    print(f"\n===== {name} =====")

    cv_results = cross_validate(
        model,
        X_train,
        y_train,
        cv=cv,
        scoring={
            "r2": "r2",
            "rmse": "neg_root_mean_squared_error",
            "mae": "neg_mean_absolute_error"
        },
        n_jobs=1
    )

    cv_r2 = cv_results["test_r2"]
    cv_rmse = -cv_results["test_rmse"]
    cv_mae = -cv_results["test_mae"]

    model.fit(X_train, y_train)
    pred = model.predict(X_test)

    test_r2 = r2_score(y_test, pred)
    test_rmse = np.sqrt(mean_squared_error(y_test, pred))
    test_mae = mean_absolute_error(y_test, pred)

    print(f"CV R2:   {cv_r2.mean():.4f} ± {cv_r2.std():.4f}")
    print(f"CV RMSE: {cv_rmse.mean():.4f} ± {cv_rmse.std():.4f}")
    print(f"CV MAE:  {cv_mae.mean():.4f} ± {cv_mae.std():.4f}")
    print(f"Test R2:   {test_r2:.4f}")
    print(f"Test RMSE: {test_rmse:.4f}")
    print(f"Test MAE:  {test_mae:.4f}")

    metrics.append({
        "model": name,
        "CV_R2_mean": cv_r2.mean(),
        "CV_R2_std": cv_r2.std(),
        "CV_RMSE_mean": cv_rmse.mean(),
        "CV_RMSE_std": cv_rmse.std(),
        "CV_MAE_mean": cv_mae.mean(),
        "CV_MAE_std": cv_mae.std(),
        "test_R2": test_r2,
        "test_RMSE": test_rmse,
        "test_MAE": test_mae
    })

    for compound_id, actual, predicted in zip(
        test[id_col], y_test, pred
    ):
        predictions.append({
            "molecule_chembl_id": compound_id,
            "actual_pIC50": actual,
            "predicted_pIC50": predicted,
            "model": name
        })

    joblib.dump(model, out_dir / f"TNIK_{name.lower().replace(' ', '_')}.pkl")

# Save additional-model metrics
metrics_df = pd.DataFrame(metrics)
metrics_df.to_csv(out_dir / "TNIK_additional_ML_metrics.csv", index=False)

# Save predictions
pred_df = pd.DataFrame(predictions)
pred_df.to_csv(out_dir / "TNIK_additional_ML_predictions.csv", index=False)

# Combine with already completed models for one comparison table
comparison = [
    {
        "model": "Random Forest",
        "CV_R2_mean": 0.5439,
        "CV_R2_std": 0.0931,
        "test_R2": 0.4914,
        "test_RMSE": 0.7976,
        "test_MAE": 0.5917
    },
    {
        "model": "ChemBERTa",
        "CV_R2_mean": np.nan,
        "CV_R2_std": np.nan,
        "test_R2": 0.0924212723,
        "test_RMSE": 1.0655537048,
        "test_MAE": 0.8626029360
    }
]

for row in metrics:
    comparison.append({
        "model": row["model"],
        "CV_R2_mean": row["CV_R2_mean"],
        "CV_R2_std": row["CV_R2_std"],
        "test_R2": row["test_R2"],
        "test_RMSE": row["test_RMSE"],
        "test_MAE": row["test_MAE"]
    })

comparison_df = pd.DataFrame(comparison)
comparison_df.to_csv(
    "09_results/tables/TNIK_all_ML_models_comparison.csv",
    index=False
)

print("\n===== ALL MODEL COMPARISON =====")
print(comparison_df.round(4).to_string(index=False))

print("\nSaved additional model metrics:")
print(out_dir / "TNIK_additional_ML_metrics.csv")
print("Saved predictions:")
print(out_dir / "TNIK_additional_ML_predictions.csv")
print("Saved final comparison:")
print("09_results/tables/TNIK_all_ML_models_comparison.csv")
print("\nAdditional ML model training completed.")
