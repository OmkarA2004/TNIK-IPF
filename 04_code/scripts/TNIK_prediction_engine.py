#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, Descriptors, rdMolDescriptors


PROJECT = Path.home() / "Desktop" / "TNIK_IPF_Project"

TRAIN_FILE = PROJECT / "05_models/splits/TNIK_train.csv"

MODEL_FILES = {
    "Random Forest": PROJECT / "05_models/rdkit_rf/TNIK_RandomForest_baseline.pkl",
    "Extra Trees": PROJECT / "05_models/additional_ml/TNIK_extra_trees.pkl",
    "Gradient Boosting": PROJECT / "05_models/additional_ml/TNIK_gradient_boosting.pkl",
    "SVR": PROJECT / "05_models/additional_ml/TNIK_svr.pkl",
}

OUTPUT_DIR = PROJECT / "05_models/predictions"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Match the Morgan fingerprint settings used in the existing
# analogue-domain and SHAP-related project scripts.
GENERATOR = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)


def load_training_data():
    if not TRAIN_FILE.exists():
        raise FileNotFoundError(f"Training data not found: {TRAIN_FILE}")

    train = pd.read_csv(TRAIN_FILE)

    excluded = {
        "molecule_chembl_id",
        "pIC50",
        "measurement_count"
    }

    feature_columns = [
        c for c in train.columns
        if c not in excluded and c.startswith("FP_")
    ]

    feature_columns = sorted(
        feature_columns,
        key=lambda x: int(x.split("_")[1])
    )

    if len(feature_columns) != 2048:
        raise ValueError(
            f"Expected 2048 fingerprint columns; found {len(feature_columns)}"
        )

    X_train = train[feature_columns].to_numpy(dtype=np.uint8)

    if X_train.shape[0] < 2:
        raise ValueError("At least two training compounds are required.")

    return train, feature_columns, X_train


def to_bitvector(array):
    bit_string = "".join(str(int(x)) for x in array)
    return DataStructs.CreateFromBitString(bit_string)


def calculate_ad_threshold(training_fps):
    """
    Calculate the 5th percentile of leave-one-out maximum
    Tanimoto similarities within the Random Forest training set.
    """
    loo_similarities = []

    for i, fp in enumerate(training_fps):
        others = training_fps[:i] + training_fps[i + 1:]

        similarities = DataStructs.BulkTanimotoSimilarity(
            fp, others
        )

        loo_similarities.append(max(similarities))

    threshold = float(np.percentile(loo_similarities, 5))

    return threshold, loo_similarities


def load_models():
    models = {}

    for name, path in MODEL_FILES.items():
        if not path.exists():
            print(f"WARNING: {name} model not found: {path}")
            continue

        models[name] = joblib.load(path)

    if "Random Forest" not in models:
        raise FileNotFoundError(
            "The primary Random Forest model could not be loaded."
        )

    return models


def predict(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(
            "Invalid SMILES. Please provide a valid chemical SMILES string."
        )

    canonical_smiles = Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True
    )

    fingerprint = GENERATOR.GetFingerprint(mol)

    query_array = np.zeros((2048,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fingerprint, query_array)

    # Query fingerprint with the same named features used during training.

    train, feature_columns, X_train = load_training_data()

    X_query = pd.DataFrame([query_array], columns=feature_columns)

    training_fps = [
        to_bitvector(row)
        for row in X_train
    ]

    threshold, loo_similarities = calculate_ad_threshold(
        training_fps
    )

    similarities = DataStructs.BulkTanimotoSimilarity(
        fingerprint,
        training_fps
    )

    max_similarity = float(max(similarities))
    nearest_index = int(np.argmax(similarities))

    nearest_compound = str(
        train.iloc[nearest_index]["molecule_chembl_id"]
    )

    if max_similarity >= threshold:
        ad_status = "Inside"
        ad_interpretation = (
            "The compound meets the training-derived similarity threshold."
        )
    else:
        ad_status = "Outside"
        ad_interpretation = (
            "The compound is below the training-derived similarity threshold. "
            "Its model prediction should be treated cautiously."
        )

    models = load_models()

    predictions = {}

    for name, model in models.items():
        value = float(model.predict(X_query)[0])
        predictions[name] = value

    rf_prediction = predictions["Random Forest"]

    # pIC50 = -log10(IC50 in mol/L)
    # IC50 in nM = 10 ** (9 - pIC50)
    predicted_ic50_nm = float(10 ** (9 - rf_prediction))

    values = list(predictions.values())

    model_summary = {
        "number_of_models": len(values),
        "prediction_mean_pIC50": float(np.mean(values)),
        "prediction_standard_deviation_pIC50": float(np.std(values)),
        "prediction_minimum_pIC50": float(np.min(values)),
        "prediction_maximum_pIC50": float(np.max(values)),
        "prediction_range_pIC50": float(np.max(values) - np.min(values))
    }

    descriptors = {
        "molecular_weight": float(Descriptors.MolWt(mol)),
        "logP": float(Descriptors.MolLogP(mol)),
        "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
        "H_bond_donors": int(rdMolDescriptors.CalcNumHBD(mol)),
        "H_bond_acceptors": int(rdMolDescriptors.CalcNumHBA(mol)),
        "rotatable_bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol))
    }

    result = {
        "compound": {
            "input_smiles": smiles,
            "canonical_smiles": canonical_smiles
        },
        "primary_model": "Random Forest",
        "predicted_pIC50": rf_prediction,
        "predicted_IC50_nM": predicted_ic50_nm,
        "all_model_predictions_pIC50": predictions,
        "model_disagreement": model_summary,
        "applicability_domain": {
            "status": ad_status,
            "maximum_training_similarity": max_similarity,
            "threshold_method": (
                "5th percentile of leave-one-out maximum "
                "Tanimoto similarities in the RF training set"
            ),
            "threshold": threshold,
            "nearest_training_compound": nearest_compound,
            "interpretation": ad_interpretation
        },
        "training_information": {
            "training_compounds": int(len(train)),
            "fingerprint_radius": 2,
            "fingerprint_bits": 2048
        },
        "molecular_descriptors": descriptors,
        "scientific_caution": (
            "This is an in-silico prediction, not experimental confirmation. "
            "Model disagreement is descriptive and is not a calibrated "
            "probability of prediction accuracy. Docking, binding affinity, "
            "functional inhibition and therapeutic efficacy are distinct "
            "evidence types."
        )
    }

    return result


def main():
    parser = argparse.ArgumentParser(
        description="TNIK machine-learning compound prediction engine"
    )

    parser.add_argument(
        "--smiles",
        required=True,
        help="Chemical SMILES string to evaluate"
    )

    parser.add_argument(
        "--json-stdout",
        action="store_true",
        help="Print only prediction JSON for website integration"
    )

    args = parser.parse_args()
    result = predict(args.smiles)
    output_json = json.dumps(result, indent=2)

    if args.json_stdout:
        print(output_json)
        return

    output_file = OUTPUT_DIR / "TNIK_latest_prediction.json"
    output_file.write_text(output_json, encoding="utf-8")

    print("\n========== TNIK PREDICTION RESULT ==========")
    print("Canonical SMILES:")
    print(result["compound"]["canonical_smiles"])

    print("\nPrimary model:")
    print(result["primary_model"])

    print(f"Predicted pIC50: {result['predicted_pIC50']:.3f}")
    print(f"Predicted IC50: {result['predicted_IC50_nM']:.3f} nM")

    print("\nModel predictions:")
    for name, value in result["all_model_predictions_pIC50"].items():
        print(f"{name}: {value:.3f}")

    ad = result["applicability_domain"]
    print("\nApplicability domain:")
    print("Status:", ad["status"])
    print("Maximum similarity:", round(ad["maximum_training_similarity"], 4))
    print("Data-derived threshold:", round(ad["threshold"], 4))
    print("Nearest training compound:", ad["nearest_training_compound"])

    print("\nModel disagreement:")
    print(json.dumps(result["model_disagreement"], indent=2))

    print("\nSaved result:")
    print(output_file)


if __name__ == "__main__":
    main()
