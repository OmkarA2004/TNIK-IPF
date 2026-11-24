from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from rdkit import DataStructs

from src.compound_intelligence import identify_compound
from src.molecular_features import smiles_to_fingerprint

ROOT = Path(__file__).resolve().parents[1]

MODEL_FILES = {
    "Random Forest": ROOT / "models/trained_models/TNIK_RandomForest_baseline.pkl",
    "SVR": ROOT / "models/trained_models/TNIK_svr.pkl",
    "Gradient Boosting": ROOT / "models/trained_models/TNIK_gradient_boosting.pkl",
    "Extra Trees": ROOT / "models/trained_models/TNIK_extra_trees.pkl",
}

REFERENCE_FILE = (
    ROOT / "data/reference_compounds/TNIK_Morgan_Fingerprints_unique.csv"
)

FEATURE_COLUMNS = [f"FP_{i}" for i in range(2048)]

_models = None
_reference = None


def load_models():
    global _models

    if _models is None:
        _models = {
            name: joblib.load(path)
            for name, path in MODEL_FILES.items()
        }

    return _models


def load_reference_data():
    global _reference

    if _reference is None:
        _reference = pd.read_csv(REFERENCE_FILE)

        missing = set(FEATURE_COLUMNS) - set(_reference.columns)

        if missing:
            raise ValueError("Reference fingerprint columns are missing.")

    return _reference


def calculate_similarity(features):
    reference = load_reference_data()

    query_fp = np.asarray(features, dtype=np.uint8)
    query = DataStructs.CreateFromBitString(
        "".join(query_fp.astype(str))
    )

    similarities = []

    for _, row in reference.iterrows():
        bits = row[FEATURE_COLUMNS].to_numpy(dtype=np.uint8)
        reference_fp = DataStructs.CreateFromBitString(
            "".join(bits.astype(str))
        )
        similarities.append(
            DataStructs.TanimotoSimilarity(query, reference_fp)
        )

    best_index = int(np.argmax(similarities))

    return {
        "maximum_similarity": float(similarities[best_index]),
        "closest_compound": reference.iloc[best_index][
            "molecule_chembl_id"
        ],
    }


def predict_compound(query):
    compound = identify_compound(query)

    # Generate the same type of Morgan fingerprint used by the models.
    _, features = smiles_to_fingerprint(
        compound["canonical_smiles"]
    )

    X = pd.DataFrame([features], columns=FEATURE_COLUMNS)

    predictions = {}

    for name, model in load_models().items():
        predicted_pic50 = float(model.predict(X)[0])

        predictions[name] = {
            "predicted_pIC50": round(predicted_pic50, 3),
            "estimated_IC50_nM": round(
                10 ** (9 - predicted_pic50), 2
            ),
        }

    similarity = calculate_similarity(features)

    return {
        "compound": compound,
        "predictions": predictions,
        "similarity": similarity,
    }
