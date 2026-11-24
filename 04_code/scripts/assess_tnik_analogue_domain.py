from pathlib import Path
import json

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

TRAIN_FILE = Path("05_models/chemberta/TNIK_ChemBERTa_train.csv")
CANDIDATE_FILE = Path(
    "08_analogues/filtered/TNIK_BRICS_all_scored_candidates.csv"
)

OUT = Path("08_analogues/filtered")
OUT.mkdir(parents=True, exist_ok=True)

train = pd.read_csv(TRAIN_FILE)
candidates = pd.read_csv(CANDIDATE_FILE)

generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, fpSize=2048
)

train_mols = [
    Chem.MolFromSmiles(str(s))
    for s in train["canonical_smiles"]
]
train_mols = [m for m in train_mols if m is not None]

train_fps = [generator.GetFingerprint(m) for m in train_mols]

# Leave-one-out nearest-neighbour similarity within the training set
loo_similarities = []

for i, fp in enumerate(train_fps):
    others = train_fps[:i] + train_fps[i + 1:]
    similarities = DataStructs.BulkTanimotoSimilarity(fp, others)
    loo_similarities.append(max(similarities))

# Use the 5th percentile as a data-derived screening threshold
threshold = float(np.percentile(loo_similarities, 5))

# Remove duplicate canonical structures before assessment
candidates = candidates.drop_duplicates(
    subset="canonical_smiles"
).copy()

# Recalculate candidate-to-training maximum similarity
candidate_max = []

for smi in candidates["canonical_smiles"]:
    mol = Chem.MolFromSmiles(str(smi))
    if mol is None:
        candidate_max.append(np.nan)
        continue

    fp = generator.GetFingerprint(mol)
    sims = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
    candidate_max.append(max(sims))

candidates["max_training_similarity"] = candidate_max
candidates["AD_threshold"] = threshold

candidates["applicability_domain"] = np.where(
    candidates["max_training_similarity"] >= threshold,
    "Inside",
    "Outside"
)

candidates = candidates.sort_values(
    "predicted_pIC50_RF", ascending=False
).reset_index(drop=True)

full_out = OUT / "TNIK_BRICS_applicability_domain.csv"
top_out = OUT / "TNIK_BRICS_top100_AD_assessed.csv"
threshold_out = OUT / "TNIK_BRICS_AD_threshold.json"

candidates.to_csv(full_out, index=False)
candidates.head(100).to_csv(top_out, index=False)

summary = {
    "training_compounds_used": len(train_fps),
    "unique_candidates_assessed": len(candidates),
    "training_LOO_similarity_5th_percentile_threshold": threshold,
    "inside_domain": int((candidates["applicability_domain"] == "Inside").sum()),
    "outside_domain": int((candidates["applicability_domain"] == "Outside").sum())
}

threshold_out.write_text(json.dumps(summary, indent=2))

print("\n===== ANALOGUE APPLICABILITY DOMAIN =====")
for key, value in summary.items():
    print(f"{key}: {value}")

print("\nTop 10 candidates with domain assessment:")
print(
    candidates[
        [
            "candidate_id",
            "predicted_pIC50_RF",
            "max_training_similarity",
            "applicability_domain",
            "QED"
        ]
    ].head(10).round(3).to_string(index=False)
)

print("\nSaved:", full_out)
print("Saved:", top_out)
print("Saved:", threshold_out)
