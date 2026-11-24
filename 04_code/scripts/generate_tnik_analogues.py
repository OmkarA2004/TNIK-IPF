from pathlib import Path
import random
import joblib

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import BRICS, AllChem, Descriptors, Crippen, Lipinski, QED, rdMolDescriptors

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

TRAIN_FILE = Path("05_models/chemberta/TNIK_ChemBERTa_train.csv")
FP_TRAIN_FILE = Path("05_models/splits/TNIK_train.csv")
MODEL_FILE = Path("05_models/rdkit_rf/TNIK_RandomForest_baseline.pkl")

OUT = Path("08_analogues")
GENERATED = OUT / "generated"
FILTERED = OUT / "filtered"
GENERATED.mkdir(parents=True, exist_ok=True)
FILTERED.mkdir(parents=True, exist_ok=True)

# Load parent compounds and fingerprint training data
parents = pd.read_csv(TRAIN_FILE)
fp_train = pd.read_csv(FP_TRAIN_FILE)

# Select the top 25% of training compounds by experimental pIC50
cutoff = parents["pIC50"].quantile(0.75)
active_parents = parents[parents["pIC50"] >= cutoff].copy()
active_parents = active_parents.sort_values("pIC50", ascending=False)

print("Selected active parent compounds:", len(active_parents))
print("Parent pIC50 threshold:", round(cutoff, 3))

# Parse parent structures and collect BRICS fragments
parent_mols = []
fragment_smiles = set()

for smi in active_parents["canonical_smiles"]:
    mol = Chem.MolFromSmiles(str(smi))
    if mol is None:
        continue
    parent_mols.append(mol)
    fragment_smiles.update(BRICS.BRICSDecompose(mol))

fragments = []
for smi in sorted(fragment_smiles):
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        fragments.append(mol)

print("Usable parent molecules:", len(parent_mols))
print("Unique BRICS fragments:", len(fragments))

if len(fragments) < 2:
    raise RuntimeError("Not enough BRICS fragments to generate analogues.")

# Recombine fragments with a limited depth
candidate_smiles = set()
builder = BRICS.BRICSBuild(
    fragments,
    onlyCompleteMols=True,
    uniquify=True,
    scrambleReagents=False,
    maxDepth=2
)

for i, mol in enumerate(builder):
    if mol is None:
        continue
    try:
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol, canonical=True)
        if smi:
            candidate_smiles.add(smi)
    except Exception:
        pass

    if len(candidate_smiles) >= 3000:
        break

print("Unique generated structures:", len(candidate_smiles))

# Prepare training molecules and Morgan fingerprints
training_smiles = parents["canonical_smiles"].dropna().astype(str).tolist()
training_mols = [
    Chem.MolFromSmiles(s) for s in training_smiles
]
training_mols = [m for m in training_mols if m is not None]

training_fps = [
    AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    for m in training_mols
]

# Load existing RF model
model = joblib.load(MODEL_FILE)

# Identify the exact fingerprint columns used by the saved training data
excluded = {"molecule_chembl_id", "pIC50", "canonical_smiles"}
feature_cols = [
    c for c in fp_train.columns
    if c not in excluded
    and pd.api.types.is_numeric_dtype(fp_train[c])
    and set(fp_train[c].dropna().unique()).issubset({0, 1})
]

if len(feature_cols) != 2048:
    raise ValueError(
        f"Expected 2048 RF features, found {len(feature_cols)}. "
        "Stopping to protect model input consistency."
    )

# Filter candidates and calculate properties
records = []

for smi in candidate_smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    qed = QED.qed(mol)

    # Basic property filters, not proof of safety or oral bioavailability
    if not (150 <= mw <= 650):
        continue
    if not (-1 <= logp <= 6):
        continue
    if tpsa > 140 or hbd > 5 or hba > 10:
        continue
    if qed < 0.35:
        continue

    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, 2, nBits=2048
    )

    max_similarity = max(
        DataStructs.BulkTanimotoSimilarity(fp, training_fps)
    )

    records.append({
        "canonical_smiles": smi,
        "molecular_weight": mw,
        "logP": logp,
        "TPSA": tpsa,
        "HBD": hbd,
        "HBA": hba,
        "QED": qed,
        "max_training_similarity": max_similarity,
        "fingerprint": list(fp)
    })

print("Candidates passing property filters:", len(records))

if not records:
    raise RuntimeError(
        "No candidates passed the filters. No model predictions were made."
    )

# Predict using the existing Random Forest
X = pd.DataFrame(
    [r["fingerprint"] for r in records],
    columns=feature_cols
)

predictions = model.predict(X)

for record, prediction in zip(records, predictions):
    record["predicted_pIC50_RF"] = float(prediction)
    record.pop("fingerprint")

result = pd.DataFrame(records)
result = result.sort_values(
    "predicted_pIC50_RF", ascending=False
).reset_index(drop=True)

result.insert(
    0, "candidate_id",
    [f"TNIK_BRICS_{i:04d}" for i in range(1, len(result) + 1)]
)

# Save all filtered candidates and the top 100 ranked candidates
all_file = FILTERED / "TNIK_BRICS_all_scored_candidates.csv"
top_file = FILTERED / "TNIK_BRICS_top100_candidates.csv"

result.to_csv(all_file, index=False)
result.head(100).to_csv(top_file, index=False)

print("\n===== ANALOGUE GENERATION SUMMARY =====")
print("Generated unique structures:", len(candidate_smiles))
print("Property-filtered and scored:", len(result))
print("Top candidates saved:", min(100, len(result)))
print("\nTop 10 candidates:")
print(
    result[
        [
            "candidate_id",
            "predicted_pIC50_RF",
            "max_training_similarity",
            "QED",
            "molecular_weight"
        ]
    ].head(10).round(3).to_string(index=False)
)

print("\nSaved all candidates:", all_file)
print("Saved top candidates:", top_file)
print("\nAnalogue generation and scoring completed.")
