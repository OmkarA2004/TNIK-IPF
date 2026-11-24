from pathlib import Path
import json
import random

# Import scikit-learn first to initialize the shared numerical libraries
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from transformers import AutoTokenizer, AutoModelForSequenceClassification

# ---------------- CONFIGURATION ----------------
SEED = 42
MODEL_NAME = "DeepChem/ChemBERTa-77M-MLM"
MAX_LENGTH = 128
BATCH_SIZE = 8
EPOCHS = 20
LEARNING_RATE = 2e-5
PATIENCE = 4

OUT = Path("05_models/chemberta")
OUT.mkdir(parents=True, exist_ok=True)

TRAIN_FILE = OUT / "TNIK_ChemBERTa_train.csv"
TEST_FILE = OUT / "TNIK_ChemBERTa_test.csv"

# ---------------- REPRODUCIBILITY ----------------
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

if torch.backends.mps.is_available():
    device = torch.device("mps")
elif torch.cuda.is_available():
    device = torch.device("cuda")
else:
    device = torch.device("cpu")

print("Device:", device)
print("Loading datasets...")

train_df = pd.read_csv(TRAIN_FILE)
test_df = pd.read_csv(TEST_FILE)

required = {"molecule_chembl_id", "canonical_smiles", "pIC50"}
for name, df in [("train", train_df), ("test", test_df)]:
    if not required.issubset(df.columns):
        raise ValueError(f"Missing required columns in {name} dataset.")
    if df[list(required)].isna().any().any():
        raise ValueError(f"Missing values in {name} dataset.")

if set(train_df.molecule_chembl_id) & set(test_df.molecule_chembl_id):
    raise ValueError("Train/test overlap detected.")

# Internal validation split; external test set remains untouched.
fit_df, val_df = train_test_split(
    train_df, test_size=0.20, random_state=SEED
)

print("Training compounds:", len(fit_df))
print("Validation compounds:", len(val_df))
print("Final test compounds:", len(test_df))

# ---------------- TOKENIZER AND MODEL ----------------
print("\nLoading pretrained ChemBERTa:", MODEL_NAME)

tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)

model = AutoModelForSequenceClassification.from_pretrained(
    MODEL_NAME,
    num_labels=1,
    problem_type="regression"
)

model.to(device)

class SmilesDataset(Dataset):
    def __init__(self, frame):
        self.smiles = frame["canonical_smiles"].astype(str).tolist()
        self.labels = frame["pIC50"].astype(float).to_numpy()

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, index):
        encoded = tokenizer(
            self.smiles[index],
            max_length=MAX_LENGTH,
            padding="max_length",
            truncation=True,
            return_tensors="pt"
        )

        return {
            "input_ids": encoded["input_ids"].squeeze(0),
            "attention_mask": encoded["attention_mask"].squeeze(0),
            "labels": torch.tensor(self.labels[index], dtype=torch.float32)
        }

fit_loader = DataLoader(
    SmilesDataset(fit_df), batch_size=BATCH_SIZE, shuffle=True
)
val_loader = DataLoader(
    SmilesDataset(val_df), batch_size=BATCH_SIZE, shuffle=False
)
test_loader = DataLoader(
    SmilesDataset(test_df), batch_size=BATCH_SIZE, shuffle=False
)

optimizer = torch.optim.AdamW(
    model.parameters(), lr=LEARNING_RATE, weight_decay=0.01
)
loss_function = torch.nn.MSELoss()

def evaluate(loader):
    model.eval()
    actual, predicted = [], []

    with torch.no_grad():
        for batch in loader:
            inputs = {
                "input_ids": batch["input_ids"].to(device),
                "attention_mask": batch["attention_mask"].to(device)
            }
            labels = batch["labels"].to(device)
            outputs = model(**inputs).logits.squeeze(-1)

            actual.extend(labels.cpu().numpy().tolist())
            predicted.extend(outputs.cpu().numpy().tolist())

    rmse = float(np.sqrt(mean_squared_error(actual, predicted)))
    mae = float(mean_absolute_error(actual, predicted))
    r2 = float(r2_score(actual, predicted))

    return {"RMSE": rmse, "MAE": mae, "R2": r2}, actual, predicted

# ---------------- FINE-TUNING ----------------
best_val_rmse = float("inf")
best_epoch = 0
epochs_without_improvement = 0
history = []

best_dir = OUT / "best_model"

print("\nStarting ChemBERTa fine-tuning...\n")

for epoch in range(1, EPOCHS + 1):
    model.train()
    losses = []

    for batch in fit_loader:
        optimizer.zero_grad()

        inputs = {
            "input_ids": batch["input_ids"].to(device),
            "attention_mask": batch["attention_mask"].to(device)
        }
        labels = batch["labels"].to(device)

        predictions = model(**inputs).logits.squeeze(-1)
        loss = loss_function(predictions, labels)

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

        losses.append(loss.item())

    val_metrics, _, _ = evaluate(val_loader)
    train_loss = float(np.mean(losses))

    history.append({
        "epoch": epoch,
        "training_loss": train_loss,
        "validation_RMSE": val_metrics["RMSE"],
        "validation_MAE": val_metrics["MAE"],
        "validation_R2": val_metrics["R2"]
    })

    print(
        f"Epoch {epoch:02d} | "
        f"Loss {train_loss:.4f} | "
        f"Val RMSE {val_metrics['RMSE']:.4f} | "
        f"Val MAE {val_metrics['MAE']:.4f} | "
        f"Val R2 {val_metrics['R2']:.4f}"
    )

    if val_metrics["RMSE"] < best_val_rmse:
        best_val_rmse = val_metrics["RMSE"]
        best_epoch = epoch
        epochs_without_improvement = 0

        model.save_pretrained(best_dir)
        tokenizer.save_pretrained(best_dir)
    else:
        epochs_without_improvement += 1

    if epochs_without_improvement >= PATIENCE:
        print("\nEarly stopping: validation RMSE did not improve.")
        break

# ---------------- FINAL TEST EVALUATION ----------------
print("\nLoading best validation checkpoint...")
model = AutoModelForSequenceClassification.from_pretrained(best_dir)
model.to(device)

test_metrics, actual, predicted = evaluate(test_loader)

predictions_df = test_df[
    ["molecule_chembl_id", "canonical_smiles", "pIC50"]
].copy()

predictions_df["predicted_pIC50"] = predicted
predictions_df["absolute_error"] = np.abs(
    predictions_df["pIC50"] - predictions_df["predicted_pIC50"]
)

predictions_df.to_csv(
    OUT / "TNIK_ChemBERTa_test_predictions.csv", index=False
)

pd.DataFrame(history).to_csv(
    OUT / "TNIK_ChemBERTa_training_history.csv", index=False
)

results = {
    "model": MODEL_NAME,
    "seed": SEED,
    "training_compounds": len(fit_df),
    "validation_compounds": len(val_df),
    "test_compounds": len(test_df),
    "best_epoch": best_epoch,
    "best_validation_RMSE": best_val_rmse,
    "test_RMSE": test_metrics["RMSE"],
    "test_MAE": test_metrics["MAE"],
    "test_R2": test_metrics["R2"],
    "device": str(device)
}

with open(OUT / "TNIK_ChemBERTa_metrics.json", "w") as f:
    json.dump(results, f, indent=4)

print("\n========== CHEMBERTA FINAL TEST RESULTS ==========")
for key, value in results.items():
    print(f"{key}: {value}")

print("\nSaved model:", best_dir)
print("Saved predictions:", OUT / "TNIK_ChemBERTa_test_predictions.csv")
print("Saved training history:", OUT / "TNIK_ChemBERTa_training_history.csv")
print("Saved metrics:", OUT / "TNIK_ChemBERTa_metrics.json")
print("\nChemBERTa training completed.")
