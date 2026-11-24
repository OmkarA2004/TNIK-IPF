from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "08_analogues/filtered/TNIK_analogues_filtered.csv"
OUT = ROOT / "08_analogues/admet/submission_files"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SOURCE)

# Keep one row per unique valid SMILES.
df = df.drop_duplicates(subset="smiles").copy()
df = df.dropna(subset=["smiles"])

# SwissADME accepts one SMILES and optional name per line.
swiss_file = OUT / "SwissADME_105_molecules.smi"

with open(swiss_file, "w") as f:
    for i, row in df.iterrows():
        name = f"{row['parent_compound']}_analogue_{i+1}"
        f.write(f"{row['smiles']} {name}\n")

# pkCSM: split into batches of at most 100.
for batch_number, start in enumerate(range(0, len(df), 100), 1):
    batch = df.iloc[start:start+100].copy()
    batch = batch.rename(columns={"smiles": "SMILES"})
    batch[["SMILES"]].to_csv(
        OUT / f"pkCSM_batch_{batch_number}.csv",
        index=False
    )

print("\n========== ADMET SUBMISSION FILES ==========")
print("Unique molecules:", len(df))
print("SwissADME file:", swiss_file)
print("pkCSM batches:", (len(df) + 99) // 100)

for file in sorted(OUT.iterdir()):
    print(file.name, "-", file.stat().st_size, "bytes")
