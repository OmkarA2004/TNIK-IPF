from rdkit import Chem
from collections import Counter
from pathlib import Path
import pandas as pd

file = "02_raw_data/approved_drugs/e-Drug3D_2197.sdf"
supplier = Chem.SDMolSupplier(file, removeHs=False)

rows = []

for mol in supplier:
    if mol is None:
        continue

    rows.append({
        "name": mol.GetProp("name") if mol.HasProp("name") else "",
        "ID": mol.GetProp("ID") if mol.HasProp("ID") else "",
        "CAS": mol.GetProp("cas") if mol.HasProp("cas") else "",
        "status": mol.GetProp("status") if mol.HasProp("status") else "",
        "SMILES": Chem.MolToSmiles(Chem.RemoveHs(mol))
    })

df = pd.DataFrame(rows)

print("Valid molecules:", len(df))
print("Unique drug names:", df["name"].nunique())
print("Unique IDs:", df["ID"].nunique())
print("Unique CAS numbers:", df["CAS"].nunique())

print("\nStatus distribution:")
print(df["status"].value_counts().to_string())

print("\nFirst 15 drug records:")
print(df[["name", "ID", "CAS", "status"]].head(15).to_string(index=False))

out = Path("02_raw_data/approved_drugs/e-Drug3D_metadata_audit.csv")
df.to_csv(out, index=False)
print("\nSaved:", out)
