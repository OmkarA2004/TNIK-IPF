from rdkit import Chem
import pandas as pd
from pathlib import Path

source = "02_raw_data/approved_drugs/e-Drug3D_2197.sdf"
supplier = Chem.SDMolSupplier(source, removeHs=False)

rows = []

for mol in supplier:
    if mol is None:
        continue

    status = mol.GetProp("status").strip().replace("\xa0", "") \
        if mol.HasProp("status") else ""

    if status == "-":
        rows.append({
            "name": mol.GetProp("name"),
            "ID": mol.GetProp("ID"),
            "CAS": mol.GetProp("cas"),
            "status": status,
            "SMILES": Chem.MolToSmiles(Chem.RemoveHs(mol))
        })

df = pd.DataFrame(rows)

out = Path("02_raw_data/approved_drugs/e-Drug3D_status_unspecified.csv")
df.to_csv(out, index=False)

print("Status-unspecified molecules:", len(df))
print("Unique drug names:", df["name"].nunique())
print("Unique IDs:", df["ID"].nunique())
print("Saved:", out)
print("\nNote: These are NOT yet verified as currently approved drugs.")
