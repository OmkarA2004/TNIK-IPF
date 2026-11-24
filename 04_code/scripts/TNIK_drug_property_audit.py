import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

df = pd.read_csv(
    "02_raw_data/approved_drugs/TNIK_provisional_valid_structures.csv"
)

results = []

for _, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["SMILES"])

    if mol is None:
        continue

    results.append({
        "ID": row["ID"],
        "SMILES": row["SMILES"],
        "MW": Descriptors.MolWt(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "Rotatable_Bonds": Lipinski.NumRotatableBonds(mol),
        "LogP": Descriptors.MolLogP(mol)
    })

props = pd.DataFrame(results)

props.to_csv(
    "02_raw_data/approved_drugs/TNIK_drug_property_audit.csv",
    index=False
)

print("Compounds audited:", len(props))
print("\nProperty summary:")
print(props[["MW", "HBD", "HBA", "Rotatable_Bonds", "LogP"]].describe().round(2))

print("\nSaved TNIK_drug_property_audit.csv")
