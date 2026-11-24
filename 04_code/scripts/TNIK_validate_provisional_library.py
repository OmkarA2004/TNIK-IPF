import pandas as pd
from rdkit import Chem

library = pd.read_csv(
    "02_raw_data/approved_drugs/TNIK_provisional_screening_library.csv",
    dtype=str
)

selected_ids = set(library["ID"].dropna())

supplier = Chem.SDMolSupplier(
    "02_raw_data/approved_drugs/e-Drug3D_2197.sdf",
    removeHs=False
)

valid = []
invalid = []
matched_ids = set()

for mol in supplier:
    if mol is None:
        continue

    drug_id = mol.GetProp("ID") if mol.HasProp("ID") else ""

    if drug_id not in selected_ids:
        continue

    matched_ids.add(drug_id)
    smiles = Chem.MolToSmiles(mol)

    if smiles:
        valid.append({
            "ID": drug_id,
            "SMILES": smiles
        })
    else:
        invalid.append(drug_id)

print("Selected IDs:", len(selected_ids))
print("Matched SDF IDs:", len(matched_ids))
print("Valid structures:", len(valid))
print("Invalid structures:", len(invalid))
print("IDs not found:", len(selected_ids - matched_ids))

pd.DataFrame(valid).to_csv(
    "02_raw_data/approved_drugs/TNIK_provisional_valid_structures.csv",
    index=False
)

print("Saved validated structures.")
