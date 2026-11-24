from rdkit import Chem
from pathlib import Path
from collections import Counter

file = Path("02_raw_data/approved_drugs/e-Drug3D_2197.sdf")

supplier = Chem.SDMolSupplier(str(file), removeHs=False)

total = 0
valid = 0
invalid = 0
names = []
properties = Counter()
missing_smiles = 0

for mol in supplier:
    total += 1

    if mol is None:
        invalid += 1
        continue

    valid += 1
    names.append(mol.GetProp("_Name") if mol.HasProp("_Name") else "Unnamed")

    for prop in mol.GetPropNames():
        properties[prop] += 1

    if not Chem.MolToSmiles(Chem.RemoveHs(mol)):
        missing_smiles += 1

print("\nDRUG LIBRARY AUDIT")
print("------------------")
print("Total records:", total)
print("Valid molecules:", valid)
print("Invalid molecules:", invalid)
print("Unique names:", len(set(names)))
print("Missing SMILES:", missing_smiles)

print("\nAvailable metadata fields:")
for prop, count in properties.most_common():
    print(f"{prop}: {count}/{valid}")

print("\nFirst 10 molecule names:")
for name in names[:10]:
    print(name)
