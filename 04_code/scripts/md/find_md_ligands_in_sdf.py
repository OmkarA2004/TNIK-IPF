from rdkit import Chem
from pathlib import Path
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.warning")

source = Path("02_raw_data/approved_drugs/e-Drug3D_2197.sdf")
output = Path("06_docking/molecular_dynamics/ligands")
output.mkdir(parents=True, exist_ok=True)

targets = ["ponatinib", "tepotinib"]
found = {name: False for name in targets}

supplier = Chem.SDMolSupplier(str(source), removeHs=False)

for index, mol in enumerate(supplier):
    if mol is None:
        continue

    if not mol.HasProp("name"):
        continue

    name = mol.GetProp("name").strip().lower()

    for target in targets:
        if name == target:
            found[target] = True

            print(f"\nFOUND: {target.upper()}")
            print("SDF record index:", index)
            print("Drug name:", mol.GetProp("name"))
            print("SDF ID:", mol.GetProp("ID"))
            print("CAS:", mol.GetProp("cas"))
            print("Atoms:", mol.GetNumAtoms())
            print("3D conformer:", mol.GetNumConformers() > 0)
            print("SMILES:", Chem.MolToSmiles(mol))

            destination = output / f"{target}_original.sdf"
            writer = Chem.SDWriter(str(destination))
            writer.write(mol)
            writer.close()

            print("Saved:", destination)

print("\n===== SEARCH SUMMARY =====")
for target in targets:
    print(f"{target}: {'FOUND' if found[target] else 'NOT FOUND'}")
