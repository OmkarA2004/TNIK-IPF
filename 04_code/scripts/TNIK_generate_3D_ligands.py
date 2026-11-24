import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

source = "02_raw_data/approved_drugs/TNIK_small_molecule_candidates.csv"
output_dir = Path("06_docking/ligand_prep/batch_sdf")
output_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(source)
success = []
failed = []

for i, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["SMILES"])

    if mol is None:
        failed.append((row["ID"], "Invalid SMILES"))
        continue

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    result = AllChem.EmbedMolecule(mol, params)

    if result != 0:
        failed.append((row["ID"], "3D generation failed"))
        continue

    if AllChem.MMFFHasAllMoleculeParams(mol):
        status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        method = "MMFF"
    else:
        try:
            status = AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            method = "UFF"
        except Exception:
            failed.append((row["ID"], "No suitable force field"))
            continue

    mol.SetProp("_Name", str(row["ID"]))
    mol.SetProp("eDrug3D_ID", str(row["ID"]))
    mol.SetProp("SMILES", str(row["SMILES"]))
    mol.SetProp("Optimization_Method", method)
    mol.SetProp("Optimization_Status", str(status))

    writer = Chem.SDWriter(
        str(output_dir / f"drug_{row['ID']}.sdf")
    )
    writer.write(mol)
    writer.close()

    success.append((row["ID"], method, status))

    if len(success) % 100 == 0:
        print("Successfully prepared:", len(success))

pd.DataFrame(
    success,
    columns=["ID", "Method", "Optimization_Status"]
).to_csv(
    "06_docking/ligand_prep/TNIK_3D_ligand_success.csv",
    index=False
)

pd.DataFrame(
    failed,
    columns=["ID", "Failure_Reason"]
).to_csv(
    "06_docking/ligand_prep/TNIK_3D_ligand_failures.csv",
    index=False
)

print("\n3D preparation finished.")
print("Successful:", len(success))
print("Failed:", len(failed))
print("Output folder:", output_dir)
