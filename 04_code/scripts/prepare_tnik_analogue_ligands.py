from pathlib import Path
import subprocess
import shutil
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

source = Path(
    "08_analogues/filtered/TNIK_BRICS_applicability_domain.csv"
)

outdir = Path("06_docking/ligand_prep/analogue_candidates")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(source)

selected = df[
    (df["applicability_domain"] == "Inside") |
    (df["candidate_id"] == "TNIK_BRICS_0001")
].drop_duplicates("candidate_id")

prepare_tool = shutil.which("mk_prepare_ligand.py")

if prepare_tool is None:
    raise SystemExit(
        "ERROR: mk_prepare_ligand.py not found. "
        "Activate the tnik_ipf environment."
    )

results = []

for _, row in selected.iterrows():
    candidate_id = row["candidate_id"]
    smiles = row["canonical_smiles"]

    print(f"\nPreparing {candidate_id}...")

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("Invalid SMILES. Skipping.")
        continue

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    status = AllChem.EmbedMolecule(mol, params)

    if status != 0:
        print("3D embedding failed. Skipping.")
        continue

    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        force_field = "MMFF"
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
        force_field = "UFF"

    mol.SetProp("_Name", candidate_id)

    sdf_file = outdir / f"{candidate_id}.sdf"
    pdbqt_file = outdir / f"{candidate_id}.pdbqt"

    writer = Chem.SDWriter(str(sdf_file))
    writer.write(mol)
    writer.close()

    command = [
        prepare_tool,
        "-i", str(sdf_file),
        "-o", str(pdbqt_file)
    ]

    run = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    if run.returncode == 0 and pdbqt_file.exists():
        print(f"PDBQT created: {pdbqt_file}")
        results.append({
            "candidate_id": candidate_id,
            "force_field": force_field,
            "status": "Prepared",
            "SDF": str(sdf_file),
            "PDBQT": str(pdbqt_file)
        })
    else:
        print("Meeko preparation failed:")
        print(run.stderr[-1500:])
        results.append({
            "candidate_id": candidate_id,
            "force_field": force_field,
            "status": "Failed",
            "SDF": str(sdf_file),
            "PDBQT": ""
        })

report = pd.DataFrame(results)
report.to_csv(
    outdir / "TNIK_analogue_ligand_preparation_report.csv",
    index=False
)

print("\n===== LIGAND PREPARATION SUMMARY =====")
print(report.to_string(index=False))
print("\nOutput directory:", outdir)
