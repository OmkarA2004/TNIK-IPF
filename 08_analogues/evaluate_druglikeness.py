from pathlib import Path
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "08_analogues/filtered/TNIK_analogues_filtered.csv"
OUT = ROOT / "08_analogues/admet"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SOURCE)
results = []

for _, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        continue

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10
    ])

    results.append({
        "parent_compound": row["parent_compound"],
        "smiles": row["smiles"],
        "tanimoto_similarity": row["tanimoto_similarity"],
        "molecular_weight": mw,
        "logP": logp,
        "HBD": hbd,
        "HBA": hba,
        "TPSA": tpsa,
        "rotatable_bonds": rotatable,
        "lipinski_violations": violations,
        "lipinski_compliant": violations <= 1,
        "veber_compliant": rotatable <= 10 and tpsa <= 140
    })

result = pd.DataFrame(results)

if result.empty:
    raise SystemExit("No valid structures available.")

result["both_guidelines_compliant"] = (
    result["lipinski_compliant"] & result["veber_compliant"]
)

path = OUT / "TNIK_analogue_druglikeness.csv"
result.to_csv(path, index=False)

print("\n========== DRUG-LIKENESS SUMMARY ==========")
print("Analogues evaluated:", len(result))
print("Lipinski compliant:", result["lipinski_compliant"].sum())
print("Veber compliant:", result["veber_compliant"].sum())
print("Both guidelines compliant:", result["both_guidelines_compliant"].sum())
print("\nResults by parent:")
print(result.groupby("parent_compound")[
    ["lipinski_compliant", "veber_compliant",
     "both_guidelines_compliant"]
].sum().to_string())
print("\nSaved to:", path)
