from pathlib import Path
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "08_analogues/generated/TNIK_generated_analogues.csv"
OUT = ROOT / "08_analogues/filtered"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SOURCE)

parents = {}
for _, row in df[df["role"] == "parent"].iterrows():
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is not None:
        parents[row["compound"]] = mol

results = []

for _, row in df[df["role"] == "generated_analogue"].iterrows():
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        continue

    parent_name = row["parent_compound"]
    parent_mol = parents.get(parent_name)

    if parent_mol is None:
        continue

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    parent_fp = AllChem.GetMorganFingerprintAsBitVect(
        parent_mol, 2, nBits=2048
    )

    similarity = DataStructs.TanimotoSimilarity(fp, parent_fp)

    results.append({
        "parent_compound": parent_name,
        "smiles": Chem.MolToSmiles(mol),
        "tanimoto_similarity": similarity,
        "molecular_weight": row["molecular_weight"],
        "logP": row["logP"],
        "TPSA": row["TPSA"],
        "HBD": row["HBD"],
        "HBA": row["HBA"],
        "rotatable_bonds": row["rotatable_bonds"]
    })

result = pd.DataFrame(results)

if result.empty:
    raise SystemExit("No valid analogue structures available for filtering.")

result = result.drop_duplicates(subset="smiles")
filtered = result[result["tanimoto_similarity"] >= 0.40].copy()

result.to_csv(OUT / "TNIK_analogue_similarity_all.csv", index=False)
filtered.to_csv(OUT / "TNIK_analogues_filtered.csv", index=False)

print("\n========== STRUCTURAL FILTERING SUMMARY ==========")
print("Valid unique candidates:", len(result))
print("Similarity threshold: 0.40")
print("Candidates retained:", len(filtered))
print("Candidates excluded:", len(result) - len(filtered))
print("\nRetained by parent:")
print(filtered.groupby("parent_compound").size().to_string())
print("\nSaved in:", OUT)
