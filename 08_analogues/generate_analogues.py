from pathlib import Path
import requests
import pandas as pd

from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "08_analogues" / "generated"
OUT.mkdir(parents=True, exist_ok=True)

PARENTS = ["Ponatinib", "Tepotinib", "Entrectinib", "Nilotinib"]

def get_pubchem_smiles(name):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"name/{name}/property/SMILES/JSON"
    )
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    properties = response.json()["PropertyTable"]["Properties"]
    if not properties or not properties[0].get("SMILES"):
        raise ValueError(f"No SMILES returned for {name}")

    return properties[0]["SMILES"]

records = []
all_candidates = {}

for name in PARENTS:
    print(f"\nProcessing {name}...")

    try:
        smiles = get_pubchem_smiles(name)
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            raise ValueError("RDKit could not parse the retrieved SMILES")

        canonical = Chem.MolToSmiles(mol)

        records.append({
            "compound": name,
            "parent_compound": name,
            "smiles": canonical,
            "role": "parent"
        })

        fragments = list(BRICS.BRICSDecompose(mol))
        fragment_mols = [
            Chem.MolFromSmiles(fragment) for fragment in fragments
        ]
        fragment_mols = [m for m in fragment_mols if m is not None]

        generated = 0

        if len(fragment_mols) >= 2:
            for candidate in BRICS.BRICSBuild(
                fragment_mols,
                maxDepth=2,
                onlyCompleteMols=True
            ):
                try:
                    Chem.SanitizeMol(candidate)
                    candidate_smiles = Chem.MolToSmiles(candidate)

                    if candidate_smiles == canonical:
                        continue

                    all_candidates.setdefault(candidate_smiles, name)
                    generated += 1

                    if generated >= 500:
                        break
                except Exception:
                    continue

        print("Structure retrieved successfully.")
        print("BRICS fragments:", len(fragments))
        print("Candidate structures collected:", generated)

    except Exception as error:
        print(f"FAILED: {name}: {error}")

if not records:
    raise SystemExit(
        "\nNo parent structures were retrieved. No analogue generation performed."
    )

for smiles, parent in all_candidates.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue

    records.append({
        "compound": f"{parent}_analogue",
        "parent_compound": parent,
        "smiles": smiles,
        "role": "generated_analogue",
        "molecular_weight": Descriptors.MolWt(mol),
        "logP": Descriptors.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol)
    })

df = pd.DataFrame(records).drop_duplicates(subset="smiles")

output_file = OUT / "TNIK_generated_analogues.csv"
df.to_csv(output_file, index=False)

print("\n========== ANALOGUE GENERATION SUMMARY ==========")
print("Parent compounds retrieved:", (df["role"] == "parent").sum())
print("Generated analogues:", (df["role"] == "generated_analogue").sum())
print("Total unique structures:", len(df))
print("Saved to:", output_file)
