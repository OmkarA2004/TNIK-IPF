from pathlib import Path
import pandas as pd

from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

source = Path("08_analogues/filtered/TNIK_BRICS_applicability_domain.csv")
outdir = Path("08_analogues/admet")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(source)

selected = df[
    (df["applicability_domain"] == "Inside") |
    (df["candidate_id"] == "TNIK_BRICS_0001")
].drop_duplicates("candidate_id")

def make_catalog(catalog_type):
    params = FilterCatalogParams()
    params.AddCatalog(catalog_type)
    return FilterCatalog(params)

pains = make_catalog(FilterCatalogParams.FilterCatalogs.PAINS)
brenk = make_catalog(FilterCatalogParams.FilterCatalogs.BRENK)

rows = []
smi_lines = []

for _, row in selected.iterrows():
    mol = Chem.MolFromSmiles(str(row["canonical_smiles"]))
    if mol is None:
        continue

    pains_alerts = sorted({
        x.GetDescription() for x in pains.GetMatches(mol)
    })
    brenk_alerts = sorted({
        x.GetDescription() for x in brenk.GetMatches(mol)
    })

    smi = Chem.MolToSmiles(mol)

    rows.append({
        "candidate_id": row["candidate_id"],
        "applicability_domain": row["applicability_domain"],
        "predicted_pIC50_RF": row["predicted_pIC50_RF"],
        "max_training_similarity": row["max_training_similarity"],
        "PAINS_alerts": "; ".join(pains_alerts) if pains_alerts else "None",
        "Brenk_alerts": "; ".join(brenk_alerts) if brenk_alerts else "None",
        "canonical_smiles": smi
    })

    smi_lines.append(f"{smi} {row['candidate_id']}")

alert_file = outdir / "TNIK_five_candidate_separate_alerts.csv"
smi_file = outdir / "TNIK_five_candidate_ADMET_input.smi"

pd.DataFrame(rows).to_csv(alert_file, index=False)
smi_file.write_text("\n".join(smi_lines) + "\n")

print("\n===== SEPARATE STRUCTURAL ALERTS =====")
print(pd.DataFrame(rows).to_string(index=False))

print("\nSMILES submission file:", smi_file)
print("Alert report:", alert_file)
