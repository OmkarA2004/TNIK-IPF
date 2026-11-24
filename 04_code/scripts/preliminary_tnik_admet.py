from pathlib import Path
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

source = Path("08_analogues/filtered/TNIK_BRICS_applicability_domain.csv")
outdir = Path("08_analogues/admet")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(source)

inside = df[df["applicability_domain"] == "Inside"].copy()
exploratory = df[df["candidate_id"] == "TNIK_BRICS_0001"].copy()

shortlist = pd.concat([inside, exploratory], ignore_index=True)
shortlist = shortlist.drop_duplicates(subset="candidate_id")

params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog(params)

results = []

for _, row in shortlist.iterrows():
    mol = Chem.MolFromSmiles(str(row["canonical_smiles"]))
    if mol is None:
        continue

    alerts = catalog.GetMatches(mol)
    alert_names = sorted({
        match.GetDescription() for match in alerts
    })

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    qed = QED.qed(mol)

    lipinski_violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10
    ])

    results.append({
        "candidate_id": row["candidate_id"],
        "applicability_domain": row["applicability_domain"],
        "predicted_pIC50_RF": row["predicted_pIC50_RF"],
        "max_training_similarity": row["max_training_similarity"],
        "molecular_weight": mw,
        "logP": logp,
        "TPSA": tpsa,
        "HBD": hbd,
        "HBA": hba,
        "rotatable_bonds": rotb,
        "fraction_CSP3": fsp3,
        "QED": qed,
        "Lipinski_violations": lipinski_violations,
        "Lipinski_pass": lipinski_violations == 0,
        "PAINS_BRENK_alert_count": len(alert_names),
        "structural_alerts": "; ".join(alert_names) if alert_names else "None detected",
        "canonical_smiles": Chem.MolToSmiles(mol)
    })

result = pd.DataFrame(results)

result["preliminary_triage"] = result.apply(
    lambda r: (
        "Review structural alerts"
        if r["PAINS_BRENK_alert_count"] > 0
        else "Review physicochemical properties"
        if not r["Lipinski_pass"]
        else "No alerts in this preliminary screen"
    ),
    axis=1
)

result = result.sort_values(
    ["applicability_domain", "predicted_pIC50_RF"],
    ascending=[True, False]
)

outfile = outdir / "TNIK_five_candidate_preliminary_ADMET.csv"
result.to_csv(outfile, index=False)

print("\n===== PRELIMINARY ADMET SCREEN =====")
print(result.to_string(index=False))
print("\nSaved:", outfile)
