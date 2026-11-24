import pandas as pd

df = pd.read_csv(
    "02_raw_data/approved_drugs/e-Drug3D_approval_audit.csv"
)

library = df[
    (df["historical_NDA_approval_evidence"] == True) &
    (df["discontinued"] == False)
].copy()

library.to_csv(
    "02_raw_data/approved_drugs/TNIK_provisional_screening_library.csv",
    index=False
)

print("Provisional screening compounds:", len(library))
print("Unique e-Drug3D IDs:", library["ID"].nunique())
print("Unique drug names:", library["name_sdf"].nunique()
      if "name_sdf" in library.columns else "Check column names")

print("Saved provisional screening library.")
