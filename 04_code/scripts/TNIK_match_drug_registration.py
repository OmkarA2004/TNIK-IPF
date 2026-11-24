import pandas as pd
from pathlib import Path

sdf_data = pd.read_csv(
    "02_raw_data/approved_drugs/e-Drug3D_metadata_audit.csv",
    dtype=str
)

registration = pd.read_csv(
    "02_raw_data/approved_drugs/e-Drug3D_registration_parsed.csv",
    dtype=str
)

merged = sdf_data.merge(
    registration,
    on="ID",
    how="left",
    suffixes=("_sdf", "_registration"),
    indicator=True
)

print("SDF molecules:", len(sdf_data))
print("Matched registration records:", (merged["_merge"] == "both").sum())
print("Unmatched structures:", (merged["_merge"] == "left_only").sum())

print("\nRegistration status:")
print(merged["status_registration"].value_counts(dropna=False).to_string())

print("\nDiscontinued records:",
      (merged["status_registration"].str.upper() == "DISCONTINUED").sum())

print("\nRecords with NDA approval evidence:",
      merged["NDA_numbers"].fillna("").str.contains(
          r"APPROVAL=\s*AP\b", case=False
      ).sum())

out = Path(
    "02_raw_data/approved_drugs/e-Drug3D_matched_registration.csv"
)
merged.to_csv(out, index=False)

print("\nSaved:", out)
