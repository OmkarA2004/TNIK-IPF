import pandas as pd

df = pd.read_csv(
    "02_raw_data/approved_drugs/e-Drug3D_matched_registration.csv"
)

df["NDA_numbers"] = df["NDA_numbers"].fillna("").astype(str)
df["status_registration"] = (
    df["status_registration"].fillna("").astype(str).str.strip().str.upper()
)

df["historical_NDA_approval_evidence"] = df["NDA_numbers"].str.contains(
    r"APPROVAL=\s*AP\b", case=False, regex=True
)

df["discontinued"] = df["status_registration"].eq("DISCONTINUED")

print("Total matched structures:", len(df))
print("\nApproval evidence and status:")
print(pd.crosstab(
    df["historical_NDA_approval_evidence"],
    df["discontinued"],
    margins=True
).to_string())

print("\nYear zero:", (df["year"].astype(str).str.strip() == "0").sum())

print("\nSaved classification audit.")
df.to_csv(
    "02_raw_data/approved_drugs/e-Drug3D_approval_audit.csv",
    index=False
)
