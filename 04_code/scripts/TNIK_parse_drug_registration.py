import pandas as pd
from pathlib import Path

source = Path(
    "02_raw_data/approved_drugs/e-Drug3D_registration_data.txt"
)

rows = []

with open(source, encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        fields = [x.strip() for x in line.split(";")]

        if len(fields) < 8:
            continue

        rows.append({
            "ID": fields[0],
            "name": fields[1],
            "CAS": fields[2],
            "year": fields[3],
            "status": fields[4],
            "indications": fields[6],
            "NDA_numbers": fields[7]
        })

df = pd.DataFrame(rows)

out = Path(
    "02_raw_data/approved_drugs/e-Drug3D_registration_parsed.csv"
)
df.to_csv(out, index=False)

print("Registration records:", len(df))
print("Unique IDs:", df["ID"].nunique())

print("\nStatus distribution:")
print(df["status"].value_counts(dropna=False).to_string())

print("\nApproval-year summary:")
print(df["year"].value_counts().head(10).to_string())

print("\nSaved:", out)
