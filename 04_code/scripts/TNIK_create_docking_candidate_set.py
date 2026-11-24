import pandas as pd

df = pd.read_csv(
    "02_raw_data/approved_drugs/TNIK_drug_property_audit.csv"
)

conditions = (
    (df["MW"] <= 800) &
    (df["HBD"] <= 10) &
    (df["HBA"] <= 15) &
    (df["Rotatable_Bonds"] <= 15)
)

candidates = df[conditions].copy()
review = df[~conditions].copy()

candidates.to_csv(
    "02_raw_data/approved_drugs/TNIK_small_molecule_candidates.csv",
    index=False
)

review.to_csv(
    "02_raw_data/approved_drugs/TNIK_large_or_complex_review.csv",
    index=False
)

print("Total compounds:", len(df))
print("Small-molecule docking candidates:", len(candidates))
print("Flagged for separate review:", len(review))
print("Total accounted for:", len(candidates) + len(review))
print("\nCandidate molecular weight range:")
print(candidates["MW"].describe().round(2).to_string())
