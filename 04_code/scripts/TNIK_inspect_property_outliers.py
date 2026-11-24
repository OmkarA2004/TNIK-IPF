import pandas as pd

df = pd.read_csv(
    "02_raw_data/approved_drugs/TNIK_drug_property_audit.csv"
)

outliers = df[
    (df["MW"] > 800) |
    (df["HBD"] > 10) |
    (df["HBA"] > 15) |
    (df["Rotatable_Bonds"] > 15) |
    (df["LogP"] > 6) |
    (df["LogP"] < -5)
].copy()

outliers = outliers.sort_values("MW", ascending=False)

outliers.to_csv(
    "02_raw_data/approved_drugs/TNIK_property_outliers.csv",
    index=False
)

print("Total compounds:", len(df))
print("Flagged compounds:", len(outliers))
print("\nTop 15 flagged compounds:")
print(
    outliers[
        ["ID", "MW", "HBD", "HBA", "Rotatable_Bonds", "LogP"]
    ].head(15).to_string(index=False)
)
print("\nSaved TNIK_property_outliers.csv")
