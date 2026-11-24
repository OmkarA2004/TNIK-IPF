import pandas as pd

outliers = pd.read_csv(
    "02_raw_data/approved_drugs/TNIK_property_outliers.csv",
    dtype={"ID": str}
)

library = pd.read_csv(
    "02_raw_data/approved_drugs/TNIK_provisional_screening_library.csv",
    dtype={"ID": str}
)

merged = outliers.merge(
    library[["ID", "name_sdf", "CAS_sdf"]],
    on="ID",
    how="left"
)

merged = merged.sort_values("MW", ascending=False)

merged.to_csv(
    "02_raw_data/approved_drugs/TNIK_property_outliers_named.csv",
    index=False
)

print("Flagged compounds:", len(merged))
print("\nTop 20 flagged compounds:")
print(
    merged[
        ["ID", "name_sdf", "CAS_sdf", "MW", "HBD", "HBA",
         "Rotatable_Bonds", "LogP"]
    ].head(20).to_string(index=False)
)

print("\nSaved TNIK_property_outliers_named.csv")
