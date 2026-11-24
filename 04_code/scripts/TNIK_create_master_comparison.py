import pandas as pd

base = "09_results/tables/"

# Read existing integrated comparison
df = pd.read_csv(base + "TNIK_four_drug_final_comparison.csv")

# Crystal ligand JWK overlap measurements
overlap = pd.DataFrame([
    {
        "Drug": "Ponatinib",
        "JWK_closest_heavy_atom_A": 0.402,
        "JWK_atom_pairs_within_2A": 60,
        "JWK_atom_pairs_within_3A": 130,
        "JWK_shared_pocket_residues": 10
    },
    {
        "Drug": "Tepotinib",
        "JWK_closest_heavy_atom_A": 0.299,
        "JWK_atom_pairs_within_2A": 72,
        "JWK_atom_pairs_within_3A": 145,
        "JWK_shared_pocket_residues": 14
    },
    {
        "Drug": "Entrectinib",
        "JWK_closest_heavy_atom_A": 0.403,
        "JWK_atom_pairs_within_2A": 63,
        "JWK_atom_pairs_within_3A": 134,
        "JWK_shared_pocket_residues": 12
    },
    {
        "Drug": "Nilotinib",
        "JWK_closest_heavy_atom_A": 0.415,
        "JWK_atom_pairs_within_2A": 60,
        "JWK_atom_pairs_within_3A": 131,
        "JWK_shared_pocket_residues": 12
    }
])

# Merge without modifying the original table
master = df.merge(overlap, on="Drug", how="left", validate="one_to_one")

output = base + "TNIK_four_drug_master_structural_comparison.csv"
master.to_csv(output, index=False)

print("Master comparison saved:", output)
print("Table dimensions:", master.shape)
print("\nMaster comparison:")
print(master.to_string(index=False))
