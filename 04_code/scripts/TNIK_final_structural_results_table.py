import pandas as pd

source = "09_results/tables/TNIK_four_drug_master_structural_comparison.csv"
output = "09_results/tables/TNIK_final_structural_results_table.csv"

df = pd.read_csv(source)

columns = {
    "Drug": "Drug",
    "Docking_score_kcal_mol": "Docking score (kcal/mol)",
    "Nearby_residues_within_4A": "Nearby residues (within 4 Å)",
    "JWK_shared_pocket_residues": "Shared JWK pocket residues",
    "Geometry_Hbond_candidates": "Potential H-bond geometry candidates",
    "Carbon_atom_pairs_within_4A": "Potential carbon contacts (pairs)",
    "Closest_heavy_atom_contact_A": "Closest heavy-atom contact (Å)"
}

result = df[list(columns)].rename(columns=columns)

result.to_csv(output, index=False)

print("\nFINAL STRUCTURAL RESULTS TABLE\n")
print(result.to_string(index=False))
print("\nSaved to:", output)
