import pandas as pd
from pathlib import Path

results = pd.DataFrame([{
    "PDB_ID": "6RA7",
    "Target": "Human TNIK",
    "Ligand": "JWK",
    "Vina_version": "1.2.7",
    "Scoring_function": "Vina",
    "Best_affinity_kcal_mol": -9.770,
    "Heavy_atoms_matched": 22,
    "Direct_RMSD_A": 0.782,
    "RMSD_threshold_A": 2.0,
    "Redocking_validation": "Successful",
    "Exhaustiveness": 32,
    "Random_seed": 42,
    "Grid_center": "10.794, 10.872, 6.540",
    "Grid_size_A": "22 x 22 x 22"
}])

out = Path("09_results/tables/TNIK_redocking_validation_summary.csv")
out.parent.mkdir(parents=True, exist_ok=True)
results.to_csv(out, index=False)

print("Redocking validation summary saved:")
print(out)
print(results.to_string(index=False))
