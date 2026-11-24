import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

analysis = Path("06_docking/molecular_dynamics/analysis")
results = Path("09_results/molecular_dynamics")
results.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(analysis / "TNIK_apo_RMSF_original_numbering.csv")

plt.figure(figsize=(11, 5))
plt.plot(
    df["Original_PDB_residue"],
    df["RMSF_nm"],
    linewidth=1.3
)

for residue, label in [(181, "TPO181"), (187, "TPO187")]:
    plt.axvline(residue, linestyle="--", linewidth=1)
    plt.text(
        residue, df["RMSF_nm"].max() * 0.95,
        label, rotation=90, va="top", ha="right"
    )

plt.xlabel("Original 6RA7 residue number")
plt.ylabel("C-alpha RMSF (nm)")
plt.title("Apo-TNIK residue-wise RMSF (600 ps)")
plt.tight_layout()

output = results / "TNIK_apo_RMSF_original_numbering.png"
plt.savefig(output, dpi=300)
plt.close()

print("Corrected RMSF figure saved:", output)
