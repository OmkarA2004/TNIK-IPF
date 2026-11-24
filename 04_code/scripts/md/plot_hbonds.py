import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

source = Path("09_results/molecular_dynamics/TNIK_apo_hbond_timeseries.csv")
output = Path("09_results/molecular_dynamics")
df = pd.read_csv(source)

plt.figure(figsize=(10, 5))
plt.plot(df["Time_ns"], df["Hydrogen_bonds"], linewidth=0.8)
plt.axhline(df["Hydrogen_bonds"].mean(), linestyle="--",
            label=f'Mean = {df["Hydrogen_bonds"].mean():.2f}')
plt.xlabel("Time (ns)")
plt.ylabel("Number of intraprotein hydrogen bonds")
plt.title("Apo-TNIK Intraprotein Hydrogen Bonds")
plt.legend()
plt.tight_layout()
plt.savefig(output / "TNIK_apo_hbond_timeseries.png",
            dpi=300, bbox_inches="tight")
plt.close()

print("Hydrogen-bond plot saved successfully.")
print(output / "TNIK_apo_hbond_timeseries.png")
