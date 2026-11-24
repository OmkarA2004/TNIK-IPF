import pandas as pd
import matplotlib.pyplot as plt

# Load Y-randomization results
input_file = "../../05_models/metrics/TNIK_Y_randomization_results.csv"

output_file = "../../09_results/figures/TNIK_Y_randomization_plot.png"

df = pd.read_csv(input_file)

# Original model performance
original_r2 = 0.5439

# Plot randomized R2 values
plt.figure(figsize=(9, 6))

plt.scatter(
    df["Randomization"],
    df["Mean_R2"],
    s=70,
    label="Y-randomized models"
)

plt.axhline(
    y=original_r2,
    linestyle="--",
    label="Original model R² = 0.5439"
)

plt.axhline(
    y=0,
    linestyle=":",
    label="R² = 0"
)

plt.xlabel("Randomization Experiment")
plt.ylabel("Mean R²")
plt.title("Y-Randomization Test - TNIK Random Forest")

plt.xticks(df["Randomization"])
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()

plt.savefig(output_file, dpi=300)
plt.show()

print("Graph saved:", output_file)