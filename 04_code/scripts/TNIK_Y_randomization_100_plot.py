

import os
import pandas as pd
import matplotlib.pyplot as plt

project_dir = os.path.expanduser("~/Desktop/TNIK_IPF_Project")

input_file = os.path.join(
    project_dir,
    "05_models/metrics/TNIK_Y_randomization_100_results.csv"
)

output_file = os.path.join(
    project_dir,
    "09_results/figures/TNIK_Y_randomization_100_labelled_plot.png"
)

df = pd.read_csv(input_file)

original_r2 = 0.5439

plt.figure(figsize=(16, 9))

plt.scatter(
    df["Randomization"],
    df["Mean_R2"],
    s=45,
    alpha=0.8,
    label="Randomized model R²"
)

# Label every point
for _, row in df.iterrows():
    plt.annotate(
        str(int(row["Randomization"])),
        (row["Randomization"], row["Mean_R2"]),
        xytext=(3, 5),
        textcoords="offset points",
        fontsize=7
    )

plt.axhline(
    original_r2,
    linestyle="--",
    linewidth=2,
    label="Original model R² = 0.5439"
)

plt.axhline(
    0,
    linestyle=":",
    linewidth=1.5,
    label="R² = 0"
)

plt.xlabel("Randomization Number", fontsize=12)
plt.ylabel("Mean R²", fontsize=12)
plt.title("Y-Randomization Test: 100 Labelled Permutations", fontsize=14)

plt.xticks(range(0, 101, 5))
plt.grid(alpha=0.25)
plt.legend()
plt.tight_layout()

plt.savefig(output_file, dpi=300)
plt.show()

print("Labelled plot saved to:")
print(output_file)
