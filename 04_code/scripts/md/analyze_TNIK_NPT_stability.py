from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.home() / "Desktop" / "TNIK_IPF_Project"

DATA_DIR = ROOT / "06_docking/molecular_dynamics/analysis/npt_stability"
PLOT_DIR = DATA_DIR / "plots"

PLOT_DIR.mkdir(parents=True, exist_ok=True)

systems = [
    "ponatinib_N32_H",
    "ponatinib_N37_H",
    "tepotinib_pH74_plus1",
]

columns = [
    "Time_ps",
    "Potential_kJ_mol",
    "Total_Energy_kJ_mol",
    "Temperature_K",
    "Pressure_bar",
    "Volume_nm3",
    "Density_kg_m3",
]

all_blocks = []
all_data = {}

for system in systems:
    file = DATA_DIR / f"{system}_npt_energy.xvg"

    if not file.exists():
        raise FileNotFoundError(f"Missing file: {file}")

    rows = []

    with open(file) as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith(("#", "@")):
                continue

            values = [float(x) for x in line.split()]

            if len(values) != 7:
                raise ValueError(
                    f"Expected 7 columns in {file}, found {len(values)}"
                )

            rows.append(values)

    df = pd.DataFrame(rows, columns=columns)

    if df.empty:
        raise ValueError(f"No numerical data found in {file}")

    all_data[system] = df

    print(f"\n{'=' * 60}")
    print(system)
    print(f"Frames: {len(df)}")
    print(f"Time range: {df.Time_ps.min():.1f}–{df.Time_ps.max():.1f} ps")

    # Calculate averages in consecutive 200 ps blocks
    for start in range(0, 1000, 200):
        end = start + 200

        block = df[
            (df.Time_ps >= start) &
            (df.Time_ps < end if end < 1000 else df.Time_ps <= end)
        ]

        if block.empty:
            continue

        record = {
            "System": system,
            "Window_ps": f"{start}-{end}",
        }

        for col in columns[1:]:
            record[f"{col}_mean"] = block[col].mean()
            record[f"{col}_SD"] = block[col].std()

        all_blocks.append(record)

    # Plot each property separately
    properties = columns[1:]

    fig, axes = plt.subplots(3, 2, figsize=(13, 11))
    axes = axes.flatten()

    for ax, prop in zip(axes, properties):
        ax.plot(df.Time_ps, df[prop], linewidth=0.8)
        ax.set_title(prop.replace("_", " "))
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(prop)
        ax.grid(alpha=0.3)

    fig.suptitle(f"NPT Stability Analysis: {system}")
    fig.tight_layout()

    fig.savefig(
        PLOT_DIR / f"{system}_npt_stability.png",
        dpi=300
    )

    plt.close(fig)

# Save block-average results
block_df = pd.DataFrame(all_blocks)
block_df.to_csv(DATA_DIR / "NPT_200ps_block_averages.csv", index=False)

# Save a combined summary
summary = []

for system, df in all_data.items():
    for col in columns[1:]:
        summary.append({
            "System": system,
            "Property": col,
            "Mean": df[col].mean(),
            "SD": df[col].std(),
            "Minimum": df[col].min(),
            "Maximum": df[col].max(),
            "First_200ps_mean": df.loc[df.Time_ps < 200, col].mean(),
            "Last_200ps_mean": df.loc[df.Time_ps >= 800, col].mean(),
        })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(DATA_DIR / "NPT_overall_summary.csv", index=False)

print("\n" + "=" * 60)
print("NPT ANALYSIS COMPLETED")
print("=" * 60)
print(f"Block averages: {DATA_DIR / 'NPT_200ps_block_averages.csv'}")
print(f"Overall summary: {DATA_DIR / 'NPT_overall_summary.csv'}")
print(f"Plots: {PLOT_DIR}")
