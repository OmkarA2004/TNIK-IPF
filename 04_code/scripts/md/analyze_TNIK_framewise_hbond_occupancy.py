from pathlib import Path
import re
import pandas as pd

ndx = Path(
    "06_docking/molecular_dynamics/analysis/npt_stability/"
    "tepotinib_hbond_pairs_per_frame.ndx"
)

output_dir = Path(
    "06_docking/molecular_dynamics/analysis/npt_stability/hbond_pairs"
)
output_dir.mkdir(parents=True, exist_ok=True)

if not ndx.exists():
    raise FileNotFoundError(f"Missing file: {ndx}")

text = ndx.read_text()

sections = re.split(r"(?m)^\[\s*(.*?)\s*\]\s*$", text)

frame_interactions = {}

for i in range(1, len(sections), 2):
    section_name = sections[i]
    section_data = sections[i + 1]

    match = re.search(r"frame_(\d+)", section_name)
    if not match:
        continue

    frame = int(match.group(1))
    entries = []

    for line in section_data.splitlines():
        values = line.split()
        if len(values) == 3 and all(v.isdigit() for v in values):
            entries.append(tuple(map(int, values)))

    frame_interactions[frame] = entries

total_frames = len(frame_interactions)

if total_frames == 0:
    raise ValueError("No frame sections were found.")

# The same interaction is identified by its protein and ligand atom indices.
interactions = {
    "GLU17_OE2--Tepotinib_N3_H14": (4895, 4932, 261),
    "CYS96_N_H--Tepotinib_O2": (1515, 1516, 4909),
}

rows = []

for name, atom_tuple in interactions.items():
    occurrence_frames = [
        frame
        for frame, entries in frame_interactions.items()
        if atom_tuple in entries
    ]

    count = len(occurrence_frames)

    rows.append({
        "interaction": name,
        "frames_present": count,
        "total_frames": total_frames,
        "occupancy_percent": 100 * count / total_frames,
        "first_frame": min(occurrence_frames) if occurrence_frames else None,
        "last_frame": max(occurrence_frames) if occurrence_frames else None,
    })

df = pd.DataFrame(rows)

out_csv = output_dir / "tepotinib_hbond_occupancy.csv"
df.to_csv(out_csv, index=False)

print(f"Frames analysed: {total_frames}")
print("\nHydrogen-bond occupancy:")
print(df.to_string(index=False))
print(f"\nSaved: {out_csv}")
