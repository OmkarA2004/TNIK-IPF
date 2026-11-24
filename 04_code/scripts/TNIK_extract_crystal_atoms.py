from pathlib import Path
import pandas as pd

pdb = Path("02_raw_data/pdb_original/JWK_ligand/JWK_original.pdb")
rows = []

for line in pdb.read_text().splitlines():
    if line.startswith("HETATM"):
        element = line[76:78].strip()
        if element.upper() != "H":
            rows.append({
                "atom_name": line[12:16].strip(),
                "element": element,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })

df = pd.DataFrame(rows)
out = Path("06_docking/redocking/JWK_crystal_atoms.csv")
df.to_csv(out, index=False)

print("Crystal heavy atoms:", len(df))
print("Unique atom names:", df["atom_name"].nunique())
print("Saved:", out)
