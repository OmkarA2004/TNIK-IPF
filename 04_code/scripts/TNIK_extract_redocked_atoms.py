from pathlib import Path
import pandas as pd

pdbqt = Path("06_docking/redocking/JWK_redocked.pdbqt")
rows = []
inside_model_1 = False

for line in pdbqt.read_text().splitlines():
    if line.startswith("MODEL"):
        inside_model_1 = line.split()[1] == "1"
        continue

    if line.startswith("ENDMDL"):
        if inside_model_1:
            break
        continue

    if inside_model_1 and line.startswith(("ATOM", "HETATM")):
        element = line[76:78].strip()
        rows.append({
            "serial": int(line[6:11]),
            "element": element,
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
        })

df = pd.DataFrame(rows)
out = Path("06_docking/redocking/JWK_redocked_model1_atoms.csv")
df.to_csv(out, index=False)

print("Model 1 atoms:", len(df))
print("Heavy atoms:", sum(df["element"].str.upper() != "H"))
print("Saved:", out)
