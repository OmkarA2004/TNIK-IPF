from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem

base = Path("06_docking/redocking")

# Read crystal and redocked coordinates
crystal = pd.read_csv(base / "JWK_crystal_atoms.csv")
docked = pd.read_csv(base / "JWK_redocked_model1_atoms.csv")

# Read the ligand SDF: its first 22 atoms are heavy atoms
mol = Chem.SDMolSupplier(
    "06_docking/ligand_prep/JWK_ideal.sdf",
    removeHs=False
)[0]

if mol is None:
    raise ValueError("Could not read ligand SDF")

# Read atom-index mapping from Model 1 remarks
pdbqt_text = Path(
    "06_docking/redocking/JWK_redocked.pdbqt"
).read_text()

mapping = {}
inside = False

for line in pdbqt_text.splitlines():
    if line.startswith("MODEL"):
        inside = line.split()[1] == "1"
    elif line.startswith("ENDMDL"):
        break
    elif inside and line.startswith("REMARK INDEX MAP"):
        values = list(map(int, line.split()[3:]))
        for i in range(0, len(values), 2):
            rdkit_index = values[i] - 1
            pdbqt_serial = values[i + 1]
            mapping[rdkit_index] = pdbqt_serial

# Match SDF heavy atoms to PDBQT heavy atoms
docked_by_serial = docked.set_index("serial")
crystal_xyz = crystal[["x", "y", "z"]].to_numpy()
docked_xyz = []

for i in range(22):
    if i not in mapping:
        raise ValueError(f"Missing atom mapping for SDF atom {i + 1}")
    serial = mapping[i]
    if serial not in docked_by_serial.index:
        raise ValueError(f"Missing PDBQT atom serial {serial}")
    docked_xyz.append(
        docked_by_serial.loc[serial, ["x", "y", "z"]].to_numpy(dtype=float)
    )

docked_xyz = np.array(docked_xyz)

if len(crystal_xyz) != 22 or len(docked_xyz) != 22:
    raise ValueError("Expected exactly 22 matched heavy atoms")

# Check element order
sdf_elements = [a.GetSymbol() for a in list(mol.GetAtoms())[:22]]
crystal_elements = crystal["element"].str.capitalize().tolist()

if sdf_elements != crystal_elements:
    raise ValueError("SDF and crystal atom order/elements do not match")

# Direct RMSD: no coordinate alignment
differences = docked_xyz - crystal_xyz
rmsd = np.sqrt(np.mean(np.sum(differences**2, axis=1)))

print("\nRedocking RMSD calculation")
print("--------------------------")
print("Matched heavy atoms:", len(docked_xyz))
print("Coordinate alignment: None")
print(f"Direct heavy-atom RMSD: {rmsd:.3f} Å")

if rmsd <= 2.0:
    print("Interpretation: Meets the commonly used 2 Å redocking criterion.")
else:
    print("Interpretation: Does not meet the commonly used 2 Å criterion.")

result = pd.DataFrame({
    "atom_name": crystal["atom_name"],
    "crystal_x": crystal_xyz[:, 0],
    "crystal_y": crystal_xyz[:, 1],
    "crystal_z": crystal_xyz[:, 2],
    "docked_x": docked_xyz[:, 0],
    "docked_y": docked_xyz[:, 1],
    "docked_z": docked_xyz[:, 2],
    "displacement_A": np.linalg.norm(differences, axis=1)
})

result.to_csv(base / "JWK_redocking_atom_comparison.csv", index=False)
print("Atom comparison saved.")
