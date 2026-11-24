from rdkit import Chem
from pathlib import Path

sdf = Path(
    "06_docking/molecular_dynamics/ligands/ponatinib_original.sdf"
)
pdbqt = Path(
    "06_docking/reproducibility/ponatinib_seed_11_out.pdbqt"
)

mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]

print("===== ORIGINAL SDF ATOMS =====")
conf = mol.GetConformer()

for atom in mol.GetAtoms():
    i = atom.GetIdx()
    p = conf.GetAtomPosition(i)
    print(
        f"{i:3d}  {atom.GetSymbol():2s}  "
        f"{p.x:10.4f} {p.y:10.4f} {p.z:10.4f}"
    )

print("\n===== FIRST DOCKING POSE =====")
in_first_model = False

for line in pdbqt.read_text().splitlines():
    if line.startswith("MODEL"):
        in_first_model = True
        continue

    if line.startswith("ENDMDL"):
        break

    if in_first_model and line.startswith(("ATOM", "HETATM")):
        name = line[12:16].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atom_type = line[77:79].strip()

        print(
            f"{name:5s} {atom_type:3s} "
            f"{x:10.4f} {y:10.4f} {z:10.4f}"
        )

print("\n===== DOCKING TORSION TREE =====")
for line in pdbqt.read_text().splitlines():
    if line.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH")):
        print(line)
