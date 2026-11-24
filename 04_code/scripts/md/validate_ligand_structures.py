from rdkit import Chem
from pathlib import Path

ligands = {
    "ponatinib": 11,
    "tepotinib": 22,
}

base = Path("06_docking/molecular_dynamics/ligands")
docking = Path("06_docking/reproducibility")

for name, seed in ligands.items():
    print(f"\n===== {name.upper()} =====")

    sdf_path = base / f"{name}_original.sdf"
    pdbqt_path = docking / f"{name}_seed_{seed}_out.pdbqt"

    mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not read {sdf_path}")

    heavy = Chem.RemoveHs(mol)
    print("SDF heavy atoms:", heavy.GetNumAtoms())
    print("SDF bonds:", heavy.GetNumBonds())
    print("SDF canonical isomeric SMILES:",
          Chem.MolToSmiles(heavy, isomericSmiles=True))

    lines = pdbqt_path.read_text().splitlines()
    models = []
    current = []

    for line in lines:
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            models.append(current)
        elif line.startswith(("ATOM", "HETATM")):
            current.append(line)

    if not models:
        raise ValueError(f"No docking models found in {pdbqt_path}")

    atoms = models[0]
    print("Selected docking seed:", seed)
    print("Models in file:", len(models))
    print("Atoms in first pose:", len(atoms))
    print("PDBQT atom types:",
          sorted(set(line[77:79].strip() for line in atoms)))

    print("First-pose atom names:")
    print(" ".join(line[12:16].strip() for line in atoms))

