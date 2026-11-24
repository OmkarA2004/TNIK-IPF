from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors
import json

base = Path("06_docking/molecular_dynamics/ligands")
report_dir = Path("06_docking/molecular_dynamics/logs")
report_dir.mkdir(parents=True, exist_ok=True)

report = {}

for name in ["ponatinib", "tepotinib"]:
    print(f"\n===== {name.upper()} =====")

    files = {
        "original": base / f"{name}_original.sdf",
        "docked": base / f"{name}_docked_first_pose_obabel.sdf"
    }

    report[name] = {}

    for label, path in files.items():
        mol = next(
            (m for m in Chem.SDMolSupplier(
                str(path), removeHs=False
            ) if m is not None),
            None
        )

        if mol is None:
            report[name][label] = {"error": "Molecule unreadable"}
            continue

        Chem.AssignStereochemistry(
            mol, cleanIt=True, force=True
        )

        chiral = Chem.FindMolChiralCenters(
            mol, includeUnassigned=True,
            includeCIP=True
        )

        nitrogens = []

        for atom in mol.GetAtoms():
            if atom.GetSymbol() != "N":
                continue

            neighbours = [
                f"{n.GetSymbol()}(index={n.GetIdx()})"
                for n in atom.GetNeighbors()
            ]

            nitrogens.append({
                "index": atom.GetIdx(),
                "formal_charge": atom.GetFormalCharge(),
                "aromatic": atom.GetIsAromatic(),
                "degree": atom.GetDegree(),
                "neighbors": neighbours
            })

        data = {
            "atoms": mol.GetNumAtoms(),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "bonds": mol.GetNumBonds(),
            "formal_charge": Chem.GetFormalCharge(mol),
            "molecular_weight": round(
                Descriptors.MolWt(mol), 3
            ),
            "canonical_smiles": Chem.MolToSmiles(
                Chem.RemoveHs(mol)
            ),
            "chiral_centers": [
                {"atom": i, "configuration": c}
                for i, c in chiral
            ],
            "nitrogen_environments": nitrogens
        }

        report[name][label] = data

        print(f"\n{label.upper()}")
        print("Atoms:", data["atoms"])
        print("Heavy atoms:", data["heavy_atoms"])
        print("Formal charge:", data["formal_charge"])
        print("Chiral centers:", data["chiral_centers"])
        print("Nitrogen atoms:", len(nitrogens))

        for n in nitrogens:
            print(
                f"  N{n['index']}: charge={n['formal_charge']}, "
                f"aromatic={n['aromatic']}, "
                f"neighbors={n['neighbors']}"
            )

out = report_dir / "TNIK_final_ligand_chemistry_check.json"
out.write_text(json.dumps(report, indent=2))

print("\n==========================================")
print("LIGAND CHEMISTRY INSPECTION COMPLETE")
print("Report:", out)
print("Original molecular files were not modified.")
