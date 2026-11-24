from pathlib import Path
from rdkit import Chem
import subprocess
import json

base = Path("06_docking/molecular_dynamics/ligands")
report_dir = Path("06_docking/molecular_dynamics/logs")
report_dir.mkdir(parents=True, exist_ok=True)

report = {}

for name in ["ponatinib", "tepotinib"]:
    print(f"\nInvestigating {name.upper()}...")

    sdf = base / f"{name}_original.sdf"
    pdbqt = base / f"{name}_docked_first_pose.pdbqt"
    mol2 = base / f"{name}_docked_first_pose_obabel.mol2"
    converted_sdf = base / f"{name}_docked_first_pose_obabel.sdf"

    original = next(
        (m for m in Chem.SDMolSupplier(str(sdf), removeHs=False)
         if m is not None), None
    )

    result = {
        "original_sdf_readable": original is not None,
        "pose_exists": pdbqt.exists(),
        "conversion": {},
    }

    if original is None or not pdbqt.exists():
        result["error"] = "Missing or unreadable input"
        report[name] = result
        continue

    for fmt, output in [("mol2", mol2), ("sdf", converted_sdf)]:
        process = subprocess.run(
            ["obabel", str(pdbqt), "-O", str(output)],
            capture_output=True, text=True
        )

        result["conversion"][fmt] = {
            "returncode": process.returncode,
            "output_exists": output.exists(),
            "message": (process.stderr + process.stdout)[-1500:]
        }

    converted = None
    if converted_sdf.exists():
        converted = next(
            (m for m in Chem.SDMolSupplier(
                str(converted_sdf), removeHs=False
            ) if m is not None), None
        )

    if converted is not None:
        result["converted_sdf_readable"] = True
        result["converted_atoms"] = converted.GetNumAtoms()
        result["converted_heavy_atoms"] = converted.GetNumHeavyAtoms()
        result["converted_bonds"] = converted.GetNumBonds()
        result["converted_smiles"] = Chem.MolToSmiles(
            Chem.RemoveHs(converted)
        )

        original_heavy = Chem.RemoveHs(original)
        converted_heavy = Chem.RemoveHs(converted)

        result["same_element_counts"] = (
            sorted(a.GetSymbol() for a in original_heavy.GetAtoms())
            == sorted(a.GetSymbol() for a in converted_heavy.GetAtoms())
        )

        result["same_graph_exact"] = (
            original_heavy.GetNumAtoms() == converted_heavy.GetNumAtoms()
            and original_heavy.GetNumBonds() == converted_heavy.GetNumBonds()
            and original_heavy.HasSubstructMatch(converted_heavy)
            and converted_heavy.HasSubstructMatch(original_heavy)
        )
    else:
        result["converted_sdf_readable"] = False

    report[name] = result
    print("  Original heavy atoms:", original.GetNumHeavyAtoms())
    print("  Converted pose readable:", result["converted_sdf_readable"])
    print("  Same element counts:", result.get("same_element_counts", "N/A"))
    print("  Exact graph match:", result.get("same_graph_exact", "NOT ESTABLISHED"))

path = report_dir / "TNIK_ligand_mapping_investigation.json"
path.write_text(json.dumps(report, indent=2))

print("\nInvestigation complete.")
print("Detailed report:", path)
print("Original SDF and PDBQT files were not modified.")
