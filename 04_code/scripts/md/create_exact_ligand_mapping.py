from pathlib import Path
from rdkit import Chem
import json

base = Path("06_docking/molecular_dynamics/ligands")
log = Path("06_docking/molecular_dynamics/logs")
log.mkdir(parents=True, exist_ok=True)

results = {}

for name in ["ponatinib", "tepotinib"]:
    original_path = base / f"{name}_original.sdf"
    pose_path = base / f"{name}_docked_first_pose_obabel.sdf"

    original = next(
        m for m in Chem.SDMolSupplier(
            str(original_path), removeHs=False
        ) if m is not None
    )
    pose = next(
        m for m in Chem.SDMolSupplier(
            str(pose_path), removeHs=False
        ) if m is not None
    )

    original = Chem.RemoveHs(original)
    pose = Chem.RemoveHs(pose)

    # Match the complete heavy-atom molecular graphs.
    original_match = original.GetSubstructMatch(
        pose, useChirality=True
    )
    pose_match = pose.GetSubstructMatch(
        original, useChirality=True
    )

    valid = (
        len(original_match) == pose.GetNumAtoms()
        and len(pose_match) == original.GetNumAtoms()
    )

    mapping = []
    if valid:
        for pose_index, original_index in enumerate(original_match):
            mapping.append({
                "pose_atom_index": pose_index,
                "original_sdf_atom_index": original_index,
                "element": pose.GetAtomWithIdx(
                    pose_index
                ).GetSymbol()
            })

    results[name] = {
        "mapping_valid": valid,
        "original_heavy_atoms": original.GetNumAtoms(),
        "pose_heavy_atoms": pose.GetNumAtoms(),
        "mapping": mapping
    }

    print(f"\n{name.upper()}")
    print("Mapping valid:", valid)
    print("Mapped atoms:", len(mapping))

    if not valid:
        print("WARNING: Exact mapping not established.")

    (base / f"{name}_atom_mapping.json").write_text(
        json.dumps(results[name], indent=2)
    )

report = log / "TNIK_exact_ligand_atom_mapping.json"
report.write_text(json.dumps(results, indent=2))

print("\nMAPPING INVESTIGATION FINISHED")
print("Report:", report)
