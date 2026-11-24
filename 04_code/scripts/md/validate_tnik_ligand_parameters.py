from pathlib import Path
import json
import subprocess

ROOT = Path.cwd()
WORK = ROOT / "06_docking/molecular_dynamics/ligand_parameterization"
LOG = ROOT / "06_docking/molecular_dynamics/logs"
LOG.mkdir(parents=True, exist_ok=True)

results = {}

for name in ["ponatinib", "tepotinib"]:
    print(f"\n===== {name.upper()} =====")
    folder = WORK / name
    mol2 = folder / f"{name}.mol2"
    frcmod = folder / f"{name}.frcmod"

    if not mol2.exists() or not frcmod.exists():
        raise FileNotFoundError(f"Missing parameter files for {name}")

    lines = mol2.read_text().splitlines()
    atom_section = False
    atoms = []
    charge_sum = 0.0

    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            atom_section = True
            continue
        if line.startswith("@<TRIPOS>") and atom_section:
            atom_section = False
        if atom_section and line.strip():
            parts = line.split()
            if len(parts) >= 9:
                atoms.append(parts)
                charge_sum += float(parts[8])

    missing = [
        line.strip() for line in frcmod.read_text().splitlines()
        if "ATTN" in line.upper()
    ]

    results[name] = {
        "mol2_atoms": len(atoms),
        "charge_sum": charge_sum,
        "frcmod_bytes": frcmod.stat().st_size,
        "missing_parameter_warnings": missing,
        "atom_types": sorted(set(
            atom[5] for atom in atoms
        ))
    }

    print("MOL2 atoms:", len(atoms))
    print("Total partial charge:", round(charge_sum, 4))
    print("FRCMOD size:", frcmod.stat().st_size, "bytes")
    print("Missing parameter warnings:", len(missing))
    print("Atom types:", ", ".join(results[name]["atom_types"]))

path = LOG / "TNIK_ligand_parameter_validation.json"
path.write_text(json.dumps(results, indent=2))

print("\n================================")
print("PARAMETER VALIDATION INSPECTION COMPLETE")
print("Report:", path)
print("No protein or simulation files modified.")
