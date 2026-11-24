from pathlib import Path
import subprocess
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

root = Path("06_docking/molecular_dynamics")
source = root / "ligand_parameterization" / "pH74_candidates"
output = root / "ligand_parameterization" / "pH74_GAFF2"
output.mkdir(parents=True, exist_ok=True)

candidates = {
    "ponatinib_N32_H": 1,
    "ponatinib_N37_H": 1,
    "tepotinib_pH74_plus1": 1,
}

def run(command, log):
    print("Running:", " ".join(map(str, command)))
    with open(log, "w") as f:
        result = subprocess.run(
            list(map(str, command)),
            stdout=f,
            stderr=subprocess.STDOUT
        )
    if result.returncode:
        print(f"ERROR: command failed. See {log}")
        sys.exit(result.returncode)

for name, expected_charge in candidates.items():
    print(f"\n========== {name} ==========")
    folder = output / name
    folder.mkdir(exist_ok=True)

    sdf = source / f"{name}.sdf"
    mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    if mol is None:
        raise RuntimeError(f"Cannot read {sdf}")

    charge = Chem.GetFormalCharge(mol)
    heavy = mol.GetNumHeavyAtoms()
    if charge != expected_charge:
        raise RuntimeError(f"{name}: expected charge +1, got {charge}")

    mol2 = folder / f"{name}.mol2"
    frcmod = folder / f"{name}.frcmod"

    run([
        "antechamber",
        "-i", sdf, "-fi", "sdf",
        "-o", mol2, "-fo", "mol2",
        "-c", "bcc",
        "-nc", str(expected_charge),
        "-at", "gaff2",
        "-s", "2"
    ], folder / "antechamber.log")

    run([
        "parmchk2",
        "-i", mol2, "-f", "mol2",
        "-o", frcmod,
        "-s", "2"
    ], folder / "parmchk2.log")

    # Read MOL2 atom count and partial charges.
    in_atoms = False
    atom_count = 0
    total_charge = 0.0
    with open(mol2) as f:
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>") and in_atoms:
                break
            if in_atoms and line.strip():
                fields = line.split()
                if len(fields) >= 9:
                    atom_count += 1
                    total_charge += float(fields[8])

    if atom_count != mol.GetNumAtoms():
        raise RuntimeError(
            f"{name}: atom count mismatch; SDF={mol.GetNumAtoms()}, MOL2={atom_count}"
        )

    if abs(total_charge - expected_charge) > 0.05:
        raise RuntimeError(
            f"{name}: partial charge sum {total_charge:.4f} "
            f"does not agree with formal charge +1"
        )

    if not frcmod.exists() or frcmod.stat().st_size == 0:
        raise RuntimeError(f"{name}: missing or empty frcmod")

    print("Formal charge:", charge)
    print("Heavy atoms:", heavy)
    print("SDF atoms:", mol.GetNumAtoms())
    print("MOL2 atoms:", atom_count)
    print("MOL2 partial-charge sum:", f"{total_charge:.4f}")
    print("FRCMOD size:", frcmod.stat().st_size, "bytes")
    print("Parameter files generated.")

print("\nALL THREE GAFF2 PARAMETER SETS GENERATED AND BASICALLY VALIDATED.")
print("Note: this is not yet full force-field or complex-topology validation.")
