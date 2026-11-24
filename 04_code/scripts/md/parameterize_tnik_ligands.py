from pathlib import Path
import subprocess
import json

ROOT = Path.cwd()
MD = ROOT / "06_docking" / "molecular_dynamics"
LIG = MD / "ligands"
LOG = MD / "logs"
WORK = MD / "ligand_parameterization"

LOG.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)

summary = {}

def run(command, cwd, logfile):
    result = subprocess.run(
        command, cwd=cwd, capture_output=True, text=True
    )
    logfile.write_text(
        "COMMAND:\n" + " ".join(command)
        + "\n\nSTDOUT:\n" + result.stdout
        + "\n\nSTDERR:\n" + result.stderr
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed. See {logfile}"
        )

for name in ["ponatinib", "tepotinib"]:
    print(f"\n===== {name.upper()} =====")

    source = LIG / f"{name}_docked_first_pose_obabel.sdf"
    folder = WORK / name
    folder.mkdir(parents=True, exist_ok=True)

    if not source.exists():
        raise FileNotFoundError(source)

    copied = folder / f"{name}.sdf"
    copied.write_bytes(source.read_bytes())

    mol2 = folder / f"{name}.mol2"
    frcmod = folder / f"{name}.frcmod"

    print("Running Antechamber...")

    run([
        "antechamber",
        "-i", str(copied),
        "-fi", "sdf",
        "-o", str(mol2),
        "-fo", "mol2",
        "-at", "gaff2",
        "-c", "bcc",
        "-nc", "0",
        "-s", "2"
    ], folder, LOG / f"{name}_antechamber.log")

    print("Running Parmchk2...")

    run([
        "parmchk2",
        "-i", str(mol2),
        "-f", "mol2",
        "-o", str(frcmod),
        "-s", "2"
    ], folder, LOG / f"{name}_parmchk2.log")

    if not mol2.exists() or not frcmod.exists():
        raise RuntimeError(f"{name}: parameter files missing")

    summary[name] = {
        "status": "PARAMETER_FILES_GENERATED",
        "mol2": str(mol2),
        "frcmod": str(frcmod),
        "mol2_bytes": mol2.stat().st_size,
        "frcmod_bytes": frcmod.stat().st_size
    }

    print("GAFF2 MOL2 generated.")
    print("FRCMOD generated.")

out = LOG / "TNIK_ligand_parameterization_summary.json"
out.write_text(json.dumps(summary, indent=2))

print("\n================================")
print("LIGAND PARAMETERIZATION FINISHED")
print("================================")
print("Summary:", out)
print("These are provisional neutral-state parameters.")
print("Final protonation and stereochemistry validation is still required.")
