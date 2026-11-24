from pathlib import Path
import subprocess
import shutil

ROOT = Path.cwd()
BASE = ROOT / "06_docking/molecular_dynamics"
OUT = BASE / "complex_systems"

PROTEIN = BASE / "structures/TNIK_6RA7_protein_raw.pdb"

SYSTEMS = {
    "ponatinib_N32_H": (
        BASE / "ligand_parameterization/pH74_GAFF2/ponatinib_N32_H/ponatinib_N32_H.mol2",
        BASE / "ligand_parameterization/pH74_GAFF2/ponatinib_N32_H/ponatinib_N32_H.frcmod",
    ),
    "ponatinib_N37_H": (
        BASE / "ligand_parameterization/pH74_GAFF2/ponatinib_N37_H/ponatinib_N37_H.mol2",
        BASE / "ligand_parameterization/pH74_GAFF2/ponatinib_N37_H/ponatinib_N37_H.frcmod",
    ),
    "tepotinib_pH74_plus1": (
        BASE / "ligand_parameterization/pH74_GAFF2/tepotinib_pH74_plus1/tepotinib_pH74_plus1.mol2",
        BASE / "ligand_parameterization/pH74_GAFF2/tepotinib_pH74_plus1/tepotinib_pH74_plus1.frcmod",
    ),
}

for required in [PROTEIN]:
    if not required.is_file():
        raise FileNotFoundError(required)

for name, (mol2, frcmod) in SYSTEMS.items():
    if not mol2.is_file() or not frcmod.is_file():
        raise FileNotFoundError(f"Missing ligand parameters for {name}")

OUT.mkdir(parents=True, exist_ok=True)

if shutil.which("tleap") is None:
    raise RuntimeError("tleap is not available in the active environment.")

for name, (mol2, frcmod) in SYSTEMS.items():
    folder = OUT / name
    folder.mkdir(parents=True, exist_ok=True)

    # Never overwrite an earlier build.
    if any(folder.iterdir()):
        raise RuntimeError(
            f"{folder} is not empty. Inspect it before rerunning; "
            "this script will not overwrite existing files."
        )

    leap_file = folder / "build.leap"
    log_file = folder / "tleap.log"

    leap_text = f"""
source leaprc.protein.ff14SB
source leaprc.phosaa14SB
source leaprc.gaff2
source leaprc.water.tip3p

loadamberparams frcmod.ionsjc_tip3p
loadamberparams {frcmod.resolve()}

protein = loadpdb {PROTEIN.resolve()}
ligand = loadmol2 {mol2.resolve()}

check protein
check ligand

complex = combine {{ protein ligand }}

check complex
charge complex

solvatebox complex TIP3PBOX 12.0

# Neutralise the complex with chloride ions.
addIonsRand complex Cl- 0

check complex
charge complex

saveamberparm complex { (folder / "complex.prmtop").resolve() } { (folder / "complex.inpcrd").resolve() }
savepdb complex { (folder / "complex.pdb").resolve() }

quit
"""
    leap_file.write_text(leap_text)

    print(f"\n===== Building {name} =====")
    result = subprocess.run(
        ["tleap", "-f", str(leap_file)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    log_file.write_text(result.stdout)

    print("tleap return code:", result.returncode)
    print("Log:", log_file)

    if result.returncode != 0:
        print(result.stdout[-5000:])
        raise RuntimeError(f"LEaP failed for {name}")

    if "Errors = 0" not in result.stdout:
        print("WARNING: Inspect LEaP log for errors or warnings.")

    expected = [
        folder / "complex.prmtop",
        folder / "complex.inpcrd",
        folder / "complex.pdb",
    ]
    missing = [str(p) for p in expected if not p.is_file()]
    if missing:
        print(result.stdout[-5000:])
        raise RuntimeError(f"Missing output files: {missing}")

    print("Topology and coordinate files created.")
    for p in expected:
        print(f"  {p.name}: {p.stat().st_size:,} bytes")

print("\nAll three LEaP build attempts completed.")
print("Inspect all tleap.log files before proceeding to GROMACS.")
