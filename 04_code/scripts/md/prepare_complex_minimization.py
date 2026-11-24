from pathlib import Path
import parmed as pmd

base = Path("06_docking/molecular_dynamics")
systems = [
    "ponatinib_N32_H",
    "ponatinib_N37_H",
    "tepotinib_pH74_plus1"
]

for name in systems:
    source = base / "complex_systems" / name
    target = base / "energy_minimization" / name
    target.mkdir(parents=True, exist_ok=True)

    print(f"\nPreparing {name}...")

    parm = pmd.load_file(
        str(source / "complex.prmtop"),
        str(source / "complex.inpcrd")
    )

    parm.save(str(target / "complex.top"), format="gromacs")
    parm.save(str(target / "complex.gro"), format="gro")

    print("Topology:", (target / "complex.top").stat().st_size, "bytes")
    print("Coordinates:", (target / "complex.gro").stat().st_size, "bytes")

print("\nAll three conversions completed.")
