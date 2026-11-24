from pathlib import Path
import subprocess

root = Path.cwd()
out = root / "06_docking/interaction_analysis/PLIP_complexes"
out.mkdir(parents=True, exist_ok=True)

receptor = (
    root / "06_docking/interaction_analysis/conversion_check/"
    "TNIK_receptor_test.pdb"
)

systems = {
    "ponatinib_seed11": (
        root / "06_docking/reproducibility/ponatinib_seed_11_out.pdbqt",
        40,
    ),
    "tepotinib_seed22": (
        root / "06_docking/reproducibility/tepotinib_seed_22_out.pdbqt",
        37,
    ),
}

if not receptor.exists():
    raise FileNotFoundError(f"Missing receptor: {receptor}")

def extract_model_one(source, destination):
    lines = source.read_text(errors="replace").splitlines()
    inside = False
    found = False
    atoms = []

    for line in lines:
        if line.startswith("MODEL"):
            model_number = line.split()[-1]
            inside = model_number == "1"
            found = found or inside
            continue

        if line.startswith("ENDMDL"):
            if inside:
                break
            inside = False
            continue

        if inside and line.startswith(("ATOM  ", "HETATM")):
            atoms.append(line)

    if not found:
        raise ValueError(f"MODEL 1 not found in {source}")

    if not atoms:
        raise ValueError(f"No atom records found in MODEL 1 of {source}")

    destination.write_text("\n".join(atoms) + "\n")
    return len(atoms)

for name, (source, expected_atoms) in systems.items():
    print(f"\nPreparing {name}")

    pose_pdbqt = out / f"{name}_model1.pdbqt"
    pose_pdb = out / f"{name}_model1.pdb"
    complex_pdb = out / f"{name}_complex.pdb"

    count = extract_model_one(source, pose_pdbqt)

    if count != expected_atoms:
        raise ValueError(
            f"{name}: expected {expected_atoms} atoms, found {count}"
        )

    subprocess.run(
        ["obabel", str(pose_pdbqt), "-O", str(pose_pdb)],
        check=True,
        capture_output=True,
        text=True,
    )

    ligand_lines = [
        line for line in pose_pdb.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]

    if len(ligand_lines) != expected_atoms:
        raise ValueError(
            f"{name}: converted atom count {len(ligand_lines)} "
            f"does not match expected {expected_atoms}"
        )

    # Make the ligand a separate HETATM residue named MOL.
    formatted_ligand = []
    for line in ligand_lines:
        line = line.ljust(80)
        line = "HETATM" + line[6:]
        line = line[:17] + "MOL" + line[20:]
        line = line[:21] + "Z" + line[22:]
        line = line[:22] + f"{1:4d}" + line[26:]
        formatted_ligand.append(line.rstrip())

    receptor_lines = [
        line for line in receptor.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]

    # Keep receptor coordinates and numbering unchanged.
    complex_lines = receptor_lines + formatted_ligand + ["TER", "END"]
    complex_pdb.write_text("\n".join(complex_lines) + "\n")

    # Validate the assembled complex.
    protein_atoms = sum(
        line.startswith("ATOM  ") for line in receptor_lines
    )
    ligand_atoms = sum(
        line.startswith("HETATM") and line[17:20] == "MOL"
        for line in complex_lines
    )
    ligand_residues = {
        (line[21], line[22:26], line[26])
        for line in complex_lines
        if line.startswith("HETATM") and line[17:20] == "MOL"
    }

    if ligand_atoms != expected_atoms or len(ligand_residues) != 1:
        raise ValueError(f"Complex validation failed for {name}")

    print("  Model extracted: 1")
    print("  Ligand atoms:", ligand_atoms)
    print("  Receptor atoms:", protein_atoms)
    print("  Ligand residues:", len(ligand_residues))
    print("  Complex:", complex_pdb)
    print("  Validation: PASSED")

print("\nAll representative docking complexes prepared.")
