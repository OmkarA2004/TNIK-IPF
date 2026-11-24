
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Project paths
project_dir = Path.home() / "Desktop" / "TNIK_IPF_Project"

input_file = (
    project_dir / "06_docking" / "receptor_prep"
    / "TNIK_6RA7_chainA_altA.pdb"
)

output_file = (
    project_dir / "06_docking" / "receptor_prep"
    / "TNIK_6RA7_repaired.pdb"
)

# Load the cleaned TNIK structure
fixer = PDBFixer(filename=str(input_file))

# Identify missing residues and atoms
fixer.findMissingResidues()
fixer.missingResidues = {}

fixer.findNonstandardResidues()

# Keep standard amino acids unchanged
fixer.findMissingAtoms()

print("Missing atoms identified:", len(fixer.missingAtoms))

for residue, atoms in fixer.missingAtoms.items():
    print(
        residue.id,
        residue.name,
        "Missing:",
        [atom.name for atom in atoms]
    )

# Add missing atoms
fixer.addMissingAtoms()

# Save repaired structure
with open(output_file, "w") as output:
    PDBFile.writeFile(
        fixer.topology,
        fixer.positions,
        output,
        keepIds=True
    )

print("\nRepaired structure saved:")
print(output_file)

