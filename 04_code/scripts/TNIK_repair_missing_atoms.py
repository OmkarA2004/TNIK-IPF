
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Project directory
project_dir = Path.home() / "Desktop" / "TNIK_IPF_Project"

# Input: cleaned TNIK protein
input_file = (
    project_dir
    / "06_docking"
    / "receptor_prep"
    / "TNIK_6RA7_chainA_altA.pdb"
)

# Output: repaired TNIK protein
output_file = (
    project_dir
    / "06_docking"
    / "receptor_prep"
    / "TNIK_6RA7_repaired.pdb"
)

print("Loading TNIK structure...")

fixer = PDBFixer(filename=str(input_file))

# Do not add missing protein residues
fixer.findMissingResidues()
fixer.missingResidues = {}

# Identify nonstandard residues
fixer.findNonstandardResidues()

print("Nonstandard residues:", fixer.nonstandardResidues)

# Identify missing atoms
fixer.findMissingAtoms()

print("\nMissing atoms identified:", len(fixer.missingAtoms))

for residue, atoms in fixer.missingAtoms.items():
    print(
        "Residue:", residue.id,
        residue.name,
        "Missing atoms:",
        [atom.name for atom in atoms]
    )

# Add missing atoms to existing residues
print("\nAdding missing atoms...")
fixer.addMissingAtoms()

# Save repaired structure
with open(output_file, "w") as output:
    PDBFile.writeFile(
        fixer.topology,
        fixer.positions,
        output,
        keepIds=True
    )

print("\nRepaired structure saved successfully:")
print(output_file)
