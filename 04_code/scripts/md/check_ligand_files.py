from pathlib import Path

folders = [
    Path("06_docking/molecular_dynamics/structures"),
    Path("06_docking/molecular_dynamics/ligands"),
    Path("06_docking/virtual_screen"),
    Path("06_docking/redocking"),
    Path("06_docking/ligand_prep"),
    Path("06_docking/receptor_prep"),
]

keywords = ("ponatinib", "tepotinib", "tnik", "6ra7")
extensions = {".pdb", ".pdbqt", ".sdf", ".mol2", ".mol", ".gro", ".top"}

print("LIGAND AND PROTEIN FILE INVENTORY")
print("=" * 45)

found = set()

for folder in folders:
    if not folder.exists():
        continue

    for file in folder.rglob("*"):
        if file.is_file() and file.suffix.lower() in extensions:
            if any(k in file.name.lower() for k in keywords):
                found.add(file)

if found:
    for file in sorted(found):
        print(f"{file}  ({file.stat().st_size:,} bytes)")
else:
    print("No matching files found in the checked folders.")

print("\nInventory complete.")
