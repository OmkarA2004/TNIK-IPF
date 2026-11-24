from pathlib import Path
import shutil
import parmed as pmd

base = Path("06_docking/molecular_dynamics")
systems = ["ponatinib_N37_H", "tepotinib_pH74_plus1"]

for name in systems:
    source = base / "complex_systems" / name
    target = base / "energy_minimization" / name

    parm = pmd.load_file(
        str(source / "complex.prmtop"),
        str(source / "complex.inpcrd")
    )

    mol = [r for r in parm.residues if r.name == "MOL"]
    if len(mol) != 1:
        raise RuntimeError(f"{name}: expected exactly one MOL residue")

    ligand = mol[0]
    old_charge = sum(a.charge for a in ligand.atoms)
    correction = 1.0 - old_charge

    # Find a hydrogen bonded to a positively charged nitrogen.
    candidates = []
    for atom in ligand.atoms:
        if atom.element != "H":
            continue
        for bond in atom.bonds:
            other = bond.atom1 if bond.atom2 is atom else bond.atom2
            if other.element == "N" and other.charge > 0:
                candidates.append((atom, other))

    if not candidates:
        raise RuntimeError(f"{name}: protonated N-H bond not found")

    hydrogen, nitrogen = max(
        candidates, key=lambda pair: pair[1].charge
    )

    # Preserve the ligand's intended +1 net charge.
    hydrogen.charge += correction

    new_charge = sum(a.charge for a in ligand.atoms)

    top = target / "complex.top"
    backup = target / "complex.top.before_charge_rounding"

    if not backup.exists():
        shutil.copy2(top, backup)

    parm.save(str(top), format="gromacs")

    print(f"\n{name}")
    print(f"Previous ligand charge: {old_charge:.6f}")
    print(f"Correction applied to H atom: {hydrogen.name}")
    print(f"Bonded nitrogen: {nitrogen.name}")
    print(f"Correction: {correction:+.6f}")
    print(f"New ligand charge: {new_charge:.6f}")
    print("Original GROMACS topology backup saved.")

print("\nCharge-rounding correction completed.")
