from pathlib import Path
import json
import numpy as np
from rdkit import Chem

ROOT = Path("06_docking/molecular_dynamics")
LIG = ROOT / "ligands"
TEST = ROOT / "ligand_parameterization" / "pH74_test"
OUT = ROOT / "ligand_parameterization" / "pH74_candidates"
OUT.mkdir(parents=True, exist_ok=True)

def load_sdf(path):
    mol = Chem.SDMolSupplier(str(path), removeHs=False)[0]
    if mol is None or mol.GetNumConformers() == 0:
        raise RuntimeError(f"Cannot read molecule/conformer: {path}")
    return mol

def heavy_atoms(mol):
    return [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]

def check_coordinates(candidate, reference, label):
    ca, ra = heavy_atoms(candidate), heavy_atoms(reference)
    if len(ca) != len(ra):
        raise RuntimeError(f"{label}: heavy-atom count mismatch")
    if [a.GetSymbol() for a in ca] != [a.GetSymbol() for a in ra]:
        raise RuntimeError(f"{label}: heavy-atom order/elements mismatch")

    cc, rc = candidate.GetConformer(), reference.GetConformer()
    errors = []
    for a, b in zip(ca, ra):
        p = np.array(cc.GetAtomPosition(a.GetIdx()))
        q = np.array(rc.GetAtomPosition(b.GetIdx()))
        errors.append(np.linalg.norm(p-q))

    maximum = max(errors, default=0.0)
    if maximum > 0.001:
        raise RuntimeError(
            f"{label}: heavy-atom coordinates changed; max difference {maximum:.6f} Å"
        )
    return maximum

def write_sdf(mol, path):
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()

# Reference docked poses
pon_ref = load_sdf(LIG / "ponatinib_docked_first_pose_obabel.sdf")
tep_ref = load_sdf(LIG / "tepotinib_docked_first_pose_obabel.sdf")

# Ponatinib doubly protonated Open Babel structure
pon_path = TEST / "ponatinib_pH74_only.sdf"
pon = load_sdf(pon_path)

if Chem.GetFormalCharge(pon) != 2:
    raise RuntimeError("Expected the Ponatinib starting structure to have charge +2")

# Docked-pose atom indices from the validated heavy-atom mapping.
# N32 = original ligand atom 33; N37 = original ligand atom 36.
nitrogen_indices = {"N32_H": 32, "N37_H": 37}

for keep_protonated, keep_idx in nitrogen_indices.items():
    mol = Chem.Mol(pon)
    rw = Chem.RWMol(mol)

    # Neutralize the OTHER piperazine nitrogen by removing its attached H.
    other_idx = 37 if keep_idx == 32 else 32
    nitrogen = rw.GetAtomWithIdx(other_idx)

    attached_h = [
        n.GetIdx() for n in nitrogen.GetNeighbors()
        if n.GetSymbol() == "H"
    ]
    if len(attached_h) != 1 or nitrogen.GetFormalCharge() != 1:
        raise RuntimeError(
            f"Expected one explicit H on protonated nitrogen {other_idx}"
        )

    rw.RemoveAtom(attached_h[0])
    # Atom indices above are original indices; removing H does not alter
    # the heavy-atom indices preceding it in this input structure.
    neutral_nitrogen = rw.GetAtomWithIdx(other_idx)
    neutral_nitrogen.SetFormalCharge(0)
    neutral_nitrogen.SetNoImplicit(True)

    result = rw.GetMol()
    Chem.SanitizeMol(result)

    charge = Chem.GetFormalCharge(result)
    if charge != 1:
        raise RuntimeError(f"{keep_protonated}: expected charge +1, got {charge}")

    max_diff = check_coordinates(
        result, pon_ref, f"Ponatinib {keep_protonated}"
    )
    filename = OUT / f"ponatinib_{keep_protonated}.sdf"
    write_sdf(result, filename)

    print(f"\n{filename}")
    print("Formal charge:", charge)
    print("Heavy atoms:", len(heavy_atoms(result)))
    print("Max heavy-atom coordinate change (Å):", f"{max_diff:.8f}")
    print("SMILES:", Chem.MolToSmiles(Chem.RemoveHs(result)))

# Validate Tepotinib +1 candidate against original docked pose.
tep = load_sdf(TEST / "tepotinib_pH74_only.sdf")
if Chem.GetFormalCharge(tep) != 1:
    raise RuntimeError("Tepotinib candidate is not +1")

tep_diff = check_coordinates(tep, tep_ref, "Tepotinib")
write_sdf(tep, OUT / "tepotinib_pH74_plus1.sdf")

print(f"\n{OUT / 'tepotinib_pH74_plus1.sdf'}")
print("Formal charge:", Chem.GetFormalCharge(tep))
print("Heavy atoms:", len(heavy_atoms(tep)))
print("Max heavy-atom coordinate change (Å):", f"{tep_diff:.8f}")

report = {
    "pH_assumption": 7.4,
    "ponatinib_candidates": [
        "ponatinib_N32_H.sdf",
        "ponatinib_N37_H.sdf"
    ],
    "tepotinib_candidate": "tepotinib_pH74_plus1.sdf",
    "status": "candidate structures generated; protonation microstate selection and parameter validation still required"
}
(OUT / "candidate_report.json").write_text(json.dumps(report, indent=2))
print("\nCANDIDATE GENERATION COMPLETE")
