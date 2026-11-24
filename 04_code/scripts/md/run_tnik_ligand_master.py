from pathlib import Path
import subprocess
import shutil
import sys
import json
from datetime import datetime

from rdkit import Chem

ROOT = Path.home() / "Desktop" / "TNIK_IPF_Project"
MD = ROOT / "06_docking" / "molecular_dynamics"
LIG = MD / "ligands"
LOG = MD / "logs"
RESULTS = ROOT / "09_results" / "molecular_dynamics"

LOG.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

report = {
    "started": datetime.now().isoformat(),
    "project": str(ROOT),
    "software": {},
    "ligands": {},
    "inputs": {},
    "warnings": [],
    "status": "INPUT_VALIDATION"
}

def check_tool(name):
    path = shutil.which(name)
    report["software"][name] = path or "NOT FOUND"
    return path is not None

def check_file(label, path):
    exists = path.exists()
    report["inputs"][label] = {
        "path": str(path),
        "exists": exists
    }
    return exists

print("\n==========================================")
print(" TNIK LIGAND-BOUND MD MASTER WORKFLOW")
print("==========================================")

print("\n[1/5] Checking software...")

tools = [
    "gmx", "obabel", "antechamber",
    "parmchk2", "tleap", "acpype"
]

for tool in tools:
    available = check_tool(tool)
    print(f"  {tool}: {'READY' if available else 'MISSING'}")

print("\n[2/5] Checking protein inputs...")

protein = MD / "structures" / "TNIK_6RA7_protein_raw.pdb"
apo_tpr = MD / "parameters" / "TNIK_npt_600ps.tpr"
apo_xtc = MD / "parameters" / "TNIK_npt_crescale.xtc"

check_file("TNIK protein structure", protein)
check_file("Apo NPT TPR", apo_tpr)
check_file("Apo NPT trajectory", apo_xtc)

print("\n[3/5] Checking ligand structures...")

for name, seed in [("ponatinib", 11), ("tepotinib", 22)]:
    sdf = LIG / f"{name}_original.sdf"
    pose = LIG / f"{name}_docked_first_pose.pdbqt"

    ligand = {
        "sdf": str(sdf),
        "pose": str(pose),
        "sdf_valid": False,
        "pose_valid": False
    }

    if sdf.exists():
        supplier = Chem.SDMolSupplier(str(sdf), removeHs=False)
        mol = next((m for m in supplier if m is not None), None)

        if mol is not None:
            ligand["sdf_valid"] = True
            ligand["atoms"] = mol.GetNumAtoms()
            ligand["heavy_atoms"] = mol.GetNumHeavyAtoms()
            ligand["bonds"] = mol.GetNumBonds()
            ligand["formal_charge"] = Chem.GetFormalCharge(mol)
            ligand["canonical_smiles"] = Chem.MolToSmiles(
                Chem.RemoveHs(mol)
            )

    if pose.exists():
        atoms = [
            line for line in pose.read_text().splitlines()
            if line.startswith(("ATOM", "HETATM"))
        ]
        ligand["pose_atoms"] = len(atoms)
        ligand["pose_valid"] = len(atoms) > 0

    report["ligands"][name] = ligand

    print(f"\n  {name.upper()}")
    print(f"  Original SDF valid: {ligand['sdf_valid']}")
    print(f"  Docked pose valid: {ligand['pose_valid']}")
    print(f"  SDF heavy atoms: {ligand.get('heavy_atoms', 'N/A')}")
    print(f"  Pose atoms: {ligand.get('pose_atoms', 'N/A')}")

    if ligand["sdf_valid"] and ligand["pose_valid"]:
        if ligand["heavy_atoms"] != ligand["pose_atoms"]:
            report["warnings"].append(
                f"{name}: heavy atom count differs from pose atom count"
            )

print("\n[4/5] Checking force-field parameter availability...")

phos = check_file(
    "Phosphothreonine force-field library",
    ROOT / "06_docking" / "molecular_dynamics" /
    "parameters" / "phosaa14SB.dat"
)

print("\n[5/5] Writing validation report...")

report["finished"] = datetime.now().isoformat()

required_tools = all(
    report["software"][x] != "NOT FOUND" for x in tools
)

required_ligands = all(
    report["ligands"].get(x, {}).get("sdf_valid", False)
    and report["ligands"].get(x, {}).get("pose_valid", False)
    for x in ["ponatinib", "tepotinib"]
)

if required_tools and required_ligands:
    report["status"] = "VALIDATION_COMPLETED_NEEDS_CHEMICAL_MAPPING"
else:
    report["status"] = "INPUTS_NEED_ATTENTION"

report_path = LOG / "TNIK_ligand_master_validation.json"
report_path.write_text(json.dumps(report, indent=2))

print("\n==========================================")
print(" MASTER WORKFLOW VALIDATION FINISHED")
print("==========================================")
print("Status:", report["status"])
print("Detailed report:", report_path)
print("\nNo original structures or apo results were modified.")
print("No MD simulation was started.")
