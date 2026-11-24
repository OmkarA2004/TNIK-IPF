from pathlib import Path
import parmed as pmd
import subprocess
import json

BASE = Path("06_docking/molecular_dynamics")
SYSTEMS = [
    "ponatinib_N32_H",
    "ponatinib_N37_H",
    "tepotinib_pH74_plus1",
]

report = {
    "project": "TNIK ligand complex MD",
    "systems": {},
    "overall": "CHECKS_IN_PROGRESS"
}

for name in SYSTEMS:
    folder = BASE / "energy_minimization" / name
    source = BASE / "complex_systems" / name

    result = {
        "complex_topology": (source / "complex.prmtop").exists(),
        "complex_coordinates": (source / "complex.inpcrd").exists(),
        "gromacs_topology": (folder / "complex.top").exists(),
        "gromacs_coordinates": (folder / "complex.gro").exists(),
        "minimization_tpr": (folder / "minimization.tpr").exists(),
        "minimized_coordinates": (folder / "minimization.gro").exists(),
        "nvt_parameters": (folder / "nvt.mdp").exists(),
        "position_restraint_file": (folder / "posre_protein.itp").exists(),
    }

    top = folder / "complex.top"
    if top.exists():
        content = top.read_text(errors="replace")
        result["topology_has_posres_reference"] = (
            "posre_protein.itp" in content
            or "posre" in content.lower()
        )
        result["topology_include_lines"] = [
            line.strip() for line in content.splitlines()
            if line.strip().startswith("#include")
        ]

    log = folder / "minimization_run.log"
    if log.exists():
        text = log.read_text(errors="replace")
        result["minimization_converged"] = (
            "Steepest Descents converged" in text
        )
        result["minimization_log_exists"] = True
    else:
        result["minimization_converged"] = False

    if result["complex_topology"] and result["complex_coordinates"]:
        try:
            parm = pmd.load_file(
                str(source / "complex.prmtop"),
                str(source / "complex.inpcrd")
            )
            result["atom_count"] = len(parm.atoms)
            result["residue_count"] = len(parm.residues)
            result["total_charge"] = round(
                sum(atom.charge for atom in parm.atoms), 6
            )
        except Exception as exc:
            result["parameter_read_error"] = str(exc)

    result["ready_for_nvt_validation"] = all([
        result["gromacs_topology"],
        result["gromacs_coordinates"],
        result["minimization_tpr"],
        result["minimized_coordinates"],
        result["nvt_parameters"],
        result["position_restraint_file"],
        result.get("topology_has_posres_reference", False),
        result.get("minimization_converged", False),
    ])

    report["systems"][name] = result

    print("\n", name)
    for key, value in result.items():
        print(f"{key}: {value}")

report["overall"] = (
    "READY_FOR_NVT_VALIDATION"
    if all(x["ready_for_nvt_validation"]
           for x in report["systems"].values())
    else "NOT_READY_CHECK_REPORTED_ITEMS"
)

out = BASE / "analysis" / "master_md_readiness.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, indent=2))

print("\n======================================")
print("OVERALL:", report["overall"])
print("Report saved:", out)
print("======================================")
