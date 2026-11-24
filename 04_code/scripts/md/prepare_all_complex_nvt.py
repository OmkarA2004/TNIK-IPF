#!/usr/bin/env python3

"""
TNIK–ligand complex NVT preparation and validation.

Systems:
  1. ponatinib_N32_H
  2. ponatinib_N37_H
  3. tepotinib_pH74_plus1

Actions:
  - Validate topology and GRO atom counts.
  - Identify protein and ligand heavy atoms.
  - Generate molecule-local position restraints.
  - Back up original topologies.
  - Insert restraints into the correct molecule types.
  - Generate NVT MDP files.
  - Run GROMACS grompp validation.
  - Write JSON and text reports.

This script does NOT run MD.
It does NOT overwrite minimization coordinates or TPR files.
"""

from pathlib import Path
import json
import shutil
import subprocess
import sys
import traceback
from datetime import datetime


# ------------------------------------------------------------
# PROJECT CONFIGURATION
# ------------------------------------------------------------

PROJECT = Path.home() / "Desktop" / "TNIK_IPF_Project"

BASE = (
    PROJECT
    / "06_docking"
    / "molecular_dynamics"
    / "energy_minimization"
)

ANALYSIS = (
    PROJECT
    / "06_docking"
    / "molecular_dynamics"
    / "analysis"
)

LOG_DIR = (
    PROJECT
    / "06_docking"
    / "molecular_dynamics"
    / "logs"
)

SYSTEMS = [
    "ponatinib_N32_H",
    "ponatinib_N37_H",
    "tepotinib_pH74_plus1",
]

# Position-restraint force constant:
# kJ mol-1 nm-2
PROTEIN_FORCE = 1000
LIGAND_FORCE = 1000

# Temperature
TEMPERATURE = 300

# NVT duration
TIMESTEP_PS = 0.002
NVT_DURATION_PS = 500

NVT_STEPS = int(NVT_DURATION_PS / TIMESTEP_PS)

# GROMACS executable
GMX = shutil.which("gmx")

# Report files
ANALYSIS.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

REPORT_JSON = ANALYSIS / "complex_nvt_preparation_report.json"
REPORT_TXT = ANALYSIS / "complex_nvt_preparation_report.txt"


# ------------------------------------------------------------
# GENERAL UTILITIES
# ------------------------------------------------------------

def log(message):
    print(message, flush=True)


def run_command(command, cwd, log_file):
    """Run a command and save its complete output."""

    log_file.parent.mkdir(parents=True, exist_ok=True)

    with open(log_file, "w") as handle:
        result = subprocess.run(
            command,
            cwd=cwd,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )

    return result.returncode


def read_gro(path):
    """Read GRO coordinates and return atom records."""

    lines = path.read_text(errors="replace").splitlines()

    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file: {path}")

    total_atoms = int(lines[1].strip())

    atom_lines = lines[2:2 + total_atoms]

    if len(atom_lines) != total_atoms:
        raise ValueError(
            f"GRO atom count mismatch in {path}"
        )

    atoms = []

    for global_index, line in enumerate(atom_lines, start=1):

        if len(line) < 20:
            raise ValueError(
                f"Malformed GRO atom line {global_index}"
            )

        atoms.append({
            "global_index": global_index,
            "residue_number": int(line[0:5]),
            "residue_name": line[5:10].strip(),
            "atom_name": line[10:15].strip(),
            "atom_number": int(line[15:20]),
            "line": line,
        })

    return total_atoms, atoms


def parse_topology_molecules(topology_path):
    """
    Parse molecule definitions and atom counts.

    ParmEd-generated topology contains separate molecule types:
      system1
      MOL
      Cl-
      WAT
    """

    lines = topology_path.read_text(errors="replace").splitlines()

    molecules = {}
    current_section = None
    current_molecule = None

    for line in lines:

        stripped = line.strip()

        if not stripped or stripped.startswith(";") or stripped.startswith("#"):
            continue

        if stripped.startswith("[") and stripped.endswith("]"):

            current_section = (
                stripped.strip("[]").strip().lower()
            )

            continue

        fields = stripped.split()

        if current_section == "moleculetype":

            if len(fields) < 2:
                raise ValueError(
                    "Could not parse moleculetype entry: " + stripped
                )

            current_molecule = fields[0]

            if current_molecule in molecules:
                raise ValueError(
                    f"Duplicate molecule type: {current_molecule}"
                )

            molecules[current_molecule] = {
                "atom_count": 0,
                "start": None,
                "end": None,
            }

        elif current_section == "atoms" and current_molecule:

            if fields[0].isdigit():
                molecules[current_molecule]["atom_count"] += 1

        elif current_section == "system":

            current_molecule = None

    return molecules, lines


def parse_molecule_order(lines):
    """Read molecule counts from the [ molecules ] section."""

    current_section = None
    result = []

    for line in lines:

        stripped = line.strip()

        if not stripped or stripped.startswith(";") or stripped.startswith("#"):
            continue

        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped.strip("[]").strip().lower()
            continue

        if current_section == "molecules":

            fields = stripped.split()

            if len(fields) >= 2 and fields[1].isdigit():
                result.append((fields[0], int(fields[1])))

    return result


def is_hydrogen(atom_name):
    """
    Identify hydrogen atoms using the atom name.

    Handles names such as:
      H
      H1
      HA
      HH22
      1H
      2HG
    """

    name = atom_name.strip().upper()

    if not name:
        return False

    if name.startswith("H"):
        return True

    if name[0].isdigit() and len(name) > 1:
        return name[1] == "H"

    return False


def is_water_or_ion(residue_name):
    name = residue_name.strip().upper()

    return name in {
        "HOH",
        "WAT",
        "SOL",
        "TIP3",
        "TIP3P",
        "CL-",
        "NA+",
        "K+",
        "MG",
        "CA",
    }


# ------------------------------------------------------------
# GRO / TOPOLOGY VALIDATION
# ------------------------------------------------------------

def validate_system(folder):

    topology = folder / "complex.top"
    coordinates = folder / "complex.gro"
    minimization_tpr = folder / "minimization.tpr"
    minimized_coordinates = folder / "minimization.gro"
    original_nvt = folder / "nvt.mdp"

    required = [
        topology,
        coordinates,
        minimization_tpr,
        minimized_coordinates,
        original_nvt,
    ]

    for path in required:
        if not path.exists():
            raise FileNotFoundError(
                f"Required file missing: {path}"
            )

    molecule_types, topology_lines = parse_topology_molecules(topology)
    molecule_order = parse_molecule_order(topology_lines)

    total_atoms, atoms = read_gro(coordinates)

    required_molecule_types = {
        "system1",
        "MOL",
        "Cl-",
        "WAT",
    }

    missing = required_molecule_types - set(molecule_types)

    if missing:
        raise ValueError(
            f"Missing molecule types: {sorted(missing)}"
        )

    expected_protein_atoms = molecule_types["system1"]["atom_count"]
    expected_ligand_atoms = molecule_types["MOL"]["atom_count"]

    if expected_protein_atoms <= 0 or expected_ligand_atoms <= 0:
        raise ValueError("Protein or ligand atom count is zero.")

    # Calculate total atom count using [ molecules ].
    calculated_total = 0

    for molecule, copies in molecule_order:

        if molecule not in molecule_types:
            raise ValueError(
                f"Molecule {molecule} has no molecule type definition."
            )

        calculated_total += (
            molecule_types[molecule]["atom_count"] * copies
        )

    if calculated_total != total_atoms:
        raise ValueError(
            f"Topology/GRO total atom mismatch: "
            f"{calculated_total} vs {total_atoms}"
        )

    # Determine GRO molecule boundaries from topology order.
    cursor = 0
    molecule_ranges = {}

    for molecule, copies in molecule_order:

        if copies != 1 and molecule in ("system1", "MOL"):
            raise ValueError(
                f"Unexpected number of {molecule} copies: {copies}"
            )

        count = molecule_types[molecule]["atom_count"] * copies

        start = cursor + 1
        end = cursor + count

        molecule_ranges[molecule] = {
            "start": start,
            "end": end,
            "atom_count": count,
        }

        cursor += count

    if cursor != total_atoms:
        raise ValueError(
            "Molecule boundaries do not match GRO atom count."
        )

    protein_range = molecule_ranges["system1"]
    ligand_range = molecule_ranges["MOL"]

    # Validate residue labels within protein and ligand ranges.
    protein_atoms = atoms[
        protein_range["start"] - 1:
        protein_range["end"]
    ]

    ligand_atoms = atoms[
        ligand_range["start"] - 1:
        ligand_range["end"]
    ]

    if not protein_atoms:
        raise ValueError("No protein atoms found.")

    if not ligand_atoms:
        raise ValueError("No ligand atoms found.")

    if any(
        atom["residue_name"].upper() == "MOL"
        for atom in protein_atoms
    ):
        raise ValueError(
            "Ligand residue unexpectedly occurs in protein molecule."
        )

    if any(
        atom["residue_name"].upper() != "MOL"
        for atom in ligand_atoms
    ):
        raise ValueError(
            "Unexpected residue name in ligand molecule."
        )

    # Identify heavy atoms.
    protein_heavy = [
        atom for atom in protein_atoms
        if not is_hydrogen(atom["atom_name"])
    ]

    ligand_heavy = [
        atom for atom in ligand_atoms
        if not is_hydrogen(atom["atom_name"])
    ]

    if not protein_heavy or not ligand_heavy:
        raise ValueError("Could not identify heavy atoms.")

    # Convert global GRO indices to local molecule indices.
    protein_local_indices = [
        atom["global_index"] - protein_range["start"] + 1
        for atom in protein_heavy
    ]

    ligand_local_indices = [
        atom["global_index"] - ligand_range["start"] + 1
        for atom in ligand_heavy
    ]

    # Ensure all local indices are valid for the molecule type.
    if max(protein_local_indices) > expected_protein_atoms:
        raise ValueError("Protein restraint index exceeds molecule atom count.")

    if max(ligand_local_indices) > expected_ligand_atoms:
        raise ValueError("Ligand restraint index exceeds molecule atom count.")

    return {
        "topology": topology,
        "coordinates": coordinates,
        "original_nvt": original_nvt,
        "molecules": molecule_types,
        "molecule_order": molecule_order,
        "molecule_ranges": molecule_ranges,
        "total_atoms": total_atoms,
        "protein_atom_count": expected_protein_atoms,
        "ligand_atom_count": expected_ligand_atoms,
        "protein_heavy_count": len(protein_heavy),
        "ligand_heavy_count": len(ligand_heavy),
        "protein_local_indices": protein_local_indices,
        "ligand_local_indices": ligand_local_indices,
        "protein_heavy_atoms": protein_heavy,
        "ligand_heavy_atoms": ligand_heavy,
    }


# ------------------------------------------------------------
# POSITION RESTRAINT GENERATION
# ------------------------------------------------------------

def write_restraint_file(path, indices, force_constant, label):

    with open(path, "w") as handle:

        handle.write("; Generated by prepare_all_complex_nvt.py\n")
        handle.write(f"; {label}\n")
        handle.write("; Atom indices are local to the molecule type.\n")
        handle.write("; functype 1: harmonic position restraints.\n\n")

        handle.write("[ position_restraints ]\n")

        for index in indices:

            handle.write(
                f"{index:8d} "
                f"1 "
                f"{force_constant:.1f} "
                f"{force_constant:.1f} "
                f"{force_constant:.1f}\n"
            )


def insert_restraints(
    topology_path,
    protein_restraint_name,
    ligand_restraint_name,
    backup_path,
):
    """
    Insert restraint sections inside the correct molecule definitions.

    system1: protein
    MOL: ligand

    The original topology is backed up first.
    """

    original = topology_path.read_text(errors="replace")
    lines = original.splitlines()

    if "TNIK_NVT_POSITION_RESTRAINTS" in original:
        raise RuntimeError(
            "Restraint marker already exists. "
            "Refusing to insert duplicate restraints."
        )

    if backup_path.exists():
        raise RuntimeError(
            f"Backup already exists: {backup_path}. "
            "Refusing to overwrite it."
        )

    shutil.copy2(topology_path, backup_path)

    # Find each molecule type block.
    molecule_starts = []

    current_section = None

    for i, line in enumerate(lines):

        stripped = line.strip()

        if stripped.startswith("[") and stripped.endswith("]"):

            current_section = stripped.strip("[]").strip().lower()

            if current_section == "moleculetype":
                molecule_starts.append(i)

    if len(molecule_starts) < 4:
        raise ValueError(
            "Expected separate system1, MOL, Cl- and WAT molecule definitions."
        )

    molecule_names = {}

    for position, start in enumerate(molecule_starts):

        next_start = (
            molecule_starts[position + 1]
            if position + 1 < len(molecule_starts)
            else len(lines)
        )

        name = None

        for j in range(start + 1, next_start):

            stripped = lines[j].strip()

            if stripped and not stripped.startswith(";") and not stripped.startswith("#"):
                name = stripped.split()[0]
                break

        if name:
            molecule_names[name] = (start, next_start)

    if "system1" not in molecule_names or "MOL" not in molecule_names:
        raise ValueError(
            "Could not locate system1 and MOL molecule definitions."
        )

    # Determine insertion points immediately before the next molecule type.
    protein_start, protein_end = molecule_names["system1"]
    ligand_start, ligand_end = molecule_names["MOL"]

    protein_insert = protein_end
    ligand_insert = ligand_end

    additions = {
        protein_insert: [
            "",
            "; TNIK_NVT_POSITION_RESTRAINTS_PROTEIN",
            '#include "' + protein_restraint_name + '"',
            "",
        ],
        ligand_insert: [
            "",
            "; TNIK_NVT_POSITION_RESTRAINTS_LIGAND",
            '#include "' + ligand_restraint_name + '"',
            "",
        ],
    }

    new_lines = []

    for i, line in enumerate(lines):

        if i in additions:
            new_lines.extend(additions[i])

        new_lines.append(line)

    # Handle any insertion at end of file (not expected here).
    if len(lines) in additions:
        new_lines.extend(additions[len(lines)])

    topology_path.write_text("\n".join(new_lines) + "\n")


# ------------------------------------------------------------
# NVT MDP GENERATION
# ------------------------------------------------------------

def create_nvt_mdp(original_mdp, output_mdp):

    text = original_mdp.read_text(errors="replace")

    # Remove any previous define=-DPOSRES line.
    # Restraints are directly included in the topology.
    cleaned = []

    for line in text.splitlines():

        stripped = line.strip().lower()

        if stripped.startswith("define") and "posres" in stripped:
            continue

        cleaned.append(line)

    text = "\n".join(cleaned).strip() + "\n"

    # Replace NVT settings with controlled values.
    replacements = {
        "integrator": "md",
        "dt": str(TIMESTEP_PS),
        "nsteps": str(NVT_STEPS),
        "tcoupl": "v-rescale",
        "tc-grps": "System",
        "tau-t": "1.0",
        "ref-t": str(TEMPERATURE),
        "pcoupl": "no",
        "constraints": "h-bonds",
        "constraint_algorithm": "lincs",
        "nstenergy": "500",
        "nstlog": "500",
        "nstxout-compressed": "500",
        "nstxout": "0",
        "nstvout": "0",
        "nstfout": "0",
        "continuation": "no",
        "gen-vel": "yes",
        "gen-temp": str(TEMPERATURE),
        "gen-seed": "-1",
    }

    existing = {}

    for line in text.splitlines():

        stripped = line.strip()

        if not stripped or stripped.startswith(";"):
            continue

        if "=" in stripped:
            key = stripped.split("=", 1)[0].strip().lower()
            existing[key] = True

    output_lines = text.splitlines()

    for key, value in replacements.items():

        if key in existing:
            output_lines = [
                (
                    f"{key} = {value}"
                    if line.strip().split("=", 1)[0].strip().lower() == key
                    and "=" in line
                    else line
                )
                for line in output_lines
            ]
        else:
            output_lines.append(f"{key} = {value}")

    # Ensure no pressure coupling remains.
    output_mdp.write_text("\n".join(output_lines) + "\n")


# ------------------------------------------------------------
# GROMACS VALIDATION
# ------------------------------------------------------------

def validate_grompp(folder, nvt_mdp, system_name):

    if not GMX:
        raise RuntimeError(
            "GROMACS executable 'gmx' not found. "
            "Activate the tnik_md environment."
        )

    output_tpr = folder / "nvt_validation.tpr"

    command = [
        GMX,
        "grompp",
        "-f", str(nvt_mdp),
        "-c", "minimization.gro",
        "-r", "minimization.gro",
        "-p", "complex.top",
        "-o", output_tpr.name,
        "-maxwarn", "1",
    ]

    log_file = LOG_DIR / f"{system_name}_nvt_grompp.log"

    return_code = run_command(
        command,
        folder,
        log_file,
    )

    if return_code != 0 or not output_tpr.exists():
        raise RuntimeError(
            f"GROMACS NVT preprocessing failed. "
            f"See {log_file}"
        )

    # Inspect log for errors.
    log_text = log_file.read_text(errors="replace")

    if "Fatal error" in log_text:
        raise RuntimeError(
            f"Fatal GROMACS error found. See {log_file}"
        )

    return {
        "return_code": return_code,
        "tpr": str(output_tpr),
        "log": str(log_file),
        "maxwarn_used": 1,
    }


# ------------------------------------------------------------
# MAIN WORKFLOW
# ------------------------------------------------------------

def main():

    if not PROJECT.exists():
        raise FileNotFoundError(
            f"Project directory not found: {PROJECT}"
        )

    if not GMX:
        raise RuntimeError(
            "GROMACS 'gmx' was not found. "
            "Run: conda activate tnik_md"
        )

    report = {
        "created": datetime.now().isoformat(),
        "project": str(PROJECT),
        "gromacs": GMX,
        "nvt_duration_ps": NVT_DURATION_PS,
        "timestep_ps": TIMESTEP_PS,
        "nvt_steps": NVT_STEPS,
        "temperature_K": TEMPERATURE,
        "protein_restraint_kj_mol_nm2": PROTEIN_FORCE,
        "ligand_restraint_kj_mol_nm2": LIGAND_FORCE,
        "systems": {},
        "overall_success": False,
    }

    text_report = []

    log("\n" + "=" * 72)
    log("TNIK COMPLEX NVT PREPARATION")
    log("=" * 72)

    for system_name in SYSTEMS:

        log("\n" + "-" * 72)
        log(f"PROCESSING: {system_name}")
        log("-" * 72)

        folder = BASE / system_name

        result = {
            "success": False,
            "system": system_name,
        }

        try:

            if not folder.exists():
                raise FileNotFoundError(
                    f"System folder missing: {folder}"
                )

            # 1. Validate topology and coordinates.
            data = validate_system(folder)

            log(f"Total atoms: {data['total_atoms']}")
            log(
                f"Protein atoms: {data['protein_atom_count']}; "
                f"heavy atoms: {data['protein_heavy_count']}"
            )
            log(
                f"Ligand atoms: {data['ligand_atom_count']}; "
                f"heavy atoms: {data['ligand_heavy_count']}"
            )

            # 2. Generate restraint files.
            protein_restraint = folder / "posre_protein_heavy.itp"
            ligand_restraint = folder / "posre_ligand_heavy.itp"

            if protein_restraint.exists() or ligand_restraint.exists():
                raise RuntimeError(
                    "Restraint file already exists. "
                    "Please inspect it before rerunning."
                )

            write_restraint_file(
                protein_restraint,
                data["protein_local_indices"],
                PROTEIN_FORCE,
                "Protein heavy atom restraints",
            )

            write_restraint_file(
                ligand_restraint,
                data["ligand_local_indices"],
                LIGAND_FORCE,
                "Ligand heavy atom restraints",
            )

            log("Generated protein and ligand restraint files.")

            # 3. Back up and modify topology.
            backup = folder / "complex.top.pre_posres.bak"

            insert_restraints(
                data["topology"],
                protein_restraint.name,
                ligand_restraint.name,
                backup,
            )

            log("Created topology backup and inserted restraint includes.")

            # 4. Create a separate NVT MDP.
            nvt_mdp = folder / "nvt_posres.mdp"

            if nvt_mdp.exists():
                raise RuntimeError(
                    f"NVT output already exists: {nvt_mdp}"
                )

            create_nvt_mdp(
                data["original_nvt"],
                nvt_mdp,
            )

            log("Created nvt_posres.mdp.")

            # 5. Validate topology restraint placement.
            topology_text = data["topology"].read_text(errors="replace")

            if "posre_protein_heavy.itp" not in topology_text:
                raise RuntimeError("Protein restraint include missing.")

            if "posre_ligand_heavy.itp" not in topology_text:
                raise RuntimeError("Ligand restraint include missing.")

            # 6. Run GROMACS preprocessing.
            grompp_result = validate_grompp(
                folder,
                nvt_mdp,
                system_name,
            )

            log("GROMACS NVT preprocessing PASSED.")

            result.update({
                "success": True,
                "total_atoms": data["total_atoms"],
                "protein_atoms": data["protein_atom_count"],
                "protein_heavy_atoms": data["protein_heavy_count"],
                "ligand_atoms": data["ligand_atom_count"],
                "ligand_heavy_atoms": data["ligand_heavy_count"],
                "protein_restraint_file": str(protein_restraint),
                "ligand_restraint_file": str(ligand_restraint),
                "topology_backup": str(backup),
                "nvt_mdp": str(nvt_mdp),
                "grompp": grompp_result,
            })

            log(f"SUCCESS: {system_name}")

        except Exception as exc:

            result["error"] = str(exc)
            result["traceback"] = traceback.format_exc()

            log(f"FAILED: {system_name}")
            log(str(exc))

        report["systems"][system_name] = result

    report["overall_success"] = all(
        report["systems"].get(name, {}).get("success", False)
        for name in SYSTEMS
    )

    # 7. Write machine-readable report.
    REPORT_JSON.write_text(
        json.dumps(report, indent=2)
    )

    # 8. Write human-readable report.
    text_report.append("TNIK COMPLEX NVT PREPARATION REPORT")
    text_report.append("=" * 60)
    text_report.append(f"Date: {report['created']}")
    text_report.append(f"GROMACS: {GMX}")
    text_report.append(f"NVT duration: {NVT_DURATION_PS} ps")
    text_report.append(f"Time step: {TIMESTEP_PS} ps")
    text_report.append(f"NVT steps: {NVT_STEPS}")
    text_report.append("")

    for name in SYSTEMS:

        item = report["systems"].get(name, {})

        text_report.append(name)
        text_report.append("-" * 50)
        text_report.append(
            "Status: " + (
                "PASSED" if item.get("success") else "FAILED"
            )
        )

        if item.get("success"):

            text_report.append(
                f"Total atoms: {item['total_atoms']}"
            )
            text_report.append(
                f"Protein heavy atoms restrained: "
                f"{item['protein_heavy_atoms']}"
            )
            text_report.append(
                f"Ligand heavy atoms restrained: "
                f"{item['ligand_heavy_atoms']}"
            )
            text_report.append(
                f"NVT TPR: {item['grompp']['tpr']}"
            )

        else:
            text_report.append(
                "Error: " + str(item.get("error"))
            )

        text_report.append("")

    text_report.append(
        "OVERALL STATUS: " +
        ("ALL SYSTEMS PASSED" if report["overall_success"]
         else "ONE OR MORE SYSTEMS FAILED")
    )

    REPORT_TXT.write_text("\n".join(text_report) + "\n")

    log("\n" + "=" * 72)
    log("FINAL SUMMARY")
    log("=" * 72)

    for name in SYSTEMS:

        item = report["systems"][name]

        status = (
            "PASSED"
            if item.get("success")
            else "FAILED"
        )

        log(f"{name}: {status}")

    log("\nJSON report:")
    log(str(REPORT_JSON))

    log("\nText report:")
    log(str(REPORT_TXT))

    if report["overall_success"]:
        log("\nALL THREE SYSTEMS PASSED NVT PREPROCESSING.")
        log("No molecular dynamics simulation was started.")
        return 0

    log("\nSOME SYSTEMS FAILED. DO NOT START NVT.")
    log("Inspect the individual GROMACS logs and JSON report.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
