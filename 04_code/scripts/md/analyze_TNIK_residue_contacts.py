from pathlib import Path
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("06_docking/molecular_dynamics")
ANALYSIS = BASE / "analysis/npt_stability"
INDEX_DIR = ANALYSIS / "index"
TRAJ_DIR = ANALYSIS / "processed_trajectories"
OUTPUT = ANALYSIS / "residue_contacts"

OUTPUT.mkdir(parents=True, exist_ok=True)

CUTOFF_NM = 0.35
EXPECTED_FRAMES = 101

SYSTEMS = {
    "Ponatinib_N32_H": {
        "tpr": BASE / "energy_minimization/ponatinib_N32_H/npt_1ns/npt_1ns.tpr",
        "gro": BASE / "energy_minimization/ponatinib_N32_H/npt_1ns/npt_1ns.gro",
        "xtc": TRAJ_DIR / "ponatinib_N32_H_whole_centered.xtc",
        "ndx": INDEX_DIR / "ponatinib_N32_H_contact_groups.ndx",
    },
    "Ponatinib_N37_H": {
        "tpr": BASE / "energy_minimization/ponatinib_N37_H/npt_1ns/npt_1ns.tpr",
        "gro": BASE / "energy_minimization/ponatinib_N37_H/npt_1ns/npt_1ns.gro",
        "xtc": TRAJ_DIR / "ponatinib_N37_H_whole_centered.xtc",
        "ndx": INDEX_DIR / "ponatinib_N37_H_contact_groups.ndx",
    },
    "Tepotinib_plus1": {
        "tpr": BASE / "energy_minimization/tepotinib_pH74_plus1/npt_1ns/npt_1ns.tpr",
        "gro": BASE / "energy_minimization/tepotinib_pH74_plus1/npt_1ns/npt_1ns.gro",
        "xtc": TRAJ_DIR / "tepotinib_pH74_plus1_whole_centered.xtc",
        "ndx": INDEX_DIR / "tepotinib_pH74_plus1_contact_groups.ndx",
    },
}

SOLVENT_IONS = {"SOL", "WAT", "HOH", "CL", "NA", "K", "MG", "CA"}


def run_command(command, logfile):
    print("\nRunning:", " ".join(map(str, command)))
    with open(logfile, "w") as log:
        result = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True
        )

    if result.returncode != 0:
        print(f"\nGROMACS error. See: {logfile}")
        print("\nLast 20 log lines:")
        try:
            print("\n".join(logfile.read_text().splitlines()[-20:]))
        except Exception:
            pass
        raise RuntimeError(f"Command failed: {command[0]}")


def read_protein_atom_map(gro):
    """Map GRO atom indices to protein residue identities."""
    atom_map = {}

    with gro.open() as f:
        lines = f.readlines()

    for line in lines[2:-1]:
        try:
            index = int(line[15:20])
            resid = int(line[0:5])
            resname = line[5:10].strip()
            atomname = line[10:15].strip()

            if resname not in SOLVENT_IONS and resname != "MOL":
                atom_map[index] = {
                    "residue_number": resid,
                    "residue_name": resname,
                    "atom_name": atomname,
                }
        except (ValueError, IndexError):
            continue

    return atom_map


def parse_gromacs_oi_dat(path, atom_map):
    """
    Parse gmx select -oi .dat output.
    Each data row is expected to contain a time followed by
    the selected atom indices for that frame.
    """
    frames = []

    for line in path.read_text().splitlines():
        line = line.strip()

        if not line or line.startswith(("#", "@")):
            continue

        tokens = line.split()

        try:
            time_ps = float(tokens[0])
        except (ValueError, IndexError):
            continue

        selected = set()

        for token in tokens[1:]:
            try:
                atom_index = int(token)
            except ValueError:
                continue

            if atom_index in atom_map:
                selected.add(atom_index)

        frames.append({
            "time_ps": time_ps,
            "protein_atoms": selected
        })

    return frames


def analyse_system(system, files):
    print("\n" + "=" * 65)
    print("ANALYSING:", system)
    print("=" * 65)

    for key, path in files.items():
        if not path.exists():
            raise FileNotFoundError(f"{system}: missing {key}: {path}")

    atom_map = read_protein_atom_map(files["gro"])
    print("Protein atoms mapped:", len(atom_map))

    raw_output = OUTPUT / f"{system}_protein_contact_atoms_per_frame.dat"
    log_file = OUTPUT / f"{system}_gmx_select.log"

    selection = (
        f'group "Protein" and '
        f'within {CUTOFF_NM} of group "resname_MOL"'
    )

    command = [
        "gmx", "select",
        "-s", str(files["tpr"]),
        "-f", str(files["xtc"]),
        "-n", str(files["ndx"]),
        "-select", selection,
        "-oi", str(raw_output),
        "-selrpos", "atom",
        "-xvg", "none",
    ]

    run_command(command, log_file)

    if not raw_output.exists():
        raise RuntimeError(f"Selection output was not created: {raw_output}")

    frames = parse_gromacs_oi_dat(raw_output, atom_map)

    if len(frames) != EXPECTED_FRAMES:
        print("Parsed frame count:", len(frames))
        print("First five parsed times:", [x["time_ps"] for x in frames[:5]])
        raise RuntimeError(
            f"{system}: expected {EXPECTED_FRAMES} frames, "
            f"but parsed {len(frames)}. Inspect {raw_output}"
        )

    residue_frames = {}

    for frame_number, frame in enumerate(frames):
        present_residues = set()

        for atom_index in frame["protein_atoms"]:
            atom = atom_map[atom_index]
            residue = (
                atom["residue_number"],
                atom["residue_name"]
            )
            present_residues.add(residue)

        for residue in present_residues:
            residue_frames.setdefault(residue, set()).add(frame_number)

    rows = []

    for (resid, resname), occurrence_frames in residue_frames.items():
        count = len(occurrence_frames)

        rows.append({
            "system": system,
            "residue_number": resid,
            "residue_name": resname,
            "frames_in_contact": count,
            "total_frames": len(frames),
            "occupancy_percent": 100 * count / len(frames),
        })

    columns = [
        "system", "residue_number", "residue_name",
        "frames_in_contact", "total_frames", "occupancy_percent"
    ]
    df = pd.DataFrame(rows, columns=columns)

    if not df.empty:
        df = df.sort_values(
            ["occupancy_percent", "residue_number"],
            ascending=[False, True]
        )

    outfile = OUTPUT / f"{system}_residue_contact_occupancy.csv"
    df.to_csv(outfile, index=False)

    print("Frames analysed:", len(frames))
    print("Frames with any protein contact:",
          sum(bool(x["protein_atoms"]) for x in frames))
    print("Residues contacting ligand:", len(df))
    print("Saved:", outfile)

    return df


def main():
    all_results = []

    for system, files in SYSTEMS.items():
        all_results.append(analyse_system(system, files))

    combined = pd.concat(all_results, ignore_index=True)

    combined_file = OUTPUT / "TNIK_all_complexes_residue_contact_occupancy.csv"
    combined.to_csv(combined_file, index=False)

    comparison = combined.pivot_table(
        index=["residue_number", "residue_name"],
        columns="system",
        values="occupancy_percent",
        fill_value=0
    ).reset_index()

    comparison_file = OUTPUT / "TNIK_residue_contact_comparison.csv"
    comparison.to_csv(comparison_file, index=False)

    # Individual top-residue plots.
    for system in SYSTEMS:
        subset = combined[combined["system"] == system].nlargest(
            15, "occupancy_percent"
        )

        if subset.empty:
            continue

        labels = (
            subset["residue_name"].astype(str)
            + subset["residue_number"].astype(str)
        )

        plt.figure(figsize=(10, 6))
        plt.barh(labels, subset["occupancy_percent"])
        plt.xlabel("Trajectory frames with residue contact (%)")
        plt.ylabel("TNIK residue")
        plt.title(f"{system}: residue contact occupancy")
        plt.xlim(0, 100)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(
            OUTPUT / f"{system}_top_residue_contacts.png",
            dpi=300
        )
        plt.close()

    # Combined comparison plot.
    if not combined.empty:
        top_residues = (
            combined.groupby(["residue_number", "residue_name"])
            ["occupancy_percent"].max()
            .nlargest(15).index
        )

        plot_data = combined[
            combined.apply(
                lambda row: (
                    row["residue_number"],
                    row["residue_name"]
                ) in top_residues,
                axis=1
            )
        ]

        if not plot_data.empty:
            pivot = plot_data.pivot_table(
                index=["residue_name", "residue_number"],
                columns="system",
                values="occupancy_percent",
                fill_value=0
            )

            pivot.index = [
                f"{name}{number}" for name, number in pivot.index
            ]

            ax = pivot.plot(kind="bar", figsize=(14, 7), width=0.8)
            ax.set_ylabel("Contact occupancy (%)")
            ax.set_xlabel("TNIK residue")
            ax.set_title("TNIK residue contact occupancy across complexes")
            ax.set_ylim(0, 100)
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            plt.savefig(
                OUTPUT / "TNIK_three_complexes_residue_contact_comparison.png",
                dpi=300
            )
            plt.close()

    summary_file = OUTPUT / "TNIK_residue_contact_summary.txt"

    with summary_file.open("w") as f:
        f.write("TNIK–LIGAND RESIDUE CONTACT ANALYSIS\n")
        f.write("=" * 55 + "\n")
        f.write(f"Contact cutoff: {CUTOFF_NM} nm\n")
        f.write(f"Frames per trajectory: {EXPECTED_FRAMES}\n")

        for system in SYSTEMS:
            subset = combined[combined["system"] == system]
            f.write(f"\n\n{system}\n")
            f.write("-" * 45 + "\n")

            if subset.empty:
                f.write("No contacting residues detected.\n")
            else:
                f.write(subset.head(15).to_string(index=False) + "\n")

        f.write(
            "\nInterpretation: Residue contact occupancy measures "
            "geometric proximity, not binding affinity or proof of "
            "biological inhibition. Protein atoms are selected by "
            "GROMACS within the specified distance of the ligand group.\n"
        )

    print("\n" + "=" * 65)
    print("ANALYSIS COMPLETED")
    print("=" * 65)
    print("Combined CSV:", combined_file)
    print("Comparison CSV:", comparison_file)
    print("Summary:", summary_file)
    print("Output directory:", OUTPUT)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("\nANALYSIS STOPPED:", exc)
        sys.exit(1)
