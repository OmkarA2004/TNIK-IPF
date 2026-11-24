from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser

ROOT = Path(".")
PDB_FILE = ROOT / "06_docking/molecular_dynamics/structures/TNIK_6RA7_protein_raw.pdb"
BASE = ROOT / "06_docking/molecular_dynamics/energy_minimization"
OUT = ROOT / "06_docking/molecular_dynamics/analysis/npt_stability/residue_contacts"

SYSTEMS = {
    "Ponatinib_N32_H": BASE / "ponatinib_N32_H/npt_1ns/npt_1ns.gro",
    "Ponatinib_N37_H": BASE / "ponatinib_N37_H/npt_1ns/npt_1ns.gro",
    "Tepotinib_plus1": BASE / "tepotinib_pH74_plus1/npt_1ns/npt_1ns.gro",
}

AA_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "TPO",
}
NORMALIZE = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "CYX": "CYS", "CYM": "CYS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS",
}

def norm(name):
    return NORMALIZE.get(name.upper(), name.upper())

def original_residues():
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("TNIK", str(PDB_FILE))
    chain = structure[0]["A"]
    residues = []

    for residue in chain:
        name = residue.resname.strip().upper()
        if name in AA_NAMES:
            hetflag, number, insertion = residue.id
            residues.append({
                "pdb_chain": "A",
                "pdb_residue_number": int(number),
                "pdb_insertion_code": insertion.strip(),
                "pdb_residue_name": name,
            })

    return residues

def gro_protein_residues(path, expected_count):
    residues = []
    seen = set()

    with path.open() as f:
        next(f)  # title
        atom_count = int(next(f).strip())

        for line in f:
            if len(line) < 20:
                continue

            try:
                resnum = int(line[0:5])
                resname = line[5:10].strip().upper()
            except ValueError:
                continue

            key = (resnum, resname)
            if key in seen:
                continue
            seen.add(key)

            # Protein residues occur first in these prepared GRO files.
            if resname not in AA_NAMES and resname not in NORMALIZE:
                if residues:
                    break
                continue

            residues.append({
                "gromacs_residue_number": resnum,
                "gromacs_residue_name": resname,
            })

            if len(residues) == expected_count:
                break

    if len(residues) != expected_count:
        raise ValueError(
            f"{path}: expected {expected_count} protein residues, "
            f"found {len(residues)}"
        )

    return residues

def main():
    if not PDB_FILE.exists():
        raise FileNotFoundError(PDB_FILE)

    OUT.mkdir(parents=True, exist_ok=True)
    pdb_res = original_residues()

    print(f"Original PDB chain A protein residues: {len(pdb_res)}")
    if len(pdb_res) != 302:
        raise ValueError("Original PDB residue count differs from expected 302.")

    all_maps = {}

    for system, gro_path in SYSTEMS.items():
        if not gro_path.exists():
            raise FileNotFoundError(gro_path)

        gro_res = gro_protein_residues(gro_path, len(pdb_res))
        mismatches = []

        for i, (pdb, gro) in enumerate(zip(pdb_res, gro_res), start=1):
            if norm(pdb["pdb_residue_name"]) != norm(gro["gromacs_residue_name"]):
                mismatches.append((i, pdb, gro))

        print(f"\n{system}:")
        print(f"  GROMACS protein residues: {len(gro_res)}")
        print(f"  Sequence comparison: {'PASS' if not mismatches else 'FAIL'}")

        if mismatches:
            for i, pdb, gro in mismatches[:10]:
                print(
                    f"  Mismatch at sequence position {i}: "
                    f"PDB {pdb['pdb_residue_name']} "
                    f"{pdb['pdb_residue_number']}{pdb['pdb_insertion_code']} "
                    f"vs GRO {gro['gromacs_residue_name']} "
                    f"{gro['gromacs_residue_number']}"
                )
            raise ValueError(
                f"Residue sequence mismatch in {system}; "
                "no updated contact tables were written."
            )

        mapping = []
        for pdb, gro in zip(pdb_res, gro_res):
            mapping.append({
                **gro,
                **pdb,
            })

        all_maps[system] = pd.DataFrame(mapping)

    # Save the validated mapping for each system.
    for system, mapping in all_maps.items():
        mapping.to_csv(
            OUT / f"{system}_validated_residue_mapping.csv",
            index=False
        )

    # Read original contact occupancy results.
    source = OUT / "TNIK_all_complexes_residue_contact_occupancy.csv"
    contacts = pd.read_csv(source)

    updated = []
    for system, mapping in all_maps.items():
        subset = contacts[contacts["system"] == system].copy()
        subset["residue_number"] = pd.to_numeric(
            subset["residue_number"], errors="raise"
        ).astype(int)

        subset = subset.merge(
            mapping[
                [
                    "gromacs_residue_number",
                    "pdb_chain",
                    "pdb_residue_number",
                    "pdb_insertion_code",
                    "pdb_residue_name",
                ]
            ],
            left_on="residue_number",
            right_on="gromacs_residue_number",
            how="left",
            validate="many_to_one",
        )

        if subset["pdb_residue_number"].isna().any():
            raise ValueError(f"Unmapped contact residue found in {system}.")

        subset.drop(columns=["gromacs_residue_number"], inplace=True)
        updated.append(subset)

    result = pd.concat(updated, ignore_index=True)
    result.to_csv(
        OUT / "TNIK_all_complexes_residue_contact_occupancy_VALIDATED_PDB.csv",
        index=False
    )

    pivot = result.pivot_table(
        index=[
            "residue_number",
            "pdb_chain",
            "pdb_residue_number",
            "pdb_insertion_code",
            "pdb_residue_name",
        ],
        columns="system",
        values="occupancy_percent",
        fill_value=0,
    ).reset_index()

    pivot.to_csv(
        OUT / "TNIK_residue_contact_comparison_VALIDATED_PDB.csv",
        index=False
    )

    print("\n" + "=" * 60)
    print("VALIDATION SUCCESSFUL")
    print("All three systems match the original PDB chain A sequence.")
    print("Validated mapping files and contact tables saved in:")
    print(OUT)
    print("\nSample validated residue mapping:")
    print(
        all_maps["Ponatinib_N32_H"][
            [
                "gromacs_residue_number",
                "gromacs_residue_name",
                "pdb_chain",
                "pdb_residue_number",
                "pdb_insertion_code",
                "pdb_residue_name",
            ]
        ].head(12).to_string(index=False)
    )
    print("\nUpdated comparison preview:")
    print(pivot.head(15).to_string(index=False))

if __name__ == "__main__":
    main()
