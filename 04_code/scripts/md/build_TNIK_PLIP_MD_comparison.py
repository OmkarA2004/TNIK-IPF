from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd

ROOT = Path.cwd()
PLIP_DIR = ROOT / "06_docking/interaction_analysis/PLIP_results"
MD_DIR = ROOT / "06_docking/molecular_dynamics/analysis/npt_stability/residue_contacts"
OUT = ROOT / "06_docking/interaction_analysis/comparative_analysis"

OUT.mkdir(parents=True, exist_ok=True)

SYSTEMS = {
    "ponatinib_seed11": "Ponatinib",
    "tepotinib_seed22": "Tepotinib",
}

# PLIP XML interaction element names
INTERACTION_TYPES = {
    "hydrophobic_interaction": "Hydrophobic",
    "hydrogen_bond": "Hydrogen bond",
    "salt_bridge": "Salt bridge",
    "pi_stack": "Pi stacking",
    "pi_cation_interaction": "Pi-cation",
    "halogen_bond": "Halogen bond",
    "water_bridge": "Water bridge",
    "metal_complex": "Metal complex",
}

def read_plip_xml(system, ligand):
    xml_path = PLIP_DIR / system / f"{system}.xml"

    if not xml_path.exists():
        raise FileNotFoundError(f"Missing PLIP XML: {xml_path}")

    root = ET.parse(xml_path).getroot()
    records = []

    for element in root.iter():
        tag = element.tag.split("}")[-1].lower()

        if tag not in INTERACTION_TYPES:
            continue

        data = dict(element.attrib)

        # Some PLIP XML versions store fields as child elements.
        for child in element:
            child_tag = child.tag.split("}")[-1]
            if child.text and child.text.strip():
                data.setdefault(child_tag, child.text.strip())

        residue_number = data.get("resnr")
        residue_name = data.get("restype")

        if residue_number is None or residue_name is None:
            continue

        records.append({
            "ligand": ligand,
            "plip_system": system,
            "interaction_type": INTERACTION_TYPES[tag],
            "pdb_chain": data.get("reschain", ""),
            "pdb_residue_number": int(float(residue_number)),
            "residue_name": residue_name,
            "distance_A": data.get(
                "dist",
                data.get("dist_d-a", "")
            ),
            "hydrogen_acceptor_distance_A": data.get("dist_h-a", ""),
            "donor_angle_deg": data.get("don_angle", ""),
            "protein_is_donor": data.get("protisdon", ""),
            "sidechain": data.get("sidechain", ""),
            "ligand_group": data.get("lig_group", ""),
            "raw_attributes": str(data),
        })

    return records

print("=" * 72)
print("TNIK PLIP–MD COMPARATIVE ANALYSIS")
print("=" * 72)

# 1. Extract all PLIP interactions.
all_records = []

for system, ligand in SYSTEMS.items():
    records = read_plip_xml(system, ligand)
    all_records.extend(records)

    print(f"{ligand}: {len(records)} interaction records extracted")

if not all_records:
    raise RuntimeError("No PLIP interactions found. Check XML structure.")

details = pd.DataFrame(all_records)
details = details.sort_values(
    ["ligand", "interaction_type", "pdb_residue_number"]
)
details_path = OUT / "TNIK_PLIP_interaction_details.csv"
details.to_csv(details_path, index=False)

# 2. Interaction summary.
summary = (
    details.groupby(["ligand", "interaction_type"])
    .agg(
        interaction_records=("interaction_type", "size"),
        unique_residues=("pdb_residue_number", "nunique"),
    )
    .reset_index()
)

summary_path = OUT / "TNIK_PLIP_interaction_summary.csv"
summary.to_csv(summary_path, index=False)

# 3. Read validated MD residue-contact comparison.
md_path = MD_DIR / "TNIK_residue_contact_comparison_VALIDATED_PDB.csv"

if not md_path.exists():
    raise FileNotFoundError(
        f"Validated MD comparison not found: {md_path}"
    )

md = pd.read_csv(md_path)

required = {
    "pdb_residue_number",
    "pdb_residue_name",
    "Ponatinib_N32_H",
    "Ponatinib_N37_H",
    "Tepotinib_plus1",
}
missing = required - set(md.columns)

if missing:
    raise ValueError(
        f"Expected MD columns missing: {sorted(missing)}. "
        f"Found: {list(md.columns)}"
    )

md_clean = md[
    [
        "pdb_residue_number",
        "pdb_residue_name",
        "Ponatinib_N32_H",
        "Ponatinib_N37_H",
        "Tepotinib_plus1",
    ]
].copy()

md_clean["pdb_residue_number"] = pd.to_numeric(
    md_clean["pdb_residue_number"], errors="raise"
).astype(int)

md_clean["md_residue_name"] = (
    md_clean["pdb_residue_name"].astype(str).str.upper()
)

md_clean = md_clean.rename(columns={
    "Ponatinib_N32_H": "Ponatinib_N32_H_contact_occupancy_percent",
    "Ponatinib_N37_H": "Ponatinib_N37_H_contact_occupancy_percent",
    "Tepotinib_plus1": "Tepotinib_plus1_contact_occupancy_percent",
})

md_clean = md_clean.drop(columns=["pdb_residue_name"])

if md_clean["pdb_residue_number"].duplicated().any():
    raise ValueError("Duplicate PDB residue numbers found in MD comparison.")

# Keep each docking ligand separate, and attach the matching MD measurements.
comparison = details.merge(
    md_clean,
    on="pdb_residue_number",
    how="left"
)

comparison["residue_name_match"] = (
    comparison["residue_name"].str.upper()
    == comparison["md_residue_name"]
)

comparison_path = OUT / "TNIK_PLIP_MD_comparison.csv"
comparison.to_csv(comparison_path, index=False)

# 4. Human-readable summary.
summary_text = []
summary_text.append("TNIK PLIP–MD COMPARATIVE ANALYSIS")
summary_text.append("=" * 60)
summary_text.append("")
summary_text.append(
    "PLIP interaction records are from representative docking poses."
)
summary_text.append(
    "MD values are residue-contact occupancies from validated PDB numbering."
)
summary_text.append(
    "Contact occupancy is not equivalent to interaction occupancy."
)
summary_text.append("")

for ligand in details["ligand"].unique():
    subset = details[details["ligand"] == ligand]
    summary_text.append(f"{ligand.upper()}")
    summary_text.append("-" * 40)

    counts = subset["interaction_type"].value_counts()
    for interaction, count in counts.items():
        summary_text.append(f"{interaction}: {count} records")

    residues = (
        subset[["residue_name", "pdb_residue_number"]]
        .drop_duplicates()
        .sort_values("pdb_residue_number")
    )
    labels = [
        f"{r.residue_name}{r.pdb_residue_number}"
        for r in residues.itertuples()
    ]
    summary_text.append("Interacting residues: " + ", ".join(labels))
    summary_text.append("")

summary_text.append("MD residue-number mapping checks")
summary_text.append("-" * 40)
summary_text.append(
    f"Rows with a matching MD residue name: "
    f"{int(comparison['residue_name_match'].sum())}/{len(comparison)}"
)
summary_text.append(
    f"Rows without an MD occupancy match: "
    f"{int(comparison['Ponatinib_N32_H_contact_occupancy_percent'].isna().sum())}"
    " (Ponatinib N32 column)"
)
summary_text.append(
    "Inspect residue_name_match=False rows before interpreting comparisons."
)

summary_path_txt = OUT / "TNIK_PLIP_MD_comparison_summary.txt"
summary_path_txt.write_text("\n".join(summary_text) + "\n")

# 5. Print concise results.
print("\nInteraction counts:")
print(summary.to_string(index=False))

print("\nResidue mapping validation:")
print(
    comparison[
        ["ligand", "residue_name", "pdb_residue_number", "residue_name_match"]
    ].drop_duplicates().to_string(index=False)
)

print("\nOutput files:")
for path in [
    details_path,
    summary_path,
    comparison_path,
    summary_path_txt,
]:
    print(f"{path.relative_to(ROOT)} ({path.stat().st_size:,} bytes)")

print("\nAnalysis completed.")
