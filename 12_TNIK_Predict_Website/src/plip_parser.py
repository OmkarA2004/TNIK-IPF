from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
PLIP_DIRECTORY = ROOT / "data" / "docking_evidence"


def get_text(element, tag, default=None):
    child = element.find(tag)
    return child.text.strip() if child is not None and child.text else default


def parse_plip_report(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()

    report = {
        "ligand": None,
        "ligand_smiles": None,
        "ligand_inchikey": None,
        "ligand_properties": {},
        "binding_site_residues": [],
        "interactions": {}
    }

    binding_site = root.find(".//bindingsite")

    if binding_site is None:
        return report

    identifiers = binding_site.find("identifiers")

    if identifiers is not None:
        report["ligand"] = get_text(identifiers, "longname")
        report["ligand_smiles"] = get_text(identifiers, "smiles")
        report["ligand_inchikey"] = get_text(identifiers, "inchikey")

    properties = binding_site.find("lig_properties")

    if properties is not None:
        for item in properties:
            try:
                value = float(item.text)
                if value.is_integer():
                    value = int(value)
            except (ValueError, TypeError):
                value = item.text

            report["ligand_properties"][item.tag] = value

    residues = binding_site.find("bs_residues")

    if residues is not None:
        for residue in residues.findall("bs_residue"):
            report["binding_site_residues"].append({
                "residue": residue.get("aa"),
                "position": residue.text.strip() if residue.text else None,
                "contact": residue.get("contact") == "True",
                "minimum_distance_angstrom": float(
                    residue.get("min_dist")
                ) if residue.get("min_dist") else None
            })

    # Extract interaction categories and their reported details.
    interaction_container = binding_site.find("interactions")

    if interaction_container is not None:
        for category in interaction_container:
            category_name = category.tag
            report["interactions"][category_name] = []

            for interaction in category:
                details = {}

                for field in interaction:
                    details[field.tag] = (
                        field.text.strip()
                        if field.text else None
                    )

                report["interactions"][category_name].append(details)

    return report


def load_all_plip_reports():
    results = {}

    for xml_file in PLIP_DIRECTORY.glob("*.xml"):
        compound = xml_file.stem.replace("_seed11", "").replace(
            "_seed22", ""
        ).lower()

        results[compound] = parse_plip_report(xml_file)

    return results


if __name__ == "__main__":
    reports = load_all_plip_reports()

    for compound, report in reports.items():
        print("\nCOMPOUND:", compound)
        print("Ligand SMILES:", report["ligand_smiles"])
        print("Residues:", len(report["binding_site_residues"]))

        for interaction_type, records in report["interactions"].items():
            print(interaction_type, ":", len(records))


def save_plip_evidence():
    import json

    output = ROOT / "data" / "compound_database" / "plip_evidence.json"
    reports = load_all_plip_reports()

    with open(output, "w") as file:
        json.dump(reports, file, indent=2)

    print(f"\nSaved PLIP evidence for {len(reports)} compounds.")
    print(f"Output: {output}")


if __name__ == "__main__":
    save_plip_evidence()
