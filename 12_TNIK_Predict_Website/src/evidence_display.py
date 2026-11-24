import streamlit as st
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

IMAGE_FILES = {
    "ponatinib": [
        ROOT / "assets/docking_images/TNIK_ponatinib_PLIP_closeup.png",
        ROOT / "assets/docking_images/TNIK_Ponatinib_PLIP.png",
    ],
    "tepotinib": [
        ROOT / "assets/docking_images/TNIK_tepotinib_PLIP_closeup.png",
        ROOT / "assets/docking_images/TNIK_Tepotinib_PLIP.png",
    ],
}


def display_plip_evidence(compound_name, plip):
    if not plip:
        st.info("No PLIP interaction report is available for this compound.")
        return

    st.subheader("Protein–ligand interaction analysis")

    residues = plip.get("binding_site_residues", [])
    interactions = plip.get("interactions", {})

    st.metric("Binding-site residues", len(residues))

    st.markdown("#### Interaction summary")

    labels = {
        "hydrophobic_interactions": "Hydrophobic",
        "hydrogen_bonds": "Hydrogen bonds",
        "water_bridges": "Water bridges",
        "salt_bridges": "Salt bridges",
        "pi_stacks": "Pi stacking",
        "pi_cation_interactions": "Pi-cation",
        "halogen_bonds": "Halogen bonds",
        "metal_complexes": "Metal complexes",
    }

    columns = st.columns(4)

    for index, (key, label) in enumerate(labels.items()):
        count = len(interactions.get(key, []))
        columns[index % 4].metric(label, count)

    st.markdown("#### Binding-site residue details")

    if residues:
        st.dataframe(residues, use_container_width=True)

    key = compound_name.strip().lower()

    st.markdown("#### Molecular interaction images")

    for image_path in IMAGE_FILES.get(key, []):
        if image_path.exists():
            st.image(
                str(image_path),
                caption=image_path.name,
                use_container_width=True
            )
