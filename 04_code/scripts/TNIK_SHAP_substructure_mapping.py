
import os
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D

# Project directory
project_dir = os.path.expanduser("~/Desktop/TNIK_IPF_Project")

# Input files
activity_file = os.path.join(
    project_dir,
    "03_processed_data/activity_clean/TNIK_IC50_all_assays.csv"
)

importance_file = os.path.join(
    project_dir,
    "07_explainability/shap/TNIK_SHAP_feature_importance.csv"
)

# Output directories
output_dir = os.path.join(
    project_dir,
    "07_explainability/shap/substructure_mapping"
)

os.makedirs(output_dir, exist_ok=True)

# Load data
activity_df = pd.read_csv(activity_file)
importance_df = pd.read_csv(importance_file)

# Select top 10 SHAP fingerprint bits
top_bits = importance_df.head(10)["Fingerprint_Bit"].tolist()

# Prepare Morgan fingerprint generator
generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)

results = []

print("Mapping top 10 SHAP fingerprint bits...")

for bit_name in top_bits:

    bit_number = int(bit_name.replace("FP_", ""))

    matching_compounds = []

    for _, row in activity_df.iterrows():

        smiles = row.get("canonical_smiles")

        if pd.isna(smiles):
            continue

        mol = Chem.MolFromSmiles(str(smiles))

        if mol is None:
            continue

        additional_output = rdFingerprintGenerator.AdditionalOutput()
        additional_output.AllocateBitInfoMap()

        fingerprint = generator.GetFingerprint(
            mol,
            additionalOutput=additional_output
        )

        bit_info = additional_output.GetBitInfoMap()

        if bit_number not in bit_info:
            continue

        # Avoid duplicate compound records
        compound_id = row["molecule_chembl_id"]

        if compound_id in matching_compounds:
            continue

        matching_compounds.append(compound_id)

        for atom_index, radius in bit_info[bit_number]:

            atom_environment = Chem.FindAtomEnvironmentOfRadiusN(
                mol,
                radius,
                atom_index
            )

            atom_indices = {atom_index}

            for bond_index in atom_environment:
                bond = mol.GetBondWithIdx(bond_index)
                atom_indices.add(bond.GetBeginAtomIdx())
                atom_indices.add(bond.GetEndAtomIdx())

            substructure_smiles = Chem.MolFragmentToSmiles(
                mol,
                atomsToUse=sorted(atom_indices),
                bondsToUse=list(atom_environment),
                canonical=True
            )

            results.append({
                "Fingerprint_Bit": bit_name,
                "Compound_ID": compound_id,
                "pIC50": row.get("pIC50"),
                "Radius": radius,
                "Center_Atom": atom_index,
                "Substructure_SMILES": substructure_smiles,
                "Canonical_SMILES": Chem.MolToSmiles(mol)
            })

    print(
        bit_name,
        "- compounds containing bit:",
        len(matching_compounds)
    )

# Save mapping results
mapping_df = pd.DataFrame(results)

mapping_file = os.path.join(
    output_dir,
    "TNIK_SHAP_substructure_mapping.csv"
)

mapping_df.to_csv(mapping_file, index=False)

# Generate representative molecular images
print("\nGenerating molecular structure images...")

for bit_name in top_bits:

    bit_rows = mapping_df[
        mapping_df["Fingerprint_Bit"] == bit_name
    ]

    if bit_rows.empty:
        print("No matching molecule found for", bit_name)
        continue

    # Select the first matching compound as a representative example
    representative = bit_rows.iloc[0]

    mol = Chem.MolFromSmiles(
        representative["Canonical_SMILES"]
    )

    if mol is None:
        continue

    additional_output = rdFingerprintGenerator.AdditionalOutput()
    additional_output.AllocateBitInfoMap()

    generator.GetFingerprint(
        mol,
        additionalOutput=additional_output
    )

    bit_info = additional_output.GetBitInfoMap()
    bit_number = int(bit_name.replace("FP_", ""))

    if bit_number not in bit_info:
        continue

    atom_indices = set()
    bond_indices = set()

    for atom_index, radius in bit_info[bit_number]:

        atom_indices.add(atom_index)

        environment = Chem.FindAtomEnvironmentOfRadiusN(
            mol,
            radius,
            atom_index
        )

        for bond_index in environment:
            bond_indices.add(bond_index)

            bond = mol.GetBondWithIdx(bond_index)
            atom_indices.add(bond.GetBeginAtomIdx())
            atom_indices.add(bond.GetEndAtomIdx())

    drawer = rdMolDraw2D.MolDraw2DCairo(600, 450)

    drawer.DrawMolecule(
        mol,
        highlightAtoms=sorted(atom_indices),
        highlightBonds=sorted(bond_indices),
        legend=f"{bit_name} - {representative['Compound_ID']}"
    )

    drawer.FinishDrawing()

    image_file = os.path.join(
        output_dir,
        f"{bit_name}_representative.png"
    )

    with open(image_file, "wb") as file:
        file.write(drawer.GetDrawingText())

print("\nSubstructure mapping completed!")

print("Total mapping records:", len(mapping_df))

print("\nOutput folder:")
print(output_dir)

print("\nMapping CSV:")
print(mapping_file)
