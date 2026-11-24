import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

FINGERPRINT_RADIUS = 2
FINGERPRINT_SIZE = 2048

generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=FINGERPRINT_RADIUS,
    fpSize=FINGERPRINT_SIZE
)


def smiles_to_fingerprint(smiles):
    """Convert a valid SMILES string into a 2048-bit Morgan fingerprint."""

    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Please provide a valid SMILES string.")

    molecule = Chem.MolFromSmiles(smiles.strip())

    if molecule is None:
        raise ValueError("Invalid SMILES. Please check the chemical structure.")

    fingerprint = generator.GetFingerprint(molecule)

    features = np.zeros((FINGERPRINT_SIZE,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fingerprint, features)

    return molecule, features


def molecule_properties(molecule):
    """Calculate basic molecular properties."""

    from rdkit.Chem import Descriptors, rdMolDescriptors

    return {
        "Molecular weight": round(Descriptors.MolWt(molecule), 2),
        "LogP": round(Descriptors.MolLogP(molecule), 2),
        "TPSA": round(Descriptors.TPSA(molecule), 2),
        "H-bond donors": rdMolDescriptors.CalcNumHBD(molecule),
        "H-bond acceptors": rdMolDescriptors.CalcNumHBA(molecule),
        "Rotatable bonds": rdMolDescriptors.CalcNumRotatableBonds(molecule),
    }
