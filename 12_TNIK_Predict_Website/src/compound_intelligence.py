import re
import time
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

HEADERS = {
    "User-Agent": "TNIK-Predict-Research-Platform/1.0"
}


def get_response(url, timeout=25):
    """Retry temporary server and network failures."""
    last_error = None

    for attempt in range(3):
        try:
            response = requests.get(
                url,
                headers=HEADERS,
                timeout=timeout
            )

            if response.status_code in (429, 500, 502, 503, 504):
                last_error = RuntimeError(
                    f"Temporary database error: HTTP {response.status_code}"
                )
                time.sleep(2 ** attempt)
                continue

            response.raise_for_status()
            return response

        except requests.RequestException as error:
            last_error = error
            time.sleep(2 ** attempt)

    raise RuntimeError(
        f"Could not retrieve compound information: {last_error}"
    )


def parse_structure(smiles):
    """Reject missing or invalid chemical structures."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("The database returned an empty SMILES structure.")

    molecule = Chem.MolFromSmiles(smiles.strip())

    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError("The chemical structure is invalid or empty.")

    return molecule


def identify_compound(query):
    if not isinstance(query, str) or not query.strip():
        raise ValueError("Enter a compound name, identifier or SMILES.")

    query = query.strip()

    if re.fullmatch(r"CHEMBL\d+", query, re.IGNORECASE):
        return identify_chembl(query.upper())

    if query.isdigit():
        return identify_pubchem_cid(query)

    molecule = Chem.MolFromSmiles(query)

    if molecule is not None and molecule.GetNumAtoms() > 0:
        return build_compound_profile(
            molecule,
            input_value=query,
            identification_source="User-provided SMILES"
        )

    return identify_pubchem_name(query)


def identify_chembl(chembl_id):
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/molecule/"
        f"{chembl_id}.json"
    )

    response = get_response(url)
    data = response.json()

    structures = data.get("molecule_structures") or {}
    smiles = structures.get("canonical_smiles")

    molecule = parse_structure(smiles)

    return build_compound_profile(
        molecule,
        input_value=chembl_id,
        name=data.get("pref_name") or chembl_id,
        chembl_id=chembl_id,
        identification_source="ChEMBL"
    )


def identify_pubchem_cid(cid):
    properties = (
        "Title,CanonicalSMILES,IsomericSMILES,"
        "MolecularFormula,InChIKey"
    )

    url = (
        f"{PUBCHEM_URL}/compound/cid/{cid}/property/"
        f"{properties}/JSON"
    )

    response = get_response(url)
    records = (
        response.json()
        .get("PropertyTable", {})
        .get("Properties", [])
    )

    if not records:
        raise ValueError("PubChem returned no compound record.")

    record = records[0]
    smiles = (
        record.get("IsomericSMILES")
        or record.get("CanonicalSMILES")
    )

    molecule = parse_structure(smiles)

    return build_compound_profile(
        molecule,
        input_value=cid,
        name=record.get("Title") or f"PubChem CID {cid}",
        pubchem_cid=str(cid),
        identification_source="PubChem"
    )


def identify_pubchem_name(name):
    url = (
        f"{PUBCHEM_URL}/compound/name/"
        f"{requests.utils.quote(name, safe='')}/cids/JSON"
    )

    response = get_response(url)
    cids = (
        response.json()
        .get("IdentifierList", {})
        .get("CID", [])
    )

    if not cids:
        raise ValueError("No matching PubChem compound was found.")

    if len(cids) > 1:
        raise ValueError(
            "Multiple PubChem records match this name. "
            "Please provide a PubChem CID, ChEMBL ID or SMILES."
        )

    return identify_pubchem_cid(str(cids[0]))


def build_compound_profile(
    molecule,
    input_value,
    name=None,
    chembl_id=None,
    pubchem_cid=None,
    identification_source="User-provided SMILES"
):
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError("Cannot create a profile without a valid structure.")

    return {
        "input": input_value,
        "compound_name": name or "Name not established",
        "chembl_id": chembl_id,
        "pubchem_cid": pubchem_cid,
        "canonical_smiles": Chem.MolToSmiles(molecule),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight": round(Descriptors.MolWt(molecule), 2),
        "logp": round(Descriptors.MolLogP(molecule), 2),
        "tpsa": round(Descriptors.TPSA(molecule), 2),
        "hbd": rdMolDescriptors.CalcNumHBD(molecule),
        "hba": rdMolDescriptors.CalcNumHBA(molecule),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(molecule),
        "identification_source": identification_source,
        "molecule": molecule
    }
