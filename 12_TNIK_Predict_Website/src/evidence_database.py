from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]

EVIDENCE_FILE = ROOT / "data" / "compound_database" / "evidence.json"
PLIP_FILE = ROOT / "data" / "compound_database" / "plip_evidence.json"


def load_json(path):
    if not path.exists():
        return {}
    with open(path, "r") as file:
        return json.load(file)


def load_evidence():
    return load_json(EVIDENCE_FILE)


def load_plip_evidence():
    return load_json(PLIP_FILE)


def get_compound_evidence(compound_name):
    key = compound_name.strip().lower()
    evidence = load_evidence()
    plip = load_plip_evidence()

    return {
        "docking": evidence.get("docking", {}).get(key),
        "interactions": evidence.get("interactions", {}).get(key),
        "shap": evidence.get("shap", {}).get(key),
        "molecular_dynamics": evidence.get("molecular_dynamics", {}).get(key),
        "plip": plip.get(key)
    }
