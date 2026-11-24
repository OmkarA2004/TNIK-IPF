from src.compound_intelligence import identify_compound

tests = [
    "Aspirin",
    "CHEMBL25",
    "CCO"
]

for query in tests:
    print("\nTesting:", query)

    try:
        result = identify_compound(query)

        print("Name:", result["compound_name"])
        print("SMILES:", result["canonical_smiles"])
        print("Formula:", result["molecular_formula"])
        print("Molecular weight:", result["molecular_weight"])
        print("Source:", result["identification_source"])

    except Exception as error:
        print("Identification issue:", error)
