from src.prediction import predict_compound

# Ethanol: a simple test molecule for checking the prediction pipeline.
result = predict_compound("CCO")

print("\n===== TNIK-PREDICT ENGINE TEST =====")

print("\nPredictions:")
for name, values in result["predictions"].items():
    print(
        name,
        "pIC50 =", round(values["predicted_pIC50"], 3),
        "estimated IC50 (nM) =",
        round(values["estimated_IC50_nM"], 2)
    )

print("\nMolecular properties:")
for key, value in result["properties"].items():
    print(key, ":", value)

print("\nClosest reference compound:", result["similarity"]["closest_compound"])
print("Maximum Tanimoto similarity:", round(result["similarity"]["maximum_similarity"], 3))

print("\nPipeline test completed.")
