import pandas as pd

pk1 = pd.read_csv("pkCSM_results_batch_1.csv", sep="\t")
pk2 = pd.read_csv("pkCSM_results_batch_2.csv", sep="\t")
pk = pd.concat([pk1, pk2], ignore_index=True)

swiss = pd.read_csv("SwissADME_results.csv")

print("pkCSM records:", len(pk))
print("SwissADME records:", len(swiss))
print("pkCSM columns:", list(pk.columns))
print("SwissADME columns:", list(swiss.columns))
