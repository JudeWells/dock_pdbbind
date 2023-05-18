import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("rmsd.csv")
df = df[df.rmsd != -1]

dists = [1,2,3,4,5,6]
proportions = []
for dist in dists:
    best_pose_within = df[(df["rmsd"] <= dist) & (df.vina_rank==0)]
    proportion = len(best_pose_within)/len(df[df.vina_rank==0])
    proportions.append(proportion)
    print(f"Proportion of ligands with best pose within {dist}A: {proportion}")
plt.plot(dists, proportions)
plt.show()

plt.clf()
plt.hist(df.rmsd.values, bins=100)
plt.show()
bp=1

