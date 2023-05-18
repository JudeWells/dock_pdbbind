import numpy as np
import pandas as pd
"""
Created by Jude Wells 2023-05-15
split the features into train, val and test sets
based on the CATH Topology split that was created in 
train_test_split_by_cath_code.py
"""

df = pd.read_csv("pdb_bind.csv")
df["mol_id"] = df.pdb_id + df.vina_rank.apply(lambda x:"_ligand_" + str(x + 1))
feature_paths = [
"/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/ECIFs/pdbbind_docked_2",
"/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/ECIFs/pdbbind_docked_1"
]
features = []
feature_mol_ids = []
for feature_path in feature_paths:
    dataset = pd.read_csv(feature_path, index_col=[0, 1])
    features.append(dataset)
    mol_ids = dataset.index.get_level_values(0).unique().values
    feature_mol_ids.append(mol_ids)
combined = pd.concat(features)
feature_mol_ids = np.concatenate(feature_mol_ids)

save_dir = "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/ECIFs/"
for subset in ["train", "val", "test"]:
    mol_ids = df[df.train_test == subset].mol_id.unique()
    mol_ids = np.intersect1d(mol_ids, feature_mol_ids)
    subset_features = combined.loc[mol_ids]
    subset_features.to_csv(f"{save_dir}/{subset}.csv")
    pass






