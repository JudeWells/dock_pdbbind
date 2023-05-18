import pandas as pd
import numpy as np
import torch

if __name__== "__main__":
    training_data_path = "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/ECIFs/pdbbind_docked_2"
    dataset = pd.read_csv(training_data_path, index_col=[0, 1])
    mol_id = dataset.index.get_level_values(0).unique()
    size = len(mol_id), len(dataset.index.get_level_values(1).unique()), dataset.shape[1]
    print(size)

    dataset_mols = dataset.values.reshape(size)
    dataset_mols = np.reshape(dataset_mols, [-1, size[1] * size[2]])
    bp=1