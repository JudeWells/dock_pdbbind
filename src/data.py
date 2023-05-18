import torch
import numpy as np
import pandas as pd

class BaseDataset(torch.utils.data.Dataset):
    def __init__(self, features_csv="/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/ECIFs/pdbbind_docked_2",
                 labels_csv="/Users/judewells/Documents/dataScienceProgramming/dock_pdbbind/pdb_bind.csv"):
        super().__init__()
        dataset = pd.read_csv(features_csv, index_col=[0, 1])
        mol_ids = dataset.index.get_level_values(0).unique()
        size = len(mol_ids), len(dataset.index.get_level_values(1).unique()), dataset.shape[1]
        print("Dataset size", size)
        self.dataset = dataset
        labels = pd.read_csv(labels_csv)
        labels["mol_id"] = labels.pdb_id + labels.vina_rank.apply(lambda x:"_ligand_" + str(x + 1))
        intersection = set(mol_ids).intersection(set(labels.mol_id))
        self.mol_ids = np.array(list(intersection))
        self.labels = labels.set_index("mol_id")
        self.labels = self.labels.loc[self.mol_ids]
        self.dataset = self.dataset.loc[self.mol_ids]

    def __len__(self):
        return len(self.mol_ids)

    def __getitem__(self, idx):
        mol_id = self.mol_ids[idx]
        x = self.dataset.loc[mol_id].values.flatten().astype(np.float32)
        y = self.labels.loc[mol_id, "rmsd"].astype(np.float32).reshape((1))
        return x, y