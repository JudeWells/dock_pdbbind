import os

import torch
import wandb
from src.data import BaseDataset

def build_data_loaders(data_config):
    if os.getcwd().startswith("/Users/judewells"):
        data_dir = "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand"
    else:
        data_dir = "/SAN/orengolab/nsp13/data_for_protein_ligand"
    train_path = data_config.get('train_path',
                                   f"{data_dir}/ECIFs/train.csv")
    val_path = data_config.get('val_path',
                                    f"{data_dir}/ECIFs/val.csv")
    test_path = data_config.get('test_path',
                                     f"{data_dir}/ECIFs/test.csv")
    label_path = data_config.get('label_path',
                                 "pdb_bind.csv")

    data_loaders = []

    for feature_path in [train_path, val_path, test_path]:
        data = BaseDataset(feature_path, label_path)
        loader = torch.utils.data.DataLoader(data, batch_size=data_config.get("batch_size", 16), shuffle=True)
        data_loaders.append(loader)

    train_loader, val_loader, test_loader = data_loaders
    return train_loader, val_loader, test_loader

if __name__ == "__main__":
    train_loader = build_data_loaders(data_config={})


