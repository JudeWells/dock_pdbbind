import torch
import wandb
from src.data import BaseDataset

def build_data_loaders(data_config):
    train_path = data_config.get('train_path',
                                   "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/ECIFs/train.csv")
    val_path = data_config.get('val_path',
                                    "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/ECIFs/val.csv")
    test_path = data_config.get('test_path',
                                     "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/ECIFs/test.csv")
    label_path = data_config.get('label_path',
                                 "/Users/judewells/Documents/dataScienceProgramming/dock_pdbbind/pdb_bind.csv")

    data_loaders = []

    for feature_path in [train_path, val_path, test_path]:
        data = BaseDataset(feature_path, label_path)
        loader = torch.utils.data.DataLoader(data, batch_size=data_config.get("batch_size", 16), shuffle=True)
        data_loaders.append(loader)

    train_loader, val_loader, test_loader = data_loaders
    return train_loader, val_loader, test_loader

if __name__ == "__main__":
    train_loader = build_data_loaders(data_config={})


