import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import glob
import os

def extract_docking_scores(filename):
    scores = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('REMARK VINA RESULT:'):
                score = float(line.split()[3])
                scores.append(score)
    return scores

def add_mol_features(df):
    """
    Adds number of heavy atoms and number of rotatable bonds to the dataframe
    :param df:
    :return:
    """
    df['num_heavy_atoms'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x).GetNumHeavyAtoms())
    df['num_rotatable_bonds'] = df['smiles'].apply(lambda x: rdMolDescriptors.CalcNumRotatableBonds(Chem.MolFromSmiles(x)))
    return df

if __name__ == "__main__":
    df = pd.read_csv("pdb_bind.csv")
    df = add_mol_features(df)
    data_dir = "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/PDBBind_processed"
    for i, dir_name in enumerate(os.listdir(data_dir)):
        dir_path = os.path.join(data_dir, dir_name)
        if os.path.isdir(dir_path):
            try:
                filename = glob.glob(f"{dir_path}/*docked_ligand.pdbqt")[0]
            except IndexError:
                continue
        scores = extract_docking_scores(filename)
        for j, score in enumerate(scores):
            df.loc[(df['pdb_id'] == dir_name) & (df['vina_rank'] == j), 'vina_score'] = score
        print(f'The docking scores are: {scores}')
    df.to_csv("pdb_bind.csv", index=False)