"""
Calculate RMS distances
between docked poses and experimental poses
"""

import os
import shutil
import glob
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import RDConfig
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'CalcLigRMSD'))
from utils.CalcLigRMSD import *

def load_sdf(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        return mol

def convert_pdbqt_to_sdf(pdbqt_file):
    subprocess.call(["obabel",
                     pdbqt_file,
                     "-O", pdbqt_file.replace(".pdbqt", ".sdf"),
                     ])
    return pdbqt_file.replace(".pdbqt", ".sdf")


def rmsd_docked(mol, dockedmol):
    """ function for calculating the rmsd between the atoms of the merged and
    docked poses
    """
    confmols=[]
    rmsds=[]
    for conf in dockedmol.GetConformers():
        molx=Chem.RWMol(dockedmol)
        molx.RemoveAllConformers()
        molx.AddConformer(conf)
        confmols.append(molx)
    for i,cm in enumerate(confmols):
        rmsd=CalcLigRMSD(cm,mol, rename_lig2 = False)
        rmsds.append(rmsd)
    return rmsds[0]

def main(data_dir):
    new_rows = []
    success_counter = 0
    bad_ligands = ["3wd1", "6cdl", "3b0w", "3tyv", "6dne",
                   "2bcd", "1fzm", "2l65", "1ibg", "4z88",
                   "6f08", "3run", "3pp7", "3u18", "4ret",
                   "3b5j", "1fkf", "3wsy", "6djc", "3upf",
                   "1fq5", "1o5p", "1fkb", "6ce2",
                   ]
    if os.path.exists("rmsd.csv"):
        df = pd.read_csv("rmsd.csv")
    else:
        df = pd.DataFrame(columns=["pdb_id", "rmsd", "vina_rank", "smiles"])
    for i, dir_name in enumerate(os.listdir(data_dir)):
        if dir_name in bad_ligands or dir_name in df["pdb_id"].values:
            continue
        dir_path = os.path.join(data_dir, dir_name)
        if os.path.isdir(dir_path):
            try:
                original_ligand_file = glob.glob(f"{dir_path}/{dir_name}_ligand.sdf")[0]
                docked_ligand_file = glob.glob(f"{dir_path}/*docked_ligand.sdf")[0]
            except IndexError:
                try:
                    docked_ligand_file = glob.glob(f"{dir_path}/*docked_ligand.pdbqt")[0]
                    docked_ligand_file = convert_pdbqt_to_sdf(docked_ligand_file)
                except IndexError:
                    continue
            original = load_sdf(original_ligand_file)
            if original is None:
                continue
            natoms = original.GetNumAtoms()
            print(dir_name, natoms)
            if natoms > 100:
                print(f"Skipping {dir_name} because it has {natoms} heavy atoms")
                continue
            # original = Chem.RemoveHs(original)
            original_smiles = Chem.MolToSmiles(original)
            suppl = Chem.SDMolSupplier(docked_ligand_file)
            for j, docked in enumerate(suppl):
                if docked is None:
                    rmsd = -1
                else:
                    rmsd = CalcLigRMSD(docked, original, rename_lig2=False)
                    if rmsd == 0:
                        print( Chem.MolToSmiles(original))
                        print( Chem.MolToSmiles(docked))
                new_row = {
                    "pdb_id": dir_name,
                    "rmsd": rmsd,
                    "vina_rank": j,
                    "smiles":original_smiles,
                }
                new_rows.append(new_row)
            success_counter += 1
            if i % 10 == 0:
                print(f"Processed {i} files, {len(new_rows)} rows in dataframe, {success_counter} successes")
        if len(new_rows) > 0:
            new_df = pd.DataFrame(new_rows)
            # combine new rows with old rows
            df = pd.concat([df, new_df])
            df.to_csv("rmsd.csv", index=False)
            new_rows = []
    return df

if __name__=="__main__":
    datadir = "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/PDBBind_processed"
    df = main(datadir)

