"""
Created by Jude Wells 2023-05-15
this script is made to run autodock vina on PDBbind data
we assume that proteins have been cleaned (remove waters etc)
"""

import os
import glob
import subprocess
import sys
import numpy as np
from rdkit import Chem

def calculate_bounding_box(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is not None:
            conf = mol.GetConformer()
            atom_positions = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
            center = atom_positions.mean(axis=0)
            sizes = atom_positions.ptp(axis=0) + 15 # peak to peak

            return center, sizes
    return None

def smiles_2_3d_mol2(smiles_string, out_path):
    with open(out_path + ".smi", "w") as f:
        f.write(smiles_string + "\n")
    # cmd_ = f"obabel {out_path}.smi -O {out_path}.mol2 --gen3d --best --canonical --conformers --weighted --nconf 50 --ff GAFF"
    # print(cmd_)
    subprocess.call(["obabel",
                     f"{out_path}.smi",
                     "-O", f"{out_path}.mol2",
                     "--gen3d",
                        "--best",
                        "--canonical",
                        "--conformers",
                        "--weighted",
                        "--nconf", "50",
                        "--ff", "GAFF"
                     ])
    return out_path + ".mol2"

def convert_sdf_to_smiles(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is not None:
            return Chem.MolToSmiles(mol)
    return None

if __name__=="__main__":
    task_index = int(sys.argv[1])
    batch_size = 100
    data_dir = "/SAN/orengolab/nsp13/dock_pdbbind/PDBbind/PDBBind_processed"
    print(f"data_dir: {data_dir}")
    for pdbid in sorted(os.listdir(data_dir))[task_index*batch_size:(task_index+1)*batch_size]:
        pdbid_path = os.path.join(data_dir, pdbid)
        if os.path.isdir(pdbid_path):
            try:
                print(f"pdbid_path: {pdbid_path}")
                mol2path = glob.glob(f"{pdbid_path}/*_ligand.sdf.mol2")
                unprocessed_ligand_sdf_file = glob.glob(f"{pdbid_path}/*_ligand.sdf")[0]

                if len(mol2path)==0:
                    print("converting sdf to mol2")
                    smiles_string = convert_sdf_to_smiles(unprocessed_ligand_sdf_file)
                    mol2path = smiles_2_3d_mol2(smiles_string, unprocessed_ligand_sdf_file)
                else:
                    mol2path = mol2path[0]
                    print(f"found existing mol2 file: {mol2path}")
                bounding_box = calculate_bounding_box(unprocessed_ligand_sdf_file)
                print(bounding_box)
                processed_ligand_file = mol2path + ".pdbqt"
                docked_ligand_file = os.path.join(pdbid_path, f"{pdbid}_docked_ligand.pdbqt")
                if os.path.exists(docked_ligand_file):
                    continue
                receptor_file = glob.glob(f"{pdbid_path}/*_processed.pdb.pdbqt")
                if len(receptor_file) == 0:
                    pdb_filepath = f"{pdbid_path}/{pdbid}_protein_processed.pdb"
                    assert os.path.exists(pdb_filepath)
                    print("preparing receptor:", pdb_filepath)
                    subprocess.call([
                        "prepare_receptor4.py",
                        "-r" , pdb_filepath,
                        "-o", f"{pdb_filepath}.pdbqt",
                        "-A", "hydrogens",
                        "-U", "nphs_lps", "-v"])
                receptor_file = glob.glob(f"{pdbid_path}/*_processed.pdb.pdbqt")
                if len(receptor_file) == 0:
                    print("failed to process receptor")
                    continue
                receptor_file = receptor_file[0]
                if not os.path.exists(processed_ligand_file):
                    print("preparing ligand:", mol2path)
                    subprocess.call([
                        "prepare_ligand4.py",
                        "-l", mol2path,
                        "-o", processed_ligand_file,
                        "-U", "nphs_lps", "-v"
                    ])
                    print(f"new ligand file: {processed_ligand_file}")
                    print(f"new ligand file exists: {os.path.exists(processed_ligand_file)}")
                if bounding_box is not None:
                    center, sizes = bounding_box
                    print("running vina")
                    subprocess.call(["vina",
                                     "--receptor", receptor_file,
                                     "--ligand", processed_ligand_file,
                                     "--out", docked_ligand_file,
                                     "--center_x", str(center[0]),
                                     "--center_y", str(center[1]),
                                     "--center_z", str(center[2]),
                                     "--size_x", str(sizes[0]),
                                     "--size_y", str(sizes[1]),
                                     "--size_z", str(sizes[2])])
            except:
                pass

