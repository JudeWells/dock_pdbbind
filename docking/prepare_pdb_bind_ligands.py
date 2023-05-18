import os
import glob
from main import convert_sdf_to_smiles, smiles_2_3d_mol2

if __name__=="__main__":
    pdbbind_dir = "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/PDBbind/refined-set"
    fail_count = 0
    for dir_name in os.listdir(pdbbind_dir):
        dir_path = os.path.join(pdbbind_dir, dir_name)

        if os.path.isdir(dir_path):
            try:
                receptor_file = glob.glob(f"{dir_path}/*_protein.pdb")[0]
                original_ligand_file = glob.glob(f"{dir_path}/*_ligand.sdf")[0]
                smiles_string = convert_sdf_to_smiles(original_ligand_file)
                mol2path = smiles_2_3d_mol2(smiles_string, original_ligand_file)
            except:
                fail_count += 1
                print(f"Failed on {dir_path}")



