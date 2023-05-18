import os
import glob
import pandas as pd

data_dir = "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/PDBBind_processed"
success_counter = 0
fail_counter = 0
for i, dir_name in enumerate(os.listdir(data_dir)):
    dir_path = os.path.join(data_dir, dir_name)
    if os.path.isdir(dir_path):
        try:
            original_ligand_file = glob.glob(f"{dir_path}/*_ligand.sdf")[0]
            docked_ligand_file = glob.glob(f"{dir_path}/*.mol2.pdbqt")[0]
            success_counter += 1
        except IndexError:
            fail_counter += 1
            continue

print(f"Success: {success_counter}")
print(f"Fail: {fail_counter}")