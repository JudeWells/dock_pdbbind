
import os
import glob
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
# from openbabel import pybel
import subprocess
# import py3Dmol

def convert_sdf_to_smiles(sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is not None:
            return Chem.MolToSmiles(mol)
    return None

def generate_3d_coords_from_smiles(smiles_string, output_sdf_file):
    mol = pybel.readstring("smi", smiles_string)
    mol.make3D()
    mol.write("sdf", output_sdf_file)



def mol2_2_pdbqt(mol2path):
    subprocess.call(["prepare_ligand4.py",
                     "-l",  mol2path,
                     "-o", mol2path + ".pdbqt",
                     "-U", "nphs_lps",
                     "-v"])
    return mol2path + ".pdbqt"


def dock_ligand(receptor_file, ligand_file, output_file, center, sizes):
    subprocess.call(["vina",
                     "--receptor", receptor_file,
                     "--ligand", ligand_file,
                     "--out", output_file,
                     "--center_x", str(center[0]),
                     "--center_y", str(center[1]),
                     "--center_z", str(center[2]),
                     "--size_x", str(sizes[0]),
                     "--size_y", str(sizes[1]),
                     "--size_z", str(sizes[2])])

def calculate_rmsd(original_ligand_file, docked_ligand_file):
    original_ligand = next(pybel.readfile("sdf", original_ligand_file))
    docked_ligand = next(pybel.readfile("sdf", docked_ligand_file))

    original_ligand.OBMol.AddHydrogens()
    docked_ligand.OBMol.AddHydrogens()

    align = pybel.ob.OBAlign(False, False)
    align.SetRefMol(original_ligand.OBMol)
    align.SetTargetMol(docked_ligand.OBMol)
    align.Align()

    return align.GetRMSD()

def main():
    pdbbind_dir = "/Users/judewells/Documents/dataScienceProgramming/data_for_protein_ligand/PDBbind/refined-set"
    results = []

    for dir_name in os.listdir(pdbbind_dir):
        dir_path = os.path.join(pdbbind_dir, dir_name)

        if os.path.isdir(dir_path):
            receptor_file = glob.glob(f"{dir_path}/*_protein.pdb")[0]
            original_ligand_file = glob.glob(f"{dir_path}/*_ligand.sdf")[0]
            docked_ligand_file = glob.glob(f"{dir_path}/*_ligand.pdbqt")[0]

            smiles_string = convert_sdf_to_smiles(original_ligand_file)

            mol2path = smiles_2_3d_mol2(smiles_string, original_ligand_file)



            generate_3d_coords_from_smiles(smiles_string, docked_ligand_file)

            dock_ligand(receptor_file, docked_ligand_file, docked_ligand_file)

            rmsd = calculate_rmsd(original_ligand_file, docked_ligand_file)

            results.append({"SMILES": smiles_string, "RMSD": rmsd})

    df = pd.DataFrame(results)
    df.to_csv("results.csv")

if __name__ == "__main__":
    main()