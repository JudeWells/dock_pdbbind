"""
Classifier that uses morgan fingerprints of ligands only
to classify the pose accuracy
Random forest regression model has R2 of 0.263 MAE of 1.54
Linear regression model had negative R2

corr coefficient between vina score and RMSD 0.017
corr coefficient between num rot bonds and RMSD 0.41
corr coefficient between num heavy atoms and RMSD 0.45
"""
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression

def prepare_data_for_classifier(train, test, RADIUS=3, NBITS=2048):
    x_train = np.array([np.array(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), RADIUS, nBits=NBITS)) for smi in train["smiles"]])
    x_test = np.array([np.array(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), RADIUS, nBits=NBITS)) for smi in test["smiles"]])
    y_train = train["rmsd"].values
    y_test = test["rmsd"].values
    return x_train, x_test, y_train, y_test

if __name__=="__main__":
    df = pd.read_csv("pdb_bind.csv")
    rmsd_vina = np.corrcoef(df.rmsd.values, df.vina_score.values)[0, 1]
    print("Corr coeff between vina score and RMSD", round(rmsd_vina,4))
    rmsd_num_rot = np.corrcoef(df.rmsd.values, df.num_rotatable_bonds)[0,1]
    print("Corr coeff between num rot bonds and RMSD", round(rmsd_num_rot,4))
    train = df[df["train_test"]=="train"]
    test = df[df["train_test"]=="test"]
    train = train.groupby("smiles").mean().reset_index()
    test = test.groupby("smiles").mean().reset_index()
    x_train, x_test, y_train, y_test = prepare_data_for_classifier(train, test)
    for model in [RandomForestRegressor(n_estimators=100, random_state=42), LinearRegression()]:
        model.fit(x_train, y_train)
        test_pred = model.predict(x_test)
        print("MSE:", np.mean((test_pred - y_test)**2))
        print("R2:", model.score(x_test, y_test))
        print("MAE:", np.mean(np.abs(test_pred - y_test)))