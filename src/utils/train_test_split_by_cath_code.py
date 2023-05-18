import pandas as pd

cath_path = "/Users/judewells/Documents/dataScienceProgramming/domdet/data/cath_chain_topology_class.pickle"

cath = pd.read_pickle(cath_path)
cath["pdb_id"] = cath.chain_id.apply(lambda x: x[:-1])
df = pd.read_csv("rmsd.csv")
df = df[df.rmsd != -1]
intersect = set(df.pdb_id).intersection(set(cath.pdb_id))
print(len(intersect), "out of", len(set(df.pdb_id)))
cat2id = {}
id2cat = {}
for i, pdb_id in enumerate(df.pdb_id.unique()):
    if i % 100 == 0:
        print(i)
    if pdb_id in intersect:
        cats = [c for sublist in cath[cath.pdb_id == pdb_id].CAT.values for c in sublist]
        cats = list(set(cats))
        for c in cats:
            if c not in cat2id:
                cat2id[c] = []
            if pdb_id not in cat2id[c]:
                cat2id[c].append(pdb_id)
        df.loc[df.pdb_id == pdb_id, "CAT"] = "|".join(cats)
        df.loc[df.pdb_id == pdb_id, "n_CATs"] = len(cats)
        id2cat[pdb_id] = cats
cat_counts = pd.DataFrame([{'CAT': k, "count": len(v)} for k,v in cat2id.items()])
cat_counts = cat_counts.sort_values("count", ascending=False)
train_cats = []
val_cats = []
test_cats = []
for i, row in cat_counts.iterrows():
    if i % 4 == 0 or i % 7 == 0:
        test_cats.append(row.CAT)
    else:
        if i % 5 == 0 and i % 10 != 0:
            val_cats.append(row.CAT)
        else:
            train_cats.append(row.CAT)

df = df[df.CAT.notnull()].copy()
print("n_train cats", len(train_cats), "n test cats", len(test_cats), "n val cats", len(val_cats))

train_pdbs = []
for c in train_cats:
    train_pdbs.extend(cat2id[c])

val_pdbs = []
for c in val_cats:
    val_pdbs.extend(cat2id[c])

test_pdbs = list(set(df.pdb_id.unique()) - (set(train_pdbs).union(set(val_pdbs))))
print("n train pdbs", len(train_pdbs), "n val pdbs", len(val_pdbs), "n test pdbs", len(test_pdbs))
for i, row in df.iterrows():
    if row.pdb_id in train_pdbs:
        df.loc[i, "train_test"] = "train"
    elif row.pdb_id in val_pdbs:
        df.loc[i, "train_test"] = "val"
    else:
        df.loc[i, "train_test"] = "test"
df.to_csv("pdb_bind2.csv", index=False)
bp=1



