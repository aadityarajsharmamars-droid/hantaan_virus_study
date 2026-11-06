import sys, os, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

csv_path = sys.argv[1]
out_sdf = sys.argv[2]

os.makedirs(os.path.dirname(out_sdf), exist_ok=True)

df = pd.read_csv(csv_path)
m = Chem.EditableMol(Chem.Mol())
coords = []
features = []
for _,row in df.iterrows():
    a = Chem.Atom(6)
    idx = m.AddAtom(a)
    coords.append((float(row['x']), float(row['y']), float(row['z'])))
    features.append(str(row['feature']))

mol = m.GetMol()
# set conformer
conf = Chem.Conformer(mol.GetNumAtoms())
for i,(x,y,z) in enumerate(coords):
    conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(x,y,z))
mol.AddConformer(conf, assignId=True)

for i,ft in enumerate(features):
    mol.SetProp(f'feat_{i}', ft)

w = Chem.SDWriter(out_sdf)
mol.SetProp('_Name', os.path.basename(csv_path).replace('.csv',''))
w.write(mol)
w.close()
print('Wrote', out_sdf)
