import sys, os, pandas as pd
from rdkit import Chem

def write_subset(df, feat, out_sdf):
    sub = df[df['feature'].str.upper()==feat]
    if sub.empty:
        return False
    m = Chem.EditableMol(Chem.Mol())
    coords=[]
    for _,row in sub.iterrows():
        idx = m.AddAtom(Chem.Atom(6))
        coords.append((float(row['x']), float(row['y']), float(row['z'])))
    mol = m.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i,(x,y,z) in enumerate(coords):
        conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(x,y,z))
    mol.AddConformer(conf, assignId=True)
    mol.SetProp('_Name', os.path.basename(out_sdf).replace('.sdf',''))
    os.makedirs(os.path.dirname(out_sdf), exist_ok=True)
    w = Chem.SDWriter(out_sdf)
    w.write(mol)
    w.close()
    return True

if __name__=='__main__':
    csv_path = sys.argv[1]
    out_prefix = sys.argv[2]  # e.g., viewer_files/6FS6_ph4_
    feats = ['DON','ACC','POS','NEG','HPH','ARO','MET']
    df = pd.read_csv(csv_path)
    made=[]
    for ft in feats:
        out_sdf = f"{out_prefix}{ft}.sdf"
        if write_subset(df, ft, out_sdf):
            made.append(out_sdf)
    print('\n'.join(made))
