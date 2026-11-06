import os, sys, csv
import numpy as np
from moleculekit.molecule import Molecule

FS6_PDB='receptors/6FS6/output.pdb'
Q99_PDB='receptors/6Q99/output.pdb'
FS6_CENTER_SDF='templates/6FS6_bxa_center.sdf'
Q99_CENTER_CSV='pharmacophore/6Q99_pocket1_center.csv'
OUTDIR='analysis/active_site'
os.makedirs(OUTDIR, exist_ok=True)

METALS=set(['ZN','MG','MN','FE','CO','NI','CU','CA'])


def read_center_from_sdf(path):
    with open(path) as f:
        for line in f:
            if line.strip().endswith('V2000'):
                coords=[]
                # next line is atom record
            if line.strip().startswith('M  END'):
                break
    # simpler: parse first atom line with floats at columns 0-30
    with open(path) as f:
        lines=f.readlines()
    for i,l in enumerate(lines):
        if 'V2000' in l:
            rec=lines[i+1]
            x=float(rec[0:10]); y=float(rec[10:20]); z=float(rec[20:30])
            return np.array([x,y,z])
    raise RuntimeError('Center not found in '+path)


def read_center_from_csv(path):
    import pandas as pd
    df=pd.read_csv(path)
    row=df.iloc[0]
    ctr=np.fromstring(row['centers'].strip('[]'), sep=',')
    return ctr


def select_pocket_residues(mol, center, radius=10.0):
    prot = mol.atomselect('protein and name CA')
    idx = np.where(prot)[0]
    coords = mol.coords[idx,:,0]
    d = np.linalg.norm(coords - center, axis=1)
    near = idx[d <= radius]
    triples = sorted(set((mol.chain[i], int(mol.resid[i]), mol.resname[i].strip()) for i in near))
    return triples


def metal_coordination(mol, maxdist=2.6):
    metal_idx = mol.atomselect('not protein and element ' + ' '.join(sorted(METALS)))
    rows=[]
    for i in np.where(metal_idx)[0]:
        mcoord = mol.coords[i,:,0]
        mname = mol.name[i].strip(); mres=mol.resname[i].strip(); mchain=mol.chain[i].strip(); mresid=int(mol.resid[i])
        # potential liganding atoms: O/N from protein side chains
        cand = mol.atomselect('protein and (element O or element N)')
        ci = np.where(cand)[0]
        for j in ci:
            d = np.linalg.norm(mol.coords[j,:,0]-mcoord)
            if d <= maxdist:
                rows.append({
                    'metal': mres, 'metal_name': mname, 'metal_chain': mchain, 'metal_resid': mresid,
                    'lig_resname': mol.resname[j].strip(), 'lig_chain': mol.chain[j].strip(), 'lig_resid': int(mol.resid[j]), 'lig_atom': mol.name[j].strip(), 'distance': float(d)
                })
    return rows


def kabsch(P, Q):
    # P -> Q
    muP=P.mean(axis=0); muQ=Q.mean(axis=0)
    P0=P-muP; Q0=Q-muQ
    H=P0.T @ Q0
    U,S,Vt=np.linalg.svd(H)
    R=Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1,:]*=-1
        R=Vt.T @ U.T
    t=muQ - R @ muP
    return R,t


def pocket_backbone_rmsd(fs6mol, q99mol, center, rad=10.0):
    sel_ca_fs6 = fs6mol.atomselect('protein and name CA')
    sel_ca_q99 = q99mol.atomselect('protein and name CA')
    i_fs6 = np.where(sel_ca_fs6)[0]
    i_q99 = np.where(sel_ca_q99)[0]
    c_fs6 = fs6mol.coords[i_fs6,:,0]
    c_q99 = q99mol.coords[i_q99,:,0]
    d_fs6 = np.linalg.norm(c_fs6 - center, axis=1)
    d_q99 = np.linalg.norm(c_q99 - center, axis=1)
    P = c_fs6[d_fs6<=rad]
    Q = c_q99[d_q99<=rad]
    # Trim to same size by nearest sampling
    n=min(len(P), len(Q))
    if n<10:
        return None
    P=P[:n]; Q=Q[:n]
    R,t=kabsch(P,Q)
    Palign=(R @ P.T).T + t
    rmsd=np.sqrt(((Palign-Q)**2).sum()/n)
    return float(rmsd), n


def main():
    fs6 = Molecule(FS6_PDB)
    q99 = Molecule(Q99_PDB)
    fs6_ctr = read_center_from_sdf(FS6_CENTER_SDF)
    q99_ctr = read_center_from_csv(Q99_CENTER_CSV)

    # Residues near pockets
    fs6_res = select_pocket_residues(fs6, fs6_ctr, radius=10.0)
    q99_res = select_pocket_residues(q99, q99_ctr, radius=10.0)

    with open(os.path.join(OUTDIR,'fs6_pocket_residues.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['chain','resid','resname']); w.writerows(fs6_res)
    with open(os.path.join(OUTDIR,'q99_pocket_residues.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['chain','resid','resname']); w.writerows(q99_res)

    # Metal coordination tables
    fs6_met = metal_coordination(fs6)
    q99_met = metal_coordination(q99)
    import pandas as pd
    pd.DataFrame(fs6_met).to_csv(os.path.join(OUTDIR,'fs6_metal_coordination.csv'), index=False)
    pd.DataFrame(q99_met).to_csv(os.path.join(OUTDIR,'q99_metal_coordination.csv'), index=False)

    # Pocket backbone RMSD around Andes center to test local structural similarity
    res = pocket_backbone_rmsd(fs6, q99, q99_ctr, rad=10.0)
    with open(os.path.join(OUTDIR,'pocket_backbone_rmsd.txt'),'w') as f:
        if res is None:
            f.write('Insufficient CA atoms within radius for RMSD')
        else:
            rmsd,n = res
            f.write(f'Pocket backbone RMSD within 10A of Andes center: {rmsd:.3f} A using {n} CA pairs')

if __name__=='__main__':
    main()
