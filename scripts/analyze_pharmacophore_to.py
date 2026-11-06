import csv
import os
import sys
from collections import defaultdict
import numpy as np
from moleculekit.molecule import Molecule

# Same as analyze_pharmacophore.py but allows specifying output directory

RES_AROMATIC = {
    'PHE': ['CG','CD1','CD2','CE1','CE2','CZ'],
    'TYR': ['CG','CD1','CD2','CE1','CE2','CZ'],
    'TRP': ['CD2','CE2','CE3','CD1','NE1','CZ2','CZ3','CH2']
}
RES_HYDROPHOBIC = set(['ALA','VAL','LEU','ILE','MET','PRO'])

DONOR_ATOMS = {
    'LYS': ['NZ'],
    'ARG': ['NE','NH1','NH2'],
    'ASN': ['ND2'],
    'GLN': ['NE2'],
    'SER': ['OG'],
    'THR': ['OG1'],
    'TYR': ['OH'],
    'TRP': ['NE1'],
    'HIS': ['ND1','NE2'],
    'CYS': ['SG']
}
ACCEPTOR_ATOMS = {
    'ASP': ['OD1','OD2'],
    'GLU': ['OE1','OE2'],
    'ASN': ['OD1'],
    'GLN': ['OE1'],
    'SER': ['OG'],
    'THR': ['OG1'],
    'TYR': ['OH'],
    'HIS': ['ND1','NE2'],
    'CYS': ['SG']
}
POSITIVE_RES = {
    'LYS': ['NZ'],
    'ARG': ['CZ']
}
NEGATIVE_RES = {
    'ASP': ['CG'],
    'GLU': ['CD']
}
METAL_NAMES = set(['ZN','MG','MN','FE','CO','NI','CU','CA'])

def ring_centroid(coords):
    return coords.mean(axis=0)

def res_center(mol, sel):
    idx = mol.atomselect(sel)
    if not np.any(idx):
        return None
    return mol.coords[idx, :, 0].mean(axis=0)

def features_for_residue(mol, chain, resid, resname):
    feats = []
    base_sel = f"chain {chain} and resid {resid}"
    for atm in DONOR_ATOMS.get(resname, []):
        c = res_center(mol, base_sel + f" and name {atm}")
        if c is not None:
            feats.append(('DON', resname, atm, c))
    for atm in ACCEPTOR_ATOMS.get(resname, []):
        c = res_center(mol, base_sel + f" and name {atm}")
        if c is not None:
            feats.append(('ACC', resname, atm, c))
    for atm in POSITIVE_RES.get(resname, []):
        c = res_center(mol, base_sel + f" and name {atm}")
        if c is not None:
            feats.append(('POS', resname, atm, c))
    for atm in NEGATIVE_RES.get(resname, []):
        c = res_center(mol, base_sel + f" and name {atm}")
        if c is not None:
            feats.append(('NEG', resname, atm, c))
    if resname in RES_HYDROPHOBIC:
        c = res_center(mol, base_sel + " and protein and element C and not name CA C O N")
        if c is not None:
            feats.append(('HPH', resname, 'SC', c))
    if resname in RES_AROMATIC:
        ring_atoms = ' '.join(RES_AROMATIC[resname])
        idx = mol.atomselect(base_sel + f" and name {ring_atoms}")
        if np.any(idx):
            c = ring_centroid(mol.coords[idx,:,0])
            feats.append(('ARO', resname, 'RING', c))
    return feats

def collect_pocket_features(pdb, centers, radius=6.0, outdir='pharmacophore'):
    os.makedirs(outdir, exist_ok=True)
    mol = Molecule(pdb)
    metal_idx = mol.atomselect('not protein and element ' + ' '.join(sorted(METAL_NAMES)))
    metal_feats = []
    if np.any(metal_idx):
        for i in np.where(metal_idx)[0]:
            metal_feats.append(('MET', mol.resname[i].strip(), mol.name[i].strip(), mol.coords[i,:,0]))
    results = []
    for pi,(score,center) in enumerate(centers, start=1):
        c = np.array(center, dtype=float)
        prot = mol.atomselect('protein and noh')
        coords = mol.coords[prot,:,0]
        idx = np.where(prot)[0]
        d = np.linalg.norm(coords - c, axis=1)
        near = idx[d <= radius]
        triples = sorted(set((mol.chain[i], int(mol.resid[i]), mol.resname[i].strip()) for i in near))
        feats = []
        for ch,resi,resn in triples:
            feats.extend(features_for_residue(mol, ch, resi, resn))
        for tag,resn,atm,mpos in metal_feats:
            if np.linalg.norm(mpos - c) <= radius:
                feats.append((tag,resn,atm,mpos))
        csvpath = os.path.join(outdir, f'pocket{pi}_features.csv')
        with open(csvpath,'w',newline='') as f:
            w = csv.writer(f)
            w.writerow(['feature','resname','atom','chain','resid','x','y','z'])
            for ft,resn,atm,pt in feats:
                chain='A'; resid=-1; bestd=1e9
                for ch,resi,resname in triples:
                    if resname==resn:
                        rc = res_center(mol, f'chain {ch} and resid {resi} and name CA')
                        if rc is None:
                            rc = res_center(mol, f'chain {ch} and resid {resi}')
                        if rc is None:
                            continue
                        dd = np.linalg.norm(rc-pt)
                        if dd<bestd:
                            bestd=dd; chain=ch; resid=resi
                w.writerow([ft,resn,atm,chain,resid,float(pt[0]),float(pt[1]),float(pt[2])])
        pdbpath = os.path.join(outdir, f'pocket{pi}_features.pdb')
        with open(pdbpath,'w') as f:
            serial=1
            for ft,resn,atm,pt in feats:
                resname = {
                    'DON':'PH_DON','ACC':'PH_ACC','POS':'PH_POS','NEG':'PH_NEG',
                    'HPH':'PH_HPH','ARO':'PH_ARO','MET':'PH_MET'
                }[ft]
                f.write("HETATM%5d  %-3s %3s A%4d    %8.3f%8.3f%8.3f  1.00 10.00          C\n" % (
                    serial, atm[:3], resname, serial, pt[0], pt[1], pt[2]))
                serial+=1
        results.append({'pocket':pi,'score':score,'n_features':len(feats),'csv':csvpath,'pdb':pdbpath})
    with open(os.path.join(outdir,'summary.csv'),'w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['pocket','score','n_features','csv','pdb'])
        for r in results:
            w.writerow([r['pocket'],r['score'],r['n_features'],r['csv'],r['pdb']])
    return results

if __name__=='__main__':
    pdb = sys.argv[1]
    csvfile = sys.argv[2]
    outdir = sys.argv[3]
    radius = float(sys.argv[4]) if len(sys.argv)>4 else 6.0
    centers=[]
    with open(csvfile) as f:
        r=csv.DictReader(f)
        for row in r:
            score=float(row['scores'])
            center = np.fromstring(row['centers'].strip('[]'),sep=',')
            centers.append((score, center))
    collect_pocket_features(pdb, centers[:1], radius=radius, outdir=outdir)
