import os
import sys
import json
import time
import math
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

CSV_PATH = sys.argv[1] if len(sys.argv) > 1 else 'user-uploads/chembl_inhibitors.csv'
OUT_SDF = sys.argv[2] if len(sys.argv) > 2 else 'ligands/ligands.sdf'
MAP_CSV = sys.argv[3] if len(sys.argv) > 3 else 'ligands/ligand_map.csv'

os.makedirs(os.path.dirname(OUT_SDF), exist_ok=True)
os.makedirs(os.path.dirname(MAP_CSV), exist_ok=True)

def fetch_smiles(chembl_id, retries=3, delay=0.2):
    url = f'https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json'
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            data = r.json()
            ms = data.get('molecule_structures') or {}
            smi = ms.get('canonical_smiles')
            return smi
        except Exception as e:
            if attempt == retries-1:
                raise
            time.sleep(delay * (2**attempt))
    return None


def build_3d_from_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.numThreads = 0
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        pass
    return mol


def main():
    df = pd.read_csv(CSV_PATH)
    # unique molecules by molecule_chembl_id
    ids = df['molecule_chembl_id'].dropna().astype(str).unique().tolist()
    ids.sort()
    writer = Chem.SDWriter(OUT_SDF)
    written = 0
    rows = []
    for cid in ids:
        try:
            smi = fetch_smiles(cid)
        except Exception as e:
            print(f'ERROR fetching {cid}: {e}', file=sys.stderr)
            continue
        if not smi:
            print(f'No SMILES for {cid}', file=sys.stderr)
            continue
        mol = build_3d_from_smiles(smi)
        if mol is None:
            print(f'Failed 3D build for {cid}', file=sys.stderr)
            continue
        mol.SetProp('_Name', cid)
        # annotate with minimal pchembl from dataset (best potency)
        sub = df[df['molecule_chembl_id'].astype(str)==cid]
        best_pchembl = pd.to_numeric(sub.get('pchembl_value'), errors='coerce').max()
        if not math.isnan(best_pchembl) if best_pchembl is not None else False:
            mol.SetProp('pchembl_max', str(best_pchembl))
        writer.write(mol)
        written += 1
        rows.append({'molecule_chembl_id': cid, 'smiles': smi})
    writer.close()
    pd.DataFrame(rows).to_csv(MAP_CSV, index=False)
    print(f'Written {written} ligands to {OUT_SDF}')

if __name__ == '__main__':
    main()
