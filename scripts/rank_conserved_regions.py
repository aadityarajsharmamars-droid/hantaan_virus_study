#!/usr/bin/env python3
import os
import requests
import pandas as pd

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
CONSDIR = os.path.join(BASE, 'data', 'conservation')
OUTDIR = os.path.join(CONSDIR, 'rankings')
os.makedirs(OUTDIR, exist_ok=True)

USER_AGENT = {'User-Agent': 'PlayMoleculeAI/1.0 (Acellera)'}

REFS = {
    'N': {'acc': 'O36307', 'csv': os.path.join(CONSDIR, 'N_deltaH.csv'), 'offset': 0},
    'GPC': {'acc': 'Q9E006', 'csv': os.path.join(CONSDIR, 'GPC_deltaH.csv'), 'offset': 0},
    'L': {'acc': 'Q9E005', 'csv': os.path.join(CONSDIR, 'L_deltaH.csv'), 'offset': 0},
}

FEATURE_TYPES = {
    'Domain', 'Region', 'Motif', 'Zinc finger', 'Topological domain', 'Transmembrane',
    'Signal', 'Site', 'Active site', 'Metal binding', 'Binding site', 'Coiled coil'
}

all_rows = []

for prot, info in REFS.items():
    acc = info['acc']
    dh = pd.read_csv(info['csv'])
    dh = dh.set_index('pos')
    url = f'https://rest.uniprot.org/uniprotkb/{acc}.json'
    uj = requests.get(url, headers=USER_AGENT, timeout=60).json()
    for feat in uj.get('features', []):
        ftype = feat.get('type')
        if ftype not in FEATURE_TYPES:
            continue
        loc = feat.get('location', {})
        b = loc.get('start', {}).get('value')
        e = loc.get('end', {}).get('value')
        if not b or not e:
            continue
        b = int(b); e = int(e)
        if e < b:
            b, e = e, b
        desc = feat.get('description') or ''
        # slice deltaH
        sub = dh.loc[(dh.index >= b) & (dh.index <= e)]
        if sub.empty:
            continue
        median = float(sub['deltaH'].median())
        mean = float(sub['deltaH'].mean())
        minv = float(sub['deltaH'].min())
        maxv = float(sub['deltaH'].max())
        all_rows.append({
            'protein': prot,
            'uniprot': acc,
            'feature_type': ftype,
            'description': desc,
            'start': b,
            'end': e,
            'length': e - b + 1,
            'deltaH_median': median,
            'deltaH_mean': mean,
            'deltaH_min': minv,
            'deltaH_max': maxv,
        })

rank_df = pd.DataFrame(all_rows)
rank_df.sort_values(['deltaH_median','deltaH_mean'], ascending=[False, False], inplace=True)
rank_csv = os.path.join(OUTDIR, 'conserved_regions_ranked.csv')
rank_df.to_csv(rank_csv, index=False)

# Also produce least conserved ranking
least_df = rank_df.sort_values(['deltaH_median','deltaH_mean'], ascending=[True, True])
least_csv = os.path.join(OUTDIR, 'variable_regions_ranked.csv')
least_df.to_csv(least_csv, index=False)

print(rank_csv)