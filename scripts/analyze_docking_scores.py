import os
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon

PA_DIR = 'docking/PA_6FS6'
ANDV_DIR = 'docking/ANDV_6Q99_pocket1'
OUT_DIR = 'analysis'
os.makedirs(OUT_DIR, exist_ok=True)


def read_scores(d):
    path = os.path.join(d, 'scores.csv')
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    # AceDock produces columns: Name, Score RDock, Path
    df = df.rename(columns={'Name':'name', 'Score RDock':'score'})
    # Derive ligand id from Path (outligXX.sdf corresponds to order of input) - keep 'name' as outligN
    return df[['name','score']]

pa = read_scores(PA_DIR)
andv = read_scores(ANDV_DIR)

if pa is None or andv is None:
    missing = ['PA' if pa is None else None, 'ANDV' if andv is None else None]
    missing = [m for m in missing if m]
    with open(os.path.join(OUT_DIR,'status.txt'),'w') as w:
        w.write('Missing: ' + ','.join(missing))
    print('Waiting for jobs to finish. Missing:', ','.join(missing))
    raise SystemExit(0)

merged = pd.merge(pa, andv, on='name', suffixes=('_PA','_ANDV'))
merged['delta'] = merged['score_ANDV'] - merged['score_PA']
merged.to_csv(os.path.join(OUT_DIR,'merged_scores.csv'), index=False)

# Paired test (lower is better).
try:
    stat, p = wilcoxon(merged['score_PA'], merged['score_ANDV'], alternative='two-sided')
except ValueError:
    stat, p = np.nan, np.nan

summary = {
    'n_common': int(merged.shape[0]),
    'median_PA': float(merged['score_PA'].median()),
    'median_ANDV': float(merged['score_ANDV'].median()),
    'median_delta_ANDV_minus_PA': float(merged['delta'].median()),
    'wilcoxon_stat': float(stat) if not np.isnan(stat) else None,
    'wilcoxon_p': float(p) if not np.isnan(p) else None
}
pd.Series(summary).to_csv(os.path.join(OUT_DIR,'summary.csv'))

# Dual top quartile binders per target
qPA = merged['score_PA'].quantile(0.25)
qAN = merged['score_ANDV'].quantile(0.25)
merged['topQ_PA'] = merged['score_PA'] <= qPA
merged['topQ_ANDV'] = merged['score_ANDV'] <= qAN
both = merged[merged['topQ_PA'] & merged['topQ_ANDV']].copy()
both.to_csv(os.path.join(OUT_DIR,'dual_top_quartile.csv'), index=False)

print('Analysis done. N common:', summary['n_common'], 'p=', summary['wilcoxon_p'])
