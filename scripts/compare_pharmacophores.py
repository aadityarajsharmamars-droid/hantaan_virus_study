import os, sys, csv
import pandas as pd
import numpy as np

csv1 = sys.argv[1]
csv2 = sys.argv[2]
outdir = sys.argv[3] if len(sys.argv)>3 else 'pharmacophore_compare'

os.makedirs(outdir, exist_ok=True)

def load(csvpath):
    df = pd.read_csv(csvpath)
    # Normalize feature codes
    df['feature'] = df['feature'].str.upper()
    return df

df1 = load(csv1)
df2 = load(csv2)

all_feats = sorted(set(df1['feature']).union(df2['feature']))

summary = []
for ft in all_feats:
    n1 = int((df1['feature']==ft).sum())
    n2 = int((df2['feature']==ft).sum())
    summary.append({'feature': ft, 'count_1': n1, 'count_2': n2, 'diff': n2 - n1})

overall = {
    'n_features_1': int(df1.shape[0]),
    'n_features_2': int(df2.shape[0]),
    'has_POS_1': bool((df1['feature']=='POS').any()),
    'has_POS_2': bool((df2['feature']=='POS').any()),
    'has_NEG_1': bool((df1['feature']=='NEG').any()),
    'has_NEG_2': bool((df2['feature']=='NEG').any()),
    'has_MET_1': bool((df1['feature']=='MET').any()),
    'has_MET_2': bool((df2['feature']=='MET').any()),
}

pd.DataFrame(summary).to_csv(os.path.join(outdir,'feature_counts.csv'), index=False)
pd.Series(overall).to_csv(os.path.join(outdir,'overall.csv'))

# Simple categorical overlap score (type-only Jaccard with multiplicities):
# overlap = sum(min(count1_ft, count2_ft)) / sum(max(count1_ft, count2_ft))
num = 0; den = 0
for ft in all_feats:
    c1 = int((df1['feature']==ft).sum()); c2 = int((df2['feature']==ft).sum())
    num += min(c1,c2); den += max(c1,c2)
score = float(num)/den if den>0 else 0.0
with open(os.path.join(outdir,'summary.txt'),'w') as f:
    f.write(f'Type-overlap score: {score:.3f} (1.0 identical by type counts)\n')
    f.write('Feature count differences (2 - 1):\n')
    for row in summary:
        f.write(f"  {row['feature']}: {row['count_2']} vs {row['count_1']} (diff {row['diff']})\n")
    f.write('\nFlags: POS/NEG/MET presence\n')
    f.write(str(overall)+'\n')
print('Wrote comparison to', outdir)
