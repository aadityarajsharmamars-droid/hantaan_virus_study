#!/usr/bin/env python3
import os
import pandas as pd

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
CONSDIR = os.path.join(BASE, 'data', 'conservation')
RANKDIR = os.path.join(CONSDIR, 'rankings')
OUTDIR = os.path.join(BASE, 'docs')
os.makedirs(OUTDIR, exist_ok=True)

# Overall mean ΔH per protein
def overall_ranking():
    rows = []
    for name in ['N','Gn','Gc','L']:
        path = os.path.join(CONSDIR, f'{name}_deltaH.csv')
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path)
        rows.append((name, float(df['deltaH'].mean()), float(df['deltaH'].median()), float(df['deltaH'].min()), len(df)))
    rows.sort(key=lambda x: x[1], reverse=True)
    return rows

# Top/Bottom 10 annotated features by ΔH median
def top_features(csv_path, top=True, k=10):
    df = pd.read_csv(csv_path)
    order = df.sort_values(['deltaH_median','deltaH_mean'], ascending=[not top, not top])
    if top:
        order = df.sort_values(['deltaH_median','deltaH_mean'], ascending=[False, False])
    else:
        order = df.sort_values(['deltaH_median','deltaH_mean'], ascending=[True, True])
    return order.head(k)

# Sliding-window minima for each protein (window=20)
def sliding_minima(name, w=20, n=3):
    path = os.path.join(CONSDIR, f'{name}_deltaH.csv')
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path)
    vals = []
    for i in range(0, len(df)-w+1):
        seg = df.iloc[i:i+w]
        vals.append((int(seg['pos'].iloc[0]), int(seg['pos'].iloc[-1]), float(seg['deltaH'].mean())))
    vals.sort(key=lambda x: x[2])
    return vals[:n]

# Build text
lines = []
lines.append('Andes orthohantavirus protein conservation ranking (ΔH, bits)')
lines.append('')

# Overall ranking
lines.append('Overall protein ranking by mean ΔH (higher = more conserved)')
rank = overall_ranking()
for i,(name,meanv,medv,minv,count) in enumerate(rank,1):
    lines.append(f"{i}) {name}: mean {meanv:.3f}, median {medv:.3f}, min {minv:.3f}, n={count}")
lines.append('')

# Top conserved annotated features
cons_csv = os.path.join(RANKDIR, 'conserved_regions_ranked.csv')
var_csv = os.path.join(RANKDIR, 'variable_regions_ranked.csv')
if os.path.exists(cons_csv):
    lines.append('Top 10 most conserved annotated features (by median ΔH)')
    topc = top_features(cons_csv, top=True, k=10)
    for _,r in topc.iterrows():
        desc = r['description'] if pd.notna(r['description']) else ''
        lines.append(f"- {r['protein']} | {r['feature_type']} | {desc} | {int(r['start'])}-{int(r['end'])} | ΔHmed={r['deltaH_median']:.3f}")
    lines.append('')

# Most variable annotated features
if os.path.exists(var_csv):
    lines.append('Top 10 most variable annotated features (lowest median ΔH)')
    bot = top_features(var_csv, top=False, k=10)
    for _,r in bot.iterrows():
        desc = r['description'] if pd.notna(r['description']) else ''
        lines.append(f"- {r['protein']} | {r['feature_type']} | {desc} | {int(r['start'])}-{int(r['end'])} | ΔHmed={r['deltaH_median']:.3f}")
    lines.append('')

# Sliding-window minima per protein
lines.append('Most variable 20-residue windows per protein (lowest mean ΔH)')
for name in ['N','Gn','Gc','L']:
    wins = sliding_minima(name, w=20, n=3)
    if not wins:
        continue
    lines.append(f"- {name}:")
    for s,e,m in wins:
        lines.append(f"    {s}-{e} | mean ΔH={m:.2f}")
lines.append('')

# Target selection rationale (brief)
lines.append('Target selection (small-molecule)')
lines.append('1) L endonuclease (N-terminus of L): essential, highly conserved catalytic site, druggable metal pocket (6Q99).')
lines.append('2) L RdRp catalytic motifs: conserved but larger/dynamic site; nucleoside/non-nucleoside strategies.')
lines.append('3) Nucleoprotein: conserved but surfaces often shallow/PPIs; possible oligomerization/RNA-binding disruptors.')
lines.append('4) Gn/Gc: better for biologics/entry blockers; small-molecule pockets less established.')

outpath = os.path.join(OUTDIR, 'ranking_summary.txt')
with open(outpath, 'w') as f:
    f.write('\n'.join(lines))
print(outpath)
