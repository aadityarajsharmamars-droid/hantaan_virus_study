#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
CONSDIR = os.path.join(BASE, 'data', 'conservation')
PLOTDIR = os.path.join(CONSDIR, 'plots')
os.makedirs(PLOTDIR, exist_ok=True)

series = {
    'N': os.path.join(CONSDIR, 'N_deltaH.csv'),
    'Gn': os.path.join(CONSDIR, 'Gn_deltaH.csv'),
    'Gc': os.path.join(CONSDIR, 'Gc_deltaH.csv'),
    'L': os.path.join(CONSDIR, 'L_deltaH.csv'),
}

for name, path in series.items():
    df = pd.read_csv(path)
    fig, ax = plt.subplots(figsize=(10,3))
    ax.plot(df['pos'], df['deltaH'], color='#1f77b4', linewidth=1.0)
    ax.set_title(f'Andes orthohantavirus {name} conservation (ΔH)')
    ax.set_xlabel('Position')
    ax.set_ylabel('ΔH (bits)')
    ax.set_ylim(0, 4.4)
    ax.grid(True, linestyle='--', alpha=0.3)
    out = os.path.join(PLOTDIR, f'{name}_deltaH.png')
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)

# Also a combined plot for Gn and Gc along GPC coordinate system
try:
    gpc = pd.read_csv(os.path.join(CONSDIR, 'GPC_deltaH.csv'))
    fig, ax = plt.subplots(figsize=(12,3))
    ax.plot(gpc['pos'], gpc['deltaH'], color='gray', linewidth=0.8, label='GPC (full)')
    gn = pd.read_csv(os.path.join(CONSDIR, 'Gn_deltaH.csv'))
    gc = pd.read_csv(os.path.join(CONSDIR, 'Gc_deltaH.csv'))
    # Map Gn/Gc back to GPC positions for overlay
    gn_start = 19
    gc_start = 652
    ax.plot(gn['pos'] + gn_start - 1, gn['deltaH'], color='blue', linewidth=1.0, label='Gn')
    ax.plot(gc['pos'] + gc_start - 1, gc['deltaH'], color='red', linewidth=1.0, label='Gc')
    ax.set_title('GPC conservation (ΔH): Gn and Gc along precursor')
    ax.set_xlabel('GPC position')
    ax.set_ylabel('ΔH (bits)')
    ax.set_ylim(0, 4.4)
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.3)
    out = os.path.join(PLOTDIR, 'GPC_overlay_deltaH.png')
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)
except Exception as e:
    pass

print(PLOTDIR)
