#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
from typing import List, Optional

CSV_PIC50_COLS = [
    'pIC50', 'pic50', 'pIC50_value', 'pIC50_M',
    'pchembl_value', 'pchembl', 'pChEMBL', 'pchembl_max', 'influenza_pchembl'
]
CSV_IC50_COLS = [
    'IC50', 'ic50', 'ic50_nM', 'standard_value', 'influenza_ic50_nM'
]

SDF_PIC50_PROPS = [
    'pIC50', 'pic50', 'pIC50_value', 'pIC50_M',
    'pchembl_value', 'pchembl', 'pChEMBL', 'pchembl_max'
]
SDF_IC50_PROPS = ['IC50', 'ic50', 'standard_value']
SDF_UNITS_PROPS = ['units', 'standard_units']
SDF_REL_PROPS = ['relation', 'standard_relation']


def _select_first_present(cols: List[str], available: List[str]) -> Optional[str]:
    for c in cols:
        if c in available:
            return c
    return None


def load_from_csv(path: Path, sample_n: Optional[int], only_type: Optional[str]=None, only_units: Optional[str]=None, require_relation: Optional[str]=None):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]

    # Optional filtering by type/units/relation if columns exist
    m = pd.Series(True, index=df.index)
    if only_type is not None and 'standard_type' in df.columns:
        m &= df['standard_type'].astype(str).str.upper().str.strip() == only_type.upper()
    if only_units is not None and 'standard_units' in df.columns:
        m &= df['standard_units'].astype(str).str.lower().str.strip() == only_units.lower()
    if require_relation is not None and 'standard_relation' in df.columns:
        m &= df['standard_relation'].astype(str).str.strip() == require_relation
    df = df.loc[m].copy()

    pic50_col = _select_first_present(CSV_PIC50_COLS, df.columns.tolist())
    if pic50_col is not None:
        vals = pd.to_numeric(df[pic50_col], errors='coerce').dropna().values
    else:
        # Try IC50 in nM -> pIC50 = 9 - log10(IC50[nM])
        ic50_col = _select_first_present(CSV_IC50_COLS, df.columns.tolist())
        if ic50_col is None:
            raise ValueError('No pIC50/pChEMBL column or IC50 column found in CSV.')
        ic50 = pd.to_numeric(df[ic50_col], errors='coerce').dropna().values.astype(float)
        ic50 = ic50[ic50 > 0]
        vals = 9 - np.log10(ic50)

    if sample_n is not None and vals.size > sample_n:
        rng = np.random.default_rng(123)
        vals = rng.choice(vals, size=sample_n, replace=False)
    return vals


def load_from_sdf(path: Path, sample_n: Optional[int]):
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(str(path), removeHs=False)
    v = []
    for mol in suppl:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        # Normalize keys case-sensitive as stored; build lower map
        keys = list(props.keys())
        pic50_key = None
        for k in SDF_PIC50_PROPS:
            if k in props:
                pic50_key = k
                break
        if pic50_key is not None:
            try:
                v.append(float(str(props[pic50_key]).strip()))
                continue
            except Exception:
                pass
        # Try IC50 values in nM
        ic50_key = None
        for k in SDF_IC50_PROPS:
            if k in props:
                ic50_key = k
                break
        if ic50_key is not None:
            try:
                val = float(str(props[ic50_key]).strip())
            except Exception:
                val = None
            if val is not None and val > 0:
                # Check units: assume nM unless indicated otherwise
                units = None
                for uk in SDF_UNITS_PROPS:
                    if uk in props:
                        units = str(props[uk]).strip()
                        break
                if units is None or units.lower() in ['nm', 'nmol/l', 'nmolar', 'nM'.lower()]:
                    v.append(9 - np.log10(val))
                elif units.lower() in ['m', 'molar']:
                    # pIC50 directly
                    v.append(-np.log10(val))
                elif units.lower() in ['um', 'µm', 'uM'.lower()]:
                    v.append(6 - np.log10(val))
                elif units.lower() in ['pm']:
                    v.append(12 - np.log10(val))
                else:
                    # unknown units; skip
                    pass
    vals = np.array([x for x in v if np.isfinite(x)])
    if sample_n is not None and vals.size > sample_n:
        rng = np.random.default_rng(123)
        vals = rng.choice(vals, size=sample_n, replace=False)
    return vals


def analyze(vals: np.ndarray, outdir: Path, prefix: str):
    if vals.size < 3:
        raise SystemExit(f'Not enough pIC50 values ({vals.size}) to fit a normal distribution.')
    mu, sigma = np.mean(vals), np.std(vals, ddof=1)
    W, pval = stats.shapiro(vals) if vals.size <= 5000 else (np.nan, np.nan)

    # Histogram + fit
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6,4))
    ax.hist(vals, bins=max(10, int(np.sqrt(vals.size))), density=True, alpha=0.6, color='C0', edgecolor='black')
    x = np.linspace(vals.min()-0.5, vals.max()+0.5, 300)
    ax.plot(x, stats.norm.pdf(x, mu, sigma), 'r-', lw=2, label=f'Normal fit (mu={mu:.2f}, sigma={sigma:.2f})')
    ax.set_xlabel('pIC50 (-log10[M])')
    ax.set_ylabel('Density')
    ax.set_title('pIC50 distribution')
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / f'{prefix}_pic50_hist.png', dpi=200)
    fig.savefig(outdir / f'{prefix}_pic50_hist.pdf')
    plt.close(fig)

    # QQ plot
    fig2, ax2 = plt.subplots(figsize=(5,5))
    stats.probplot(vals, dist='norm', sparams=(mu, sigma), plot=ax2)
    ax2.set_title('Q-Q plot of pIC50 vs Normal')
    fig2.tight_layout()
    fig2.savefig(outdir / f'{prefix}_pic50_qq.png', dpi=200)
    plt.close(fig2)

    # Save stats
    stats_df = pd.DataFrame({
        'n':[vals.size],
        'mu':[mu],
        'sigma':[sigma],
        'shapiro_W':[W],
        'shapiro_p':[pval],
    })
    stats_df.to_csv(outdir / f'{prefix}_pic50_stats.csv', index=False)


def main():
    ap = argparse.ArgumentParser(description='Make normal distribution of pIC50 values from CSV or SDF.')
    ap.add_argument('--input', type=Path, required=True, help='Input CSV or SDF file containing pIC50/pChEMBL or IC50 values')
    ap.add_argument('--sample', type=int, default=None, help='Randomly sample this many molecules if more are available (e.g., 300)')
    ap.add_argument('--outdir', type=Path, default=Path('plots'))
    ap.add_argument('--prefix', type=str, default='dataset')
    ap.add_argument('--csv_only_type', type=str, default=None, help='Filter CSV rows to this standard_type (e.g., IC50)')
    ap.add_argument('--csv_only_units', type=str, default=None, help='Filter CSV rows to this standard_units (e.g., nM)')
    ap.add_argument('--csv_require_relation', type=str, default=None, help='Require this standard_relation (e.g., =)')
    args = ap.parse_args()

    ext = args.input.suffix.lower()
    if ext == '.csv':
        vals = load_from_csv(args.input, args.sample,
                             only_type=args.csv_only_type,
                             only_units=args.csv_only_units,
                             require_relation=args.csv_require_relation)
    elif ext == '.sdf':
        vals = load_from_sdf(args.input, args.sample)
    else:
        raise SystemExit('Unsupported input format. Use CSV or SDF.')
    analyze(vals, args.outdir, args.prefix)
    print(f'Processed {vals.size} pIC50 values from {args.input}')
    print(f'- Plots and stats saved to {args.outdir} with prefix {args.prefix}_pic50_*')


if __name__ == '__main__':
    main()
