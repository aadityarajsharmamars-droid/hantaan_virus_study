#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path


def load_ic50_values(path: Path,
                      type_col: str = 'standard_type',
                      value_col: str = 'standard_value',
                      units_col: str = 'standard_units',
                      relation_col: str = 'standard_relation',
                      ic50_label: str = 'IC50',
                      allowed_units=('nM','nm','nM '),
                      require_equality=True):
    df = pd.read_csv(path)
    # Normalize column names
    df.columns = [c.strip() for c in df.columns]
    # Filter to IC50 rows
    m = df[type_col].astype(str).str.upper().str.strip() == ic50_label.upper()
    if units_col in df.columns:
        m &= df[units_col].astype(str).str.lower().str.strip().isin([u.lower() for u in allowed_units])
    if require_equality and relation_col in df.columns:
        m &= df[relation_col].astype(str).str.strip() == '='
    sub = df.loc[m].copy()
    # Coerce numeric
    sub[value_col] = pd.to_numeric(sub[value_col], errors='coerce')
    vals = sub[value_col].dropna().values.astype(float)
    # Remove non-positive
    vals = vals[vals > 0]
    return vals, sub


def main():
    p = argparse.ArgumentParser(description='Analyze IC50 distribution and fit normal in log10-space.')
    p.add_argument('--input', type=Path, default=Path('data/summary/pocket1_ic50_activities.csv'), help='Input CSV file with IC50 data')
    p.add_argument('--outdir', type=Path, default=Path('plots'), help='Output directory for plots')
    p.add_argument('--stats_out', type=Path, default=Path('tables/ic50_distribution_stats.csv'), help='Output CSV for summary stats')
    p.add_argument('--units', type=str, default='nM', help='Units of IC50 values to assume/filter (default: nM)')
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.stats_out.parent.mkdir(parents=True, exist_ok=True)

    vals, sub = load_ic50_values(args.input, allowed_units=(args.units,))
    if vals.size == 0:
        raise SystemExit('No IC50 values found after filtering. Check column names/units/relations.')

    # Log10-transform (nM) to approximate normal
    log_ic50 = np.log10(vals)
    mu, sigma = np.mean(log_ic50), np.std(log_ic50, ddof=1)

    # Shapiro-Wilk normality test
    W, pval = stats.shapiro(log_ic50) if log_ic50.size >= 3 else (np.nan, np.nan)

    # Back-transform summary: geometric mean (GM) and geometric SD (GSD)
    GM = 10**mu
    GSD = 10**sigma

    # 95% CI for the mean in log-space, then back-transform to fold-range around GM
    se = sigma/np.sqrt(log_ic50.size)
    ci95_mu = stats.t.interval(0.95, df=log_ic50.size-1, loc=mu, scale=se) if log_ic50.size >= 2 else (np.nan, np.nan)
    ci95_low, ci95_high = (10**ci95_mu[0], 10**ci95_mu[1]) if all(np.isfinite(ci95_mu)) else (np.nan, np.nan)

    # Histogram with fitted normal PDF in log-space
    fig, ax = plt.subplots(figsize=(6,4))
    n, bins, patches = ax.hist(log_ic50, bins=max(5, int(np.sqrt(log_ic50.size))), density=True, alpha=0.6, color='C0', edgecolor='black')
    x = np.linspace(min(bins), max(bins), 200)
    ax.plot(x, stats.norm.pdf(x, mu, sigma), 'r-', lw=2, label=f'Normal fit (mu={mu:.2f}, sigma={sigma:.2f})')
    ax.set_xlabel('log10(IC50 [nM])')
    ax.set_ylabel('Density')
    ax.set_title('IC50 distribution (log10-space)')
    ax.legend()
    fig.tight_layout()
    hist_path_png = args.outdir / 'ic50_log10nM_hist.png'
    hist_path_pdf = args.outdir / 'ic50_log10nM_hist.pdf'
    fig.savefig(hist_path_png, dpi=200)
    fig.savefig(hist_path_pdf)
    plt.close(fig)

    # QQ plot
    fig2, ax2 = plt.subplots(figsize=(5,5))
    stats.probplot(log_ic50, dist="norm", sparams=(mu, sigma), plot=ax2)
    ax2.set_title('Q-Q plot of log10(IC50 [nM]) vs Normal')
    fig2.tight_layout()
    qq_path = args.outdir / 'ic50_log10nM_qq.png'
    fig2.savefig(qq_path, dpi=200)
    plt.close(fig2)

    # Linear-scale histogram (optional)
    fig3, ax3 = plt.subplots(figsize=(6,4))
    ax3.hist(vals, bins=max(5, int(np.sqrt(vals.size))), color='C1', edgecolor='black')
    ax3.set_xlabel(f'IC50 [{args.units}]')
    ax3.set_ylabel('Count')
    ax3.set_title('IC50 distribution (linear scale)')
    fig3.tight_layout()
    lin_hist_path = args.outdir / 'ic50_linear_hist.png'
    fig3.savefig(lin_hist_path, dpi=200)
    plt.close(fig3)

    # Save stats
    stats_df = pd.DataFrame({
        'n': [log_ic50.size],
        'mu_log10_nM': [mu],
        'sigma_log10_nM': [sigma],
        'shapiro_W': [W],
        'shapiro_p': [pval],
        'GM_nM': [GM],
        'GSD_fold': [GSD],
        'GM_95CI_low_nM': [ci95_low],
        'GM_95CI_high_nM': [ci95_high],
    })
    stats_df.to_csv(args.stats_out, index=False)

    print(f'Processed {log_ic50.size} IC50 values from {args.input}')
    print(f'- Saved histogram (log10-space) to {hist_path_png}')
    print(f'- Saved Q-Q plot to {qq_path}')
    print(f'- Saved histogram (linear) to {lin_hist_path}')
    print(f'- Saved summary stats to {args.stats_out}')


if __name__ == '__main__':
    main()
