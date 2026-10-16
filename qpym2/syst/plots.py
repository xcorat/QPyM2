"""Plot helpers for the syst pipeline.

These mirror the figure styles already used in the syst notebooks
(``06-bmmc_ndbd.ipynb``, ``06-syst-rate-new.ipynb``) so the refactored
pipeline produces the same look-and-feel.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from qpym2.utils.mathy import get_mode
from qpym2.plotting.syst_posterior import plot_posterior_comparison


def plot_syst_comparison(
    nominal_bc: np.ndarray,
    syst_post: np.ndarray,
    delta: float | None = None,
    deltas: list = [],
    ax=None,
    bins: int = 60,
    syst_label: str = 'syst',
    nom_label: str = 'nominal',
    xlabel: str = r'$n_{0\nu\beta\beta}$',
    title: str = '',
    interval: str = 'hdi',
):
    plt.rcParams.update({
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 16,
        'legend.title_fontsize': 18,
        'font.size': 18,
        'font.family': 'serif',
        'font.serif': ['Computer Modern Roman', 'DejaVu Serif'],
        'font.sans-serif': ['Computer Modern Sans Serif', 'DejaVu Sans'],
        'mathtext.fontset': 'cm',
        'mathtext.default': 'regular',
        'axes.unicode_minus': True,
        'text.usetex': False,
    })
    """Overlay nominal vs systematic bias-corrected posteriors.

    Thin wrapper around :func:`qpym2.plotting.syst_posterior.plot_posterior_comparison`
    that adds the delta value to the title when supplied.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5))
    if delta is not None and title:
        title = f'{title}  (Δmode = {delta:+.1f})'
    if len(deltas) and title:
        title  = f'{title}\n $\Delta$[left, right, mode]  = [{deltas[0]:+.1f}, {deltas[1]:+.1f}, {deltas[2]:+.1f}]'
    elif delta is not None:
        title = f'Δmode = {delta:+.1f}'

    plot_posterior_comparison(
        syst_post, nominal_bc, ax=ax, bins=bins,
        syst_label=syst_label, nom_label=nom_label,
        xlabel=xlabel, title=title, interval=interval,
        ylabel=None
    )
    ax.legend() #fontsize=9)
    return ax.figure, ax


def plot_delta_table(
    syst_deltas,
    sigma_right: float,
    sigma_left: float,
    *,
    ax=None,
):
    """Tabular view of per-syst contributions and the quad-sum totals."""
    rows = [
        {'Source': d.name,
         '+σ (right)': float(d.delta_right),
         '-σ (left)': float(d.delta_left)}
        for d in syst_deltas
    ]
    rows.append({'Source': 'TOTAL (quad)',
                 '+σ (right)': float(sigma_right),
                 '-σ (left)': -float(sigma_left)})
    df = pd.DataFrame(rows)

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 0.4 * (len(rows) + 1) + 1))
    ax.axis('off')
    fmt = lambda v: f'{v:+.4f}' if isinstance(v, float) else str(v)
    cells = [[fmt(v) for v in row] for row in df.values]
    tbl = ax.table(cellText=cells, colLabels=df.columns,
                   loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.3)
    # bold the totals row
    for j in range(len(df.columns)):
        tbl[(len(df), j)].set_text_props(weight='bold')
    ax.set_title(
        f'Syst contributions   σ_right={sigma_right:.3f}   σ_left={sigma_left:.3f}',
        fontsize=10,
    )
    return ax.figure, ax, df


def plot_combined_posterior(
    marg_post: np.ndarray,
    nominal_bc: np.ndarray,
    sigma_right: float,
    sigma_left: float,
    *,
    kde_samples: np.ndarray | None = None,
):
    """Two-panel diagnostic: full posterior overlay + near-zero KDE zoom.

    Mirrors the figure in ``06-syst-rate-new.ipynb`` (cell with the
    side-by-side ``axes[0]`` / ``axes[1]`` plot).
    """
    marg_clip = np.clip(marg_post, 0, None)
    nom_clip = np.clip(nominal_bc, 0, None)
    nom_mode = get_mode(nom_clip)
    nom_p90 = float(np.percentile(nom_clip, 90))
    marg_mode = get_mode(marg_clip)
    marg_p90 = float(np.percentile(marg_clip, 90))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    xmax = max(np.percentile(marg_clip, 99.5), np.percentile(nom_clip, 99.5))
    bins = np.linspace(0, xmax, 100)
    ax.hist(nom_clip, bins=bins, density=True, alpha=0.35, color='tab:green',
            linewidth=0,
            label=f'nominal BC — mode={nom_mode:.2f}, p90={nom_p90:.2f}')
    ax.hist(marg_clip, bins=bins, density=True, histtype='step', lw=2,
            color='tab:red',
            label=f'marginalised — mode={marg_mode:.2f}, p90={marg_p90:.2f}')
    ax.axvline(nom_mode, color='tab:green', lw=1.2, ls='--')
    ax.axvline(marg_mode, color='tab:red', lw=1.2, ls='--')
    ax.axvline(nom_p90, color='tab:green', lw=1.0, ls=':')
    ax.axvline(marg_p90, color='tab:red', lw=1.0, ls=':')
    ax.set_xlabel('ndbd counts')
    ax.set_ylabel('density')
    ax.set_title('Nominal BC vs marginalised')
    # ax.legend(fontsize=8)

    ax2 = axes[1]
    zhi = np.percentile(marg_clip, 65)
    bins_z = np.linspace(-sigma_left * 3, zhi, 80)
    if kde_samples is not None:
        ax2.hist(kde_samples, bins=bins_z, density=True, alpha=0.3,
                 color='tab:blue', linewidth=0,
                 label='KDE samples (sub-zero OK)')
    ax2.hist(marg_post, bins=bins_z, density=True, histtype='step', lw=1.5,
             color='tab:orange', label='smeared (unclipped)')
    ax2.hist(marg_clip, bins=bins_z, density=True, histtype='step', lw=1.5,
             ls='--', color='tab:red', label='smeared (clipped)')
    ax2.axvline(0, color='k', lw=0.8, ls='--')
    ax2.set_xlabel('ndbd counts')
    ax2.set_ylabel('density')
    ax2.set_title('Near-zero zoom — KDE tail')
    ax2.legend(fontsize=8)

    fig.tight_layout()
    return fig, axes


def plot_rate_posterior(
    rate: np.ndarray,
    stats: dict,
    *,
    rate_no_psa: np.ndarray | None = None,
    stats_no_psa: dict | None = None,
    ax=None,
    title: str = 'Final rate posterior (marginalised + PSA)',
):
    """Plot a rate posterior with mode / p90 / HDI annotations.

    Mirrors the rate plot from ``06-syst-rate-new.ipynb``.  If a second
    (no-PSA) rate posterior is supplied, it is overlaid for comparison.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5))

    rate_xmax = float(np.percentile(rate, 99.8))
    ax.hist(rate, bins=120, range=(0, rate_xmax), density=True,
            color='#2b6ea3', alpha=0.25, linewidth=0,
            label=f'+ PSA — mode={stats["mode"]:.3e}, p90={stats["p90"]:.3e}')
    ax.axvline(stats['mode'], color='#c0392b', lw=1.6,
               label=f"mode = {stats['mode']:.3e} yr$^{{-1}}$")
    ax.axvline(stats['p90'], color='#27ae60', ls='--', lw=1.4,
               label=f"90% UL = {stats['p90']:.3e} yr$^{{-1}}$")
    lo, hi = stats['hdi']
    ax.axvspan(lo, hi, alpha=0.08, color='#2b6ea3',
               label=f"{int(stats['hdi_prob']*100)}% HDI [{lo:.2e}, {hi:.2e}]")

    if rate_no_psa is not None and stats_no_psa is not None:
        bins_cmp = np.linspace(0, rate_xmax, 120)
        ax.hist(rate_no_psa, bins=bins_cmp, density=True, histtype='step',
                lw=1.5, color='tab:green',
                label=(f'no PSA — mode={stats_no_psa["mode"]:.3e}, '
                       f'p90={stats_no_psa["p90"]:.3e}'))
        ax.axvline(stats_no_psa['p90'], color='tab:green', ls=':', lw=1.0)

    ax.set_xlabel(r'decay rate [yr$^{-1}$]')
    ax.set_ylabel('density')
    ax.set_title(title)
    # ax.legend(fontsize=8)
    ax.set_xlim(left=0)
    return ax.figure, ax


# ── Per-dataset HDI-edge barplots ─────────────────────────────────────────────

def _hatch_bars(bars, vals, *, keep_when_negative: bool):
    """Hatch + fade bars whose sign does not contribute to the quadrature sum.

    Parameters
    ----------
    keep_when_negative:
        True  → bar contributes only when value ≤ 0 (Δlo direction).
        False → bar contributes only when value ≥ 0 (Δhi direction).
    """
    for v, b in zip(vals, bars):
        non_contributing = (v > 0) if keep_when_negative else (v < 0)
        if non_contributing:
            b.set_hatch('//')
            b.set_alpha(0.35)


def _add_ds_separators(ax, x):
    """Add dotted grey vertical lines between dataset ticks."""
    for xi in x[:-1]:
        ax.axvline(xi + 0.5, color='grey', lw=0.6, ls=':', alpha=0.5)


def plot_lss_sscl_barplot(
    ds_ids,
    dlo: np.ndarray,
    dhi: np.ndarray,
    deltas: tuple,
    *,
    ax=None,
    co_shift: float = 1.0,
    figsize: tuple = (12, 5),
):
    """Single-panel barplot of per-dataset HDI-edge deltas for the sscl systematic.

    Parameters
    ----------
    ds_ids : sequence
        Dataset labels (plotted on x-axis).
    dlo : ndarray
        Per-ds left-edge shift  (Δlo = HDI_lo_syst − HDI_lo_nom).
    dhi : ndarray
        Per-ds right-edge shift (Δhi = HDI_hi_syst − HDI_hi_nom).
    deltas : (sigma_left, sigma_right)
        Quadrature-summed result from ``compute_1sigma_edge_delta``.
    ax : matplotlib Axes, optional
        Axes to draw on; a new figure is created if not supplied.
    co_shift : float
        Co-60 shift label shown in the title.
    figsize : tuple
        Figure size when a new figure is created.

    Returns
    -------
    fig, ax
    """
    from matplotlib.patches import Patch

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure

    x  = np.arange(len(ds_ids))
    bw = 0.35

    bars_lo = ax.bar(x - bw / 2, dlo, bw, color='C0',         alpha=0.85)
    bars_hi = ax.bar(x + bw / 2, dhi, bw, color='darkorange', alpha=0.85)

    _hatch_bars(bars_lo, dlo, keep_when_negative=True)
    _hatch_bars(bars_hi, dhi, keep_when_negative=False)

    ax.axhline(0, color='k', lw=0.8, ls='--')
    _add_ds_separators(ax, x)
    ax.set_xticks(x)
    ax.set_xticklabels([str(d) for d in ds_ids], rotation=45, ha='right', fontsize=8)
    ax.set_xlabel('Dataset')
    ax.set_ylabel(r'$\Delta$ HDI edge (counts)')
    ax.set_title(
        rf'sscl per-ds HDI-edge $\Delta$  (cs={co_shift})'
        rf'  — Δ=[{deltas[0]:+.3f}, {deltas[1]:+.3f}]'
    )
    legend_handles = [
        Patch(facecolor='C0',        alpha=0.85, label=r'$\Delta$ HDI-low (contributes when $\leq 0$)'),
        Patch(facecolor='darkorange', alpha=0.85, label=r'$\Delta$ HDI-high (contributes when $\geq 0$)'),
        Patch(facecolor='grey',       alpha=0.35, hatch='//', label='non-contributing (clipped to 0)'),
    ]
    ax.legend(handles=legend_handles, fontsize=8)
    return fig, ax


def plot_lss_dq_barplot(
    ds_ids,
    dlo_1r: np.ndarray,
    dhi_1r: np.ndarray,
    dlo_1l: np.ndarray,
    dhi_1l: np.ndarray,
    deltas: tuple,
    *,
    ax=None,
    co_shift: float = 1.0,
    figsize: tuple = (12, 5),
):
    """Single-panel barplot of per-dataset HDI-edge deltas for the dQ systematic.

    Four bars per dataset: ``[1r_lo | 1l_lo | 1r_hi | 1l_hi]``.
    Non-contributing bars (wrong sign) are hatched + faded.

    Parameters
    ----------
    ds_ids : sequence
        Dataset labels.
    dlo_1r, dhi_1r : ndarray
        Per-ds Δlo / Δhi for the ``dqres_1r`` variation.
    dlo_1l, dhi_1l : ndarray
        Per-ds Δlo / Δhi for the ``dqres_1l`` variation.
    deltas : (sigma_left, sigma_right)
        Quadrature-summed result from ``compute_1sigma_edge_delta``.
    ax : matplotlib Axes, optional
    co_shift : float
    figsize : tuple

    Returns
    -------
    fig, ax
    """
    from matplotlib.patches import Patch

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure

    x  = np.arange(len(ds_ids))
    bw = 0.2

    bars_1r_lo = ax.bar(x - 1.5*bw, dlo_1r, bw, color='C0',             alpha=0.85)
    bars_1l_lo = ax.bar(x - 0.5*bw, dlo_1l, bw, color='mediumseagreen', alpha=0.85)
    bars_1r_hi = ax.bar(x + 0.5*bw, dhi_1r, bw, color='darkorange',     alpha=0.85)
    bars_1l_hi = ax.bar(x + 1.5*bw, dhi_1l, bw, color='C3',             alpha=0.85)

    _hatch_bars(bars_1r_lo, dlo_1r, keep_when_negative=True)
    _hatch_bars(bars_1l_lo, dlo_1l, keep_when_negative=True)
    _hatch_bars(bars_1r_hi, dhi_1r, keep_when_negative=False)
    _hatch_bars(bars_1l_hi, dhi_1l, keep_when_negative=False)

    ax.axhline(0, color='k', lw=0.8, ls='--')
    _add_ds_separators(ax, x)
    ax.set_xticks(x)
    ax.set_xticklabels([str(d) for d in ds_ids], rotation=45, ha='right', fontsize=8)
    ax.set_xlabel('Dataset')
    ax.set_ylabel(r'$\Delta$ HDI edge (counts)')
    ax.set_title(
        rf'dQ per-ds HDI-edge $\Delta$  (cs={co_shift})'
        rf'  — Δ=[{deltas[0]:+.3f}, {deltas[1]:+.3f}]'
    )
    legend_handles = [
        Patch(facecolor='C0',             alpha=0.85, label=r'dQ$_{1r}$ $\Delta$lo'),
        Patch(facecolor='mediumseagreen',  alpha=0.85, label=r'dQ$_{1l}$ $\Delta$lo'),
        Patch(facecolor='darkorange',      alpha=0.85, label=r'dQ$_{1r}$ $\Delta$hi'),
        Patch(facecolor='C3',              alpha=0.85, label=r'dQ$_{1l}$ $\Delta$hi'),
        Patch(facecolor='grey',            alpha=0.35, hatch='//', label='non-contributing (clipped to 0)'),
    ]
    ax.legend(handles=legend_handles, fontsize=8, loc='upper right')
    return fig, ax
