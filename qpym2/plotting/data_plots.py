"""Plotting functions for Chapter 2 data visualizations.

Produces:
  - 2D ROI heatmap in u-v coordinates
  - Full expected spectra in (E1, E2) coordinates
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from qpym2.utils.hists import get_hist_settings


# Default colors for the four merged background groups
COMP_COLORS = {
    'co60':   '#1f77b4',
    'u238':   '#ff7f0e',
    'th232':  '#2ca02c',
    'others': '#9467bd',
    'ndbd':   '#d62728',
}


def plot_2d_heatmap(hist_2d, hm, ax=None, title='', cmap='viridis',
                    log_scale=True, show_threshold=True, colorbar=True):
    """Plot a 2D histogram as a heatmap in u-v coordinates.

    Args:
        hist_2d: 2D numpy array (nu × nv bins).
        hm: TempConfig histogram model with umin/umax/vmin/vmax.
        ax: matplotlib Axes (created if None).
        title: plot title.
        cmap: colormap name.
        log_scale: use logarithmic colour scale.
        show_threshold: draw a horizontal line at the threshold v-cut.
        colorbar: add a colorbar.

    Returns:
        matplotlib Axes.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))

    nbins, ranges = get_hist_settings(hm)
    u_edges = np.linspace(ranges[0][0], ranges[0][1], nbins[0] + 1)
    v_edges = np.linspace(ranges[1][0], ranges[1][1], nbins[1] + 1)

    plot_data = hist_2d.T  # pcolormesh expects (v, u) indexing
    norm = mcolors.LogNorm(vmin=max(plot_data[plot_data > 0].min(), 1e-2),
                           vmax=plot_data.max()) if log_scale else None

    mesh = ax.pcolormesh(u_edges, v_edges, plot_data, cmap=cmap, norm=norm)
    if colorbar:
        plt.colorbar(mesh, ax=ax, label='Counts')

    ax.set_xlabel(r'$u = E_{\mathrm{sum}}/\sqrt{2}$ [keV/$\sqrt{2}$]')
    ax.set_ylabel(r'$v = \Delta E/\sqrt{2}$ [keV/$\sqrt{2}$]')
    if title:
        ax.set_title(title)

    if show_threshold:
        threshold_v = (hm.esum_max - 2 * hm.threshold) / math.sqrt(2)
        ax.axhline(y=threshold_v, color='white', ls='--', lw=1.2,
                   label=f'threshold v={threshold_v:.0f}')

    return ax


def plot_expected_spectra_e1e2(hist_2d, hm, ax=None, title='Expected M2 Spectra',
                                cmap='viridis', log_scale=True, colorbar=True):
    """Plot the expected 2D spectra in rotated-back (E1, E2) coordinates.

    Takes a u-v histogram (e.g. summed MC background + signal) and displays
    it in the original (E1, E2) energy space using a non-rectangular mesh.

    Args:
        hist_2d: 2D numpy array (nu × nv bins) in u-v space.
        hm: TempConfig histogram model.
        ax: matplotlib Axes (created if None).
        title: plot title.
        cmap: colormap.
        log_scale: logarithmic colour normalisation.
        colorbar: show colorbar.

    Returns:
        matplotlib Axes.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))

    nbins, ranges = get_hist_settings(hm)
    u_edges = np.linspace(ranges[0][0], ranges[0][1], nbins[0] + 1)
    v_edges = np.linspace(ranges[1][0], ranges[1][1], nbins[1] + 1)

    # Build mesh grid in u-v, then rotate to e1-e2
    U, V = np.meshgrid(u_edges, v_edges, indexing='ij')
    s2 = math.sqrt(2)
    E1 = (U + V) / s2
    E2 = (U - V) / s2

    plot_data = hist_2d.copy().astype(float)
    plot_data[plot_data <= 0] = np.nan  # hide empty bins

    norm = None
    if log_scale:
        valid = plot_data[np.isfinite(plot_data) & (plot_data > 0)]
        if len(valid):
            norm = mcolors.LogNorm(vmin=max(valid.min(), 1e-2), vmax=valid.max())

    mesh = ax.pcolormesh(E1, E2, plot_data, cmap=cmap, norm=norm)
    if colorbar:
        plt.colorbar(mesh, ax=ax, label='Counts')

    ax.set_xlabel(r'$E_1$ [keV]')
    ax.set_ylabel(r'$E_2$ [keV]')
    if title:
        ax.set_title(title)
    ax.set_aspect('equal')

    return ax


def make_total_expected_hist(mctable, signal_table=None, comp_cfg=None, smooth=None):
    """Build total expected 2D histogram from MC components.

    Sums background components (weighted by BM mode) and optionally adds
    a signal histogram (unweighted, i.e. raw integral).

    Args:
        mctable: DataFrame from hdf.read_bkg_model() with 'mode'/'mean'
                 and 'mchist' columns.
        signal_table: DataFrame from hdf.read_signal() (optional).
        comp_cfg: component grouping config list (optional).  If provided,
                  groups are merged first via make_bkgcomps_table().
        smooth: smoothing tuple or False.

    Returns:
        2D numpy array of expected counts.
    """
    if comp_cfg is not None:
        from qpym2.models.m2 import make_bkgcomps_table
        comps = make_bkgcomps_table(comp_cfg, mctable, smooth=smooth)
        total = np.sum([row['integral'] * row['pdf_hist']
                        for _, row in comps.iterrows()], axis=0)
    else:
        norm_col = 'mode' if 'mode' in mctable.columns else 'mean'
        total = np.sum(
            [row[norm_col] * row['mchist'] / np.sum(row['mchist'])
             for _, row in mctable.iterrows()
             if np.sum(row['mchist']) > 0],
            axis=0
        )

    if signal_table is not None:
        sig_hist = signal_table['mchist'].iloc[0]
        total = total + sig_hist

    return total
