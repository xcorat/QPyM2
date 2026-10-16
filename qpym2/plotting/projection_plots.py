"""1D projection plots with pull panels for fit validation.

Provides numpy-based plotting functions for comparing data and fit spectra
derived from 2D histograms in u-v (rotated energy) space.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from qpym2.utils.hists import get_h2proj


def plot_proj_with_pull(h1_data, h1_fit, left_edges, xlabel, xlim=None, title=None, figsize=(12, 8), rebin=1, log_scale=True):
    """2-panel plot: 1D spectrum (top) + pull panel (bottom).

    Args:
        h1_data:     1D array of data counts per bin
        h1_fit:      1D array of fit counts per bin
        left_edges:  left edge of each bin (keV)  — length == len(h1_data)
        xlabel:      x-axis label string
        xlim:        (xmin, xmax) tuple, or None for auto scaling
        title:       optional title string
        figsize:     figure size tuple
        rebin:       integer — merge this many consecutive bins (default 1)
        log_scale:   boolean — use log scale on y-axis (default True)
    Returns:
        (fig, ax_spec, ax_pull)
    """
    if rebin > 1:
        n          = (len(h1_data) // rebin) * rebin
        h1_data    = h1_data[:n].reshape(-1, rebin).sum(axis=1)
        h1_fit     = h1_fit[:n].reshape(-1, rebin).sum(axis=1)
        left_edges = left_edges[:n:rebin]

    bin_w    = float(np.median(np.diff(left_edges)))
    centers  = left_edges + bin_w / 2
    pull     = (h1_data - h1_fit) / np.sqrt(np.maximum(h1_fit, 1))

    fig = plt.figure(figsize=figsize)
    gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax_spec = fig.add_subplot(gs[0])
    ax_pull = fig.add_subplot(gs[1], sharex=ax_spec)

    # --- Spectrum ---
    bar_color = '#4e8ab5'
    yerr_data = np.maximum(np.sqrt(h1_data), 1)
    ax_spec.bar(left_edges, h1_data, width=bin_w, align='edge',
                color=bar_color, alpha=0.6, linewidth=0, label='Data')
    ax_spec.errorbar(centers, h1_data, yerr=yerr_data,
                     fmt='none', color=bar_color, linewidth=1.0, capsize=2)
    ax_spec.plot(centers, h1_fit, 'o-',
                 color='C3', linewidth=1.8, markersize=3.5, label='Fit (mode)')
    ax_spec.set_ylabel('Counts / bin')
    if log_scale:
        ax_spec.set_yscale('log')
    if xlim is not None:
        ax_spec.set_xlim(*xlim)
    ax_spec.legend(frameon=False)
    ax_spec.minorticks_on()
    plt.setp(ax_spec.get_xticklabels(), visible=False)
    if title:
        ax_spec.set_title(title)

    # --- Pull ---
    ax_pull.bar(left_edges, pull, width=bin_w, align='edge',
                color='grey', alpha=0.7, linewidth=0)
    ax_pull.axhline(0,  color='black', linewidth=0.8)
    ax_pull.axhline( 1, color='C2', linewidth=0.8, linestyle='--', label=r'$\pm$1')
    ax_pull.axhline(-1, color='C2', linewidth=0.8, linestyle='--')
    ax_pull.axhline( 2, color='C1', linewidth=0.8, linestyle=':', label=r'$\pm$2')
    ax_pull.axhline(-2, color='C1', linewidth=0.8, linestyle=':')
    ax_pull.set_xlabel(xlabel)
    ax_pull.set_ylabel(r'$(d - f)\,/\,\sqrt{f}$')
    ax_pull.set_ylim(-4, 4)
    ax_pull.minorticks_on()
    # ax_pull.legend(frameon=False, fontsize=8, loc='upper right')
    # ax_pull.legend(frameon=False, fontsize=18, loc='lower right')
    # ax.legend(fontsize=18, loc='lower right') 
    plt.tight_layout()
    return fig, ax_spec, ax_pull


def plot_proj_h2(h2_data, h2_fit, xedges, yedges, hm, axis='esum', rebin=1, title=None, figsize=(12, 8), xlim=None, log_scale=True):
    """Project 2D histograms and draw spectrum + pull panel.

    Args:
        h2_data:   2D data histogram, shape (nu, nv)
        h2_fit:    2D fit histogram, shape (nu, nv)
        xedges:    u-axis bin edges (length nu+1)
        yedges:    v-axis bin edges (length nv+1)
        hm:        Histogram model (TempConfig)
        axis:      'esum' or 'ediff' — which axis to project
        rebin:     merge this many consecutive bins
        title:     plot title (auto-generated if None)
        figsize:   figure size tuple
        xlim:      optional tuple (xmin, xmax) to override default xlim
        log_scale: boolean — use log scale on y-axis (default True)
    Returns:
        (fig, ax_spec, ax_pull)
    """
    if axis == 'esum':
        h1_data    = get_h2proj(h2_data, axis=1)
        h1_fit     = get_h2proj(h2_fit,  axis=1)
        left_edges = xedges[:-1] * np.sqrt(2)
        xlim       = xlim or (hm.esum_min, hm.esum_max)
        xlabel     = r'E$_\mathrm{sum}$ [keV]'
        title      = title or 'ESum projection'
    elif axis == 'ediff':
        h1_data    = get_h2proj(h2_data, axis=0)
        h1_fit     = get_h2proj(h2_fit,  axis=0)
        left_edges = yedges[:-1] * np.sqrt(2)
        # Derive from hm: ediff_max = vmax * √2 = (esum_max - 2*threshold)
        ediff_max = (hm.esum_max - 2 * hm.threshold)
        xlim       = xlim or (0, ediff_max)
        xlabel     = r'E$_\mathrm{diff}$ [keV]'
        title      = title or 'EDiff projection'
    else:
        raise ValueError(f"axis must be 'esum' or 'ediff', got {axis!r}")

    return plot_proj_with_pull(h1_data, h1_fit, left_edges, xlabel,
                               xlim=xlim, title=title, figsize=figsize, rebin=rebin, log_scale=log_scale)
