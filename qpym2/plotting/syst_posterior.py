"""Plotting functions for systematic study posterior comparisons."""

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats as spstats
import matplotlib.patches as mpatches

from qpym2.utils.mathy import get_mode, calc_intervals
from qpym2.plotting.posterior_plots import plot_1d_posterior


def _boundary_kde(samples, y):
    """KDE with reflection at 0 to correct boundary bias for non-negative posteriors."""
    reflected = np.concatenate([samples, -samples])
    return spstats.gaussian_kde(reflected, bw_method='scott')(y) * 2


def plot_prior_posterior_violin(
    nom_trace, alt_trace, nom_comps, var_names,
    nom_label='nominal', alt_label='alternative',
    title=None, vwidth=0.38,
):
    """Split-violin prior-posterior comparison for all fit components.

    Left half-violin = nominal posterior (blue), right = alternative (red).
    Dashed navy line = prior shape (truncated normal or flat vertical line).
    Dotted horizontal lines = posterior modes.

    Parameters
    ----------
    nom_trace : az.InferenceData
        Nominal fit trace.
    alt_trace : az.InferenceData
        Alternative / systematic fit trace.
    nom_comps : pd.DataFrame
        Nominal fit components table (must have 'prior' column with
        ``(ptype, pparams)`` tuples or ``'flat'``).
    var_names : list of str
        Component names to plot (must exist in both traces and nom_comps.index).
    nom_label, alt_label : str
        Legend labels for the two posteriors.
    title : str, optional
        Figure suptitle.
    vwidth : float
        Max half-width of each violin side (in x-axis units).

    Returns
    -------
    fig, axes : matplotlib Figure and array of Axes.
    """
    n = len(var_names)
    fig, axes = plt.subplots(1, n, figsize=(2.8 * n, 5.5), squeeze=False)
    axes = axes[0]

    for ax, var in zip(axes, var_names):
        s_nom = nom_trace.posterior[var].values.flatten()
        s_alt = alt_trace.posterior[var].values.flatten()

        ylo = 0.0
        yhi = max(np.percentile(s_nom, 99.8), np.percentile(s_alt, 99.8))
        y = np.linspace(ylo, yhi, 400)

        # Boundary-corrected posterior KDEs
        p_nom = _boundary_kde(s_nom, y)
        p_nom = p_nom / p_nom.max() * vwidth
        p_alt = _boundary_kde(s_alt, y)
        p_alt = p_alt / p_alt.max() * vwidth

        # Half-violins
        ax.fill_betweenx(y, -p_nom, 0, alpha=0.45, color='tab:blue')
        ax.fill_betweenx(y, 0, p_alt, alpha=0.45, color='tab:red')
        ax.plot(-p_nom, y, color='tab:blue', lw=1.2)
        ax.plot(p_alt, y, color='tab:red', lw=1.2)

        # Prior (from nominal fit comps)
        pr = nom_comps.at[var, 'prior']
        ptype, pparams = pr if isinstance(pr, tuple) else ('flat', {})
        if ptype == 'flat':
            ax.plot([-vwidth, -vwidth], [ylo, yhi],
                    color='navy', lw=1.5, ls='--', alpha=0.8)
        else:
            mu = pparams.get('mu', 0.0)
            sigma = pparams.get('sigma', 1.0)
            a_bound = (0.0 - mu) / sigma
            prior_pdf = spstats.truncnorm.pdf(
                y, a=a_bound, b=np.inf, loc=mu, scale=sigma)
            prior_norm = (prior_pdf / prior_pdf.max() * vwidth
                          if prior_pdf.max() > 0 else prior_pdf)
            ax.plot(-prior_norm, y, color='navy', lw=1.5, ls='--', alpha=0.8)

        # Mode lines
        ax.axhline(get_mode(s_nom), color='tab:blue', ls=':', lw=1.4)
        ax.axhline(get_mode(s_alt), color='tab:red', ls=':', lw=1.4)
        ax.axvline(0, color='k', lw=0.6)

        ax.set_title(var, fontsize=10)
        ax.set_xlim(-vwidth * 1.25, vwidth * 1.25)
        ax.set_xticks([])
        ax.set_xlabel('')

    axes[0].set_ylabel('counts')
    for ax in axes[1:]:
        ax.set_ylabel('')

    legend_handles = [
        mpatches.Patch(color='tab:blue', alpha=0.6, label=f'{nom_label} posterior'),
        mpatches.Patch(color='tab:red', alpha=0.6, label=f'{alt_label} posterior'),
        plt.Line2D([0], [0], color='navy', ls='--', lw=1.5, label='prior'),
        plt.Line2D([0], [0], color='gray', ls=':', lw=1.4, label='mode'),
    ]
    fig.legend(handles=legend_handles, loc='upper right', fontsize=9,
               bbox_to_anchor=(1.0, 1.0), frameon=True)
    if title:
        fig.suptitle(title, fontsize=12)
    plt.tight_layout()

    return fig, axes


def plot_posterior_comparison(syst_samples, nom_samples, ax, bins=60,
                               syst_color='steelblue', nom_color='C3',
                               syst_label='syst', nom_label='nominal',
                               xlabel='ndbd (events)', ylabel='density', title='',
                               interval='hdi', **kwargs):
    """Overlay syst (filled) and nominal (line-only) posteriors with mode + 1σ interval markers.

    Parameters
    ----------
    syst_samples : array-like
        Systematic posterior samples (drawn with filled histogram).
    nom_samples : array-like
        Nominal posterior samples (step line only, no fill).
    ax : matplotlib.axes.Axes
    bins : int, optional
    syst_color, nom_color : str, optional
    syst_label, nom_label : str, optional
    xlabel, ylabel, title : str, optional
    interval : {'hdi', 'half_interval', 'percentiles', 'stddev'}, optional
        Method used to compute the 1σ interval bounds (default ``'hdi'``).
        Passed to :func:`qpym2.utils.mathy.calc_intervals`.
    """
    syst_samples = np.asarray(syst_samples).flatten()
    nom_samples  = np.asarray(nom_samples).flatten()

    (s_mode, s_low, s_high), = calc_intervals([syst_samples], method=interval)
    (n_mode, n_low, n_high), = calc_intervals([nom_samples],  method=interval)

    print(ylabel)
    # syst: step histogram + semi-transparent fill
    plot_1d_posterior(syst_samples, ax=ax, bins=bins, color=syst_color, lw=1.2,
                      label=syst_label,
                      lines=[
                          dict(x=s_mode, color=syst_color, ls='--', lw=1.5,
                               label=f'{syst_label}: mode={s_mode:.2f}'),
                          dict(x=s_low,  color=syst_color, ls=':',  lw=1.0,
                               label=f'{syst_label}: [{interval}]=[{s_low:.2f}, {s_high:.2f}]'),
                          dict(x=s_high, color=syst_color, ls=':',  lw=1.0),
                      ],
                      xlabel=xlabel, ylabel=ylabel, title=title, **kwargs)
    h, xe = np.histogram(syst_samples, bins=bins, density=True)
    ax.fill_between(0.5 * (xe[:-1] + xe[1:]), h, alpha=0.2, color=syst_color, step='mid')

    # nominal: step line only (no fill)
    plot_1d_posterior(nom_samples, ax=ax, bins=bins, color=nom_color, lw=1.5,
                      label=nom_label,
                      lines=[
                          dict(x=n_mode, color=nom_color, ls='--', lw=1.5,
                               label=f'{nom_label}: mode={n_mode:.2f}'),
                          dict(x=n_low,  color=nom_color, ls=':',  lw=1.0,
                               label=f'{nom_label}: [{interval}]=[{n_low:.2f}, {n_high:.2f}]'),
                          dict(x=n_high, color=nom_color, ls=':',  lw=1.0),
                      ], 
                      xlabel=xlabel, ylabel=ylabel, title=title,
                        **kwargs)


def plot_multi_posterior_overlay(
    samples_list, labels, ax,
    bins=80, palette='tab10',
    xlabel='counts', ylabel='density', title='',
    highlight_idx=None, highlight_lw=2.2, default_lw=1.2,
    show_mode=True, show_p16p84=False,
):
    """Overlay multiple posteriors as step histograms on a shared axis.

    Each posterior is drawn as a step line using a distinct colour from
    ``palette``.  Optionally one entry can be highlighted with a thicker line.
    Mode (and optionally p16/p84) markers are drawn as vertical dotted lines
    in the matching colour.

    Parameters
    ----------
    samples_list : list of array-like
        One 1-D array of posterior samples per series.
    labels : list of str
        Legend label for each series (same length as ``samples_list``).
    ax : matplotlib.axes.Axes
    bins : int or array-like, optional
        Histogram bins shared across all series.
    palette : str, optional
        Seaborn / matplotlib colour palette name (default ``'tab10'``).
    xlabel, ylabel, title : str, optional
    highlight_idx : int or None, optional
        Index of the series to draw with a thicker line.
    highlight_lw, default_lw : float, optional
        Line widths for the highlighted and regular series.
    show_mode : bool, optional
        Draw a vertical dashed line at each posterior mode (default ``True``).
    show_p16p84 : bool, optional
        Also draw dotted verticals at p16 and p84 for each series.
    """
    try:
        from matplotlib import colormaps
        get_cmap = colormaps.get_cmap
    except (ImportError, AttributeError):
        from matplotlib.cm import get_cmap
    import matplotlib.colors as mcolors

    try:
        import seaborn as _sns
        colors = _sns.color_palette(palette, n_colors=len(samples_list))
    except (ImportError, Exception):
        cmap = get_cmap(palette)
        colors = [cmap(i / max(len(samples_list) - 1, 1))
                  for i in range(len(samples_list))]

    # Build shared bin edges from the pooled range
    if isinstance(bins, int):
        all_concat = np.concatenate([np.asarray(s).flatten() for s in samples_list])
        xmax = np.percentile(all_concat, 99.5)
        bins = np.linspace(0.0, xmax, bins + 1)

    for i, (samples, label, color) in enumerate(zip(samples_list, labels, colors)):
        samples = np.asarray(samples).flatten()
        lw = highlight_lw if i == highlight_idx else default_lw
        h, xe = np.histogram(samples, bins=bins, density=True)
        xc = 0.5 * (xe[:-1] + xe[1:])
        ax.step(xc, h, where='mid', color=color, lw=lw, label=label)

        if show_mode:
            mode_val = get_mode(samples)
            ax.axvline(mode_val, color=color, ls='--', lw=0.9, alpha=0.7)

        if show_p16p84:
            p16_val, p84_val = np.percentile(samples, [16, 84])
            ax.axvline(p16_val, color=color, ls=':', lw=0.8, alpha=0.6)
            ax.axvline(p84_val, color=color, ls=':', lw=0.8, alpha=0.6)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
