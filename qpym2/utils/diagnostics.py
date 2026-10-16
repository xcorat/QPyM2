"""Fit diagnostics utilities.

Core functions for comparing fits via log-likelihood samples (pik).

In QPyM2 the ``pik`` CustomDist uses a scalar logp function (sum over all
events), so ``trace.log_likelihood.pik`` has shape ``(chains, draws)`` — one
total log-likelihood value per MCMC draw.  Flattening gives a 1-D array of
per-draw log-likelihood samples that acts as a goodness-of-fit score for the
fit (consistent with usage in ``qpym2.utils.syst_rwp``).

Comparing two fits therefore means comparing their distributions of total
log-likelihood across draws.  The mode of this distribution is used as the
primary quality metric, following the convention in ``syst_rwp.combine_posteriors``.
"""

import numpy as np


def get_pik_samples(trace):
    """Return per-draw total log-likelihood samples from a saved trace.

    Equivalent to ``trace.log_likelihood.pik.values.flatten()`` as used
    throughout the QPyM2 codebase (e.g. ``syst_rwp.py``, ``results.py``).

    Parameters
    ----------
    trace : arviz.InferenceData
        Trace containing a ``log_likelihood.pik`` group.

    Returns
    -------
    pik : np.ndarray, shape (n_chains * n_draws,)
        Per-draw total log-likelihood values.
    """
    return trace.log_likelihood.pik.values.flatten()


def compare_pik(trace_a, trace_b, label_a='fit_a', label_b='fit_b'):
    """Compare log-likelihood distributions between two fits.

    Computes summary statistics of the per-draw total log-likelihood for each
    fit and quantifies the shift between them.

    Parameters
    ----------
    trace_a, trace_b : arviz.InferenceData
        Two traces each containing ``log_likelihood.pik``.
    label_a, label_b : str
        Human-readable labels for printing / plotting.

    Returns
    -------
    result : dict with keys:
        ``pik_a``        – 1-D array of per-draw log-likelihoods for trace_a
        ``pik_b``        – 1-D array of per-draw log-likelihoods for trace_b
        ``mean_a``       – mean of pik_a
        ``mean_b``       – mean of pik_b
        ``delta_mean``   – mean_b - mean_a
        ``frac_diff``    – delta_mean / |mean_a|  (fractional shift)
        ``ks_stat``      – KS test statistic (scipy)
        ``ks_pvalue``    – KS test p-value
        ``label_a``      – label_a
        ``label_b``      – label_b
    """
    from scipy.stats import ks_2samp

    pik_a = get_pik_samples(trace_a)
    pik_b = get_pik_samples(trace_b)

    mean_a = float(pik_a.mean())
    mean_b = float(pik_b.mean())
    delta_mean = mean_b - mean_a
    frac_diff = delta_mean / abs(mean_a) if mean_a != 0.0 else float('nan')

    ks_stat, ks_pvalue = ks_2samp(pik_a, pik_b)

    return {
        'pik_a': pik_a,
        'pik_b': pik_b,
        'mean_a': mean_a,
        'mean_b': mean_b,
        'delta_mean': delta_mean,
        'frac_diff': frac_diff,
        'ks_stat': float(ks_stat),
        'ks_pvalue': float(ks_pvalue),
        'label_a': label_a,
        'label_b': label_b,
    }


def plot_pik_hist(result, ax=None, bins=80):
    """Overlay histogram of per-draw pik distributions for two fits.

    Parameters
    ----------
    result : dict
        Output of :func:`compare_pik`.
    ax : matplotlib.axes.Axes or None
        Axes to plot on.  A new figure is created when *None*.
    bins : int
        Number of histogram bins.

    Returns
    -------
    fig, ax
    """
    import matplotlib.pyplot as plt

    pik_a = result['pik_a']
    pik_b = result['pik_b']
    label_a = result['label_a']
    label_b = result['label_b']

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.figure

    lo = min(pik_a.min(), pik_b.min())
    hi = max(pik_a.max(), pik_b.max())
    bin_edges = np.linspace(lo, hi, bins + 1)

    ax.hist(pik_a, bins=bin_edges, density=True, alpha=0.4,
            color='tab:blue', label=f'{label_a}  (mean={result["mean_a"]:.2f})')
    ax.hist(pik_b, bins=bin_edges, density=True, histtype='step',
            lw=2, color='tab:red',
            label=f'{label_b}  (mean={result["mean_b"]:.2f})')

    ax.axvline(result['mean_a'], color='tab:blue', ls='--', lw=1.2)
    ax.axvline(result['mean_b'], color='tab:red',  ls='--', lw=1.2)

    ax.set_xlabel('total log-likelihood per draw  (pik)')
    ax.set_ylabel('density')
    ax.set_title(
        f'pik distribution comparison\n'
        f'Δmean={result["delta_mean"]:+.3f}  '
        f'({result["frac_diff"]*100:+.3f}%),  '
        f'KS p={result["ks_pvalue"]:.3f}'
    )
    ax.legend(fontsize=9)

    return fig, ax
