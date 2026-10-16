"""syst_lss.py — Modular LSS systematic fit functions.

Layered API for running lineshape systematic (LSS) fits:

    **Transform** (per-component):
        transform_sscl, transform_dq

    **Histogram** (all components → mctable):
        build_sscl_mctable, build_dq_mctable

    **Fit** (shared):
        run_lss_fit

    **Post-processing** (shared):
        extract_posterior_stats, apply_bias_correction

    **Full-pipeline convenience**:
        run_sscl_fit, run_dq_fit

All public functions take explicit, pre-indexed physics parameters.
No ``ds_start``, ``scfg``, or ``lss_params`` blobs.
"""

import numpy as np
import pandas as pd

from qpym2.io.staging_transforms import sigma_scaling_np, residual_sampling_np
from qpym2.fit.syst import _hist_from_np, _empty_hist


# ---------------------------------------------------------------------------
# Transform layer (per-component)
# ---------------------------------------------------------------------------

def transform_sscl(mc_events, target_ds, ch_sigma0, delta_p1, rng=None):
    """Apply sigma-scaling transform to a single component's events.

    Parameters
    ----------
    mc_events : dict or None
        Single-component event arrays:
        ``{'u0': arr, 'v0': arr, 'ch1': arr, 'ch2': arr, 'ds': arr}``.
    target_ds : int
        Dataset number to transform (events with other ds are kept as-is).
    ch_sigma0 : ndarray, shape (n_ch,)
        Per-channel relative sigma for this dataset.
    delta_p1 : float
        Delta-p1 value for this dataset.
    rng : np.random.Generator or None
        Random generator for stochastic smearing.

    Returns
    -------
    dict or None
        ``{'u': arr, 'v': arr}`` with transformed coordinates,
        or None if *mc_events* is None.
    """
    if mc_events is None:
        return None
    u0, v0 = mc_events['u0'], mc_events['v0']
    ch1, ch2 = mc_events['ch1'], mc_events['ch2']
    ds_arr = mc_events['ds']

    u_out = np.array(u0, dtype=float)
    v_out = np.array(v0, dtype=float)
    mask = np.asarray(ds_arr) == target_ds
    if np.any(mask):
        u_out[mask], v_out[mask] = sigma_scaling_np(
            u0[mask], v0[mask], ch1[mask], ch2[mask],
            ch_sigma0, delta_p1, rng=rng,
        )
    return {'u': u_out, 'v': v_out}


def transform_dq(mc_events, target_ds, q_pars_mean, q_pars_sample):
    """Apply Q-residual polynomial transform to a single component's events.

    Parameters
    ----------
    mc_events : dict or None
        Single-component event arrays:
        ``{'u0': arr, 'v0': arr, 'ds': arr}``.
    target_ds : int
        Dataset number to transform (events with other ds are kept as-is).
    q_pars_mean : array-like, shape (>=3,)
        Mean Q-polynomial parameters for this dataset.
    q_pars_sample : array-like, shape (>=3,)
        Sampled Q-polynomial parameters for this dataset.

    Returns
    -------
    dict or None
        ``{'u': arr, 'v': arr}`` with transformed coordinates,
        or None if *mc_events* is None.
    """
    if mc_events is None:
        return None
    u0, v0 = mc_events['u0'], mc_events['v0']
    ds_arr = mc_events['ds']

    u_out = np.array(u0, dtype=float)
    v_out = np.array(v0, dtype=float)
    mask = np.asarray(ds_arr) == target_ds
    if np.any(mask):
        u_out[mask], v_out[mask] = residual_sampling_np(
            u0[mask], v0[mask], q_pars_mean, q_pars_sample,
        )
    return {'u': u_out, 'v': v_out}


# ---------------------------------------------------------------------------
# Histogramming layer (all components → mctable)
# ---------------------------------------------------------------------------

def _build_mctable(mc_events_dict, mcfit_table, hm, transform_comp,
                   cuts_mc=None, cuts_ndbd=None):
    """Shared histogramming logic for build_sscl_mctable / build_dq_mctable.

    Parameters
    ----------
    transform_comp : callable
        ``(mc_events) -> {'u': arr, 'v': arr} | None``
    cuts_mc : list[str] or None
        Filter expressions for background components.
    cuts_ndbd : list[str] or None
        Filter expressions for the ndbd component.
    """
    mctable = mcfit_table.copy()
    mchists = [
        _hist_from_np(
            transform_comp(mc_events_dict.get(comp)),
            hm, filters=cuts_mc,
        )
        for comp in mctable.index
    ]
    mctable['mchist'] = mchists
    mctable['integral'] = [float(h.sum()) for h in mchists]

    ndbd_h = _hist_from_np(
        transform_comp(mc_events_dict.get('ndbd')),
        hm, filters=cuts_ndbd,
    )
    signal_df = pd.DataFrame(
        {'mchist': [ndbd_h], 'integral': [float(ndbd_h.sum())]},
        index=['ndbd'],
    )
    return mctable, signal_df


def build_sscl_mctable(mc_events_dict, mcfit_table, hm, target_ds,
                        ch_sigma0, delta_p1, cuts_mc=None, cuts_ndbd=None,
                        rng=None):
    """Transform all components with sigma-scaling and histogram.

    Parameters
    ----------
    mc_events_dict : dict[str, dict | None]
        ``{comp_name: {'u0','v0','ch1','ch2','ds'} | None}``.
    mcfit_table : DataFrame
        Background MC fit table (indexed by component name).
    hm : HistType
        Histogram model (binning, ROI, threshold).
    target_ds : int
        Dataset number to transform.
    ch_sigma0 : ndarray, shape (n_ch,)
        Per-channel relative sigma for this dataset.
    delta_p1 : float
        Delta-p1 value for this dataset.
    cuts_mc : list[str] or None
        Filter expressions for background components.
    cuts_ndbd : list[str] or None
        Filter expressions for ndbd component.
    rng : np.random.Generator or None

    Returns
    -------
    (mctable, signal_df)
    """
    def _transform(mc_events):
        return transform_sscl(mc_events, target_ds, ch_sigma0, delta_p1,
                              rng=rng)
    return _build_mctable(mc_events_dict, mcfit_table, hm, _transform,
                          cuts_mc=cuts_mc, cuts_ndbd=cuts_ndbd)


def build_dq_mctable(mc_events_dict, mcfit_table, hm, target_ds,
                     q_pars_mean, q_pars_sample,
                     cuts_mc=None, cuts_ndbd=None):
    """Transform all components with Q-residual and histogram.

    Parameters
    ----------
    mc_events_dict : dict[str, dict | None]
        ``{comp_name: {'u0','v0','ds'} | None}``.
    mcfit_table : DataFrame
        Background MC fit table (indexed by component name).
    hm : HistType
        Histogram model (binning, ROI, threshold).
    target_ds : int
        Dataset number to transform.
    q_pars_mean : array-like, shape (>=3,)
        Mean Q-polynomial parameters for this dataset.
    q_pars_sample : array-like, shape (>=3,)
        Sampled Q-polynomial parameters for this dataset.
    cuts_mc : list[str] or None
        Filter expressions for background components.
    cuts_ndbd : list[str] or None
        Filter expressions for ndbd component.

    Returns
    -------
    (mctable, signal_df)
    """
    def _transform(mc_events):
        return transform_dq(mc_events, target_ds, q_pars_mean, q_pars_sample)
    return _build_mctable(mc_events_dict, mcfit_table, hm, _transform,
                          cuts_mc=cuts_mc, cuts_ndbd=cuts_ndbd)


# ---------------------------------------------------------------------------
# Fit layer
# ---------------------------------------------------------------------------

def run_lss_fit(mctable, signal_df, data, hm, comp_cfg, nominal_priors,
                chains, draws):
    """Run a single MCMC fit on pre-built LSS mctable.

    Uses priors read from the nominal fit rather than recomputing them.

    Parameters
    ----------
    mctable : DataFrame
        Background MC table with ``mchist`` column (from build_*_mctable).
    signal_df : DataFrame
        Signal (ndbd) table with ``mchist`` column.
    data : ndarray, shape (n_obs, 2)
        Observed event coordinates ``[[u, v], ...]``.
    hm : HistType
        Histogram model.
    comp_cfg : list[dict]
        Component configuration (filter_regex for grouping).
    nominal_priors : list[tuple]
        Prior tuples from the nominal fit's comps table, e.g.
        ``[('halfnormal', {'mu': 150, 'sigma': 15}), ..., ('flat', {'lower': 0, 'upper': N})]``.
        Length must equal ``len(comp_cfg) + 1`` (background components + ndbd).
    chains : int
    draws : int

    Returns
    -------
    trace : InferenceData
        ArviZ trace object.
    """
    from qpym2.models.m2 import make_bkgcomps_table, make_fik, m2_model
    from qpym2.utils import hists as hist_tools

    # Build background components table (groups MC sources, applies smoothing)
    fit_comps = make_bkgcomps_table(comp_cfg, mctable, smooth=hm.smooth)

    # Prepare signal row
    sig_hist = signal_df['mchist'].iloc[0]
    if hm.smooth:
        sig_hist = hist_tools.smooth_nph2(sig_hist, smooth=hm.smooth)
    signal_df = signal_df.copy()
    signal_df['pdf_hist'] = [sig_hist / signal_df['integral'].iloc[0]]

    # Concat and set nominal priors
    comps = pd.concat([fit_comps, signal_df], join='inner')
    comps['prior'] = nominal_priors
    comps['varname'] = comps.index

    # Build likelihood matrix and model
    _, xe, ye = hist_tools.get_empty_hist(hm, return_numpy=True)
    fik = make_fik(comps, data, xe, ye)
    model = m2_model(comps, fik)

    # Run MCMC
    import pymc as pm
    with model:
        trace = pm.sample(
            draws=draws, chains=chains, tune=1000,
            return_inferencedata=True, target_accept=0.9,
            idata_kwargs={"log_likelihood": True},
            progressbar=False, cores=1,
        )
    return trace


# ---------------------------------------------------------------------------
# Post-processing layer
# ---------------------------------------------------------------------------

def extract_posterior_stats(trace):
    """Extract summary statistics from an MCMC trace.

    Returns
    -------
    dict with keys:
        component_modes : dict[str, float]
        ndbd_samples : ndarray
        p16 : float
        p90 : float
    """
    from qpym2.utils.mathy import get_mode

    component_modes = {}
    for var_name in trace.posterior.data_vars:
        samples = trace.posterior[var_name].values.flatten()
        component_modes[var_name] = float(get_mode(samples))

    ndbd_samples = trace.posterior['ndbd'].values.flatten()
    p16 = float(np.percentile(ndbd_samples, 16))
    p90 = float(np.percentile(ndbd_samples, 90))

    return {
        'component_modes': component_modes,
        'ndbd_samples': ndbd_samples,
        'p16': p16,
        'p90': p90,
    }


def apply_bias_correction(ndbd_samples, bias_popt, bias_pcov):
    """Apply bias correction to ndbd posterior samples.

    Returns
    -------
    dict with keys: mode_bc, p16_bc, p90_bc
    """
    from qpym2.fit.bias import apply_bias_correction
    from qpym2.utils.mathy import get_mode

    corrected = apply_bias_correction(ndbd_samples, bias_popt)
    return {
        'mode_bc': float(get_mode(corrected)),
        'p16_bc': float(np.percentile(corrected, 16)),
        'p90_bc': float(np.percentile(corrected, 90)),
    }


# ---------------------------------------------------------------------------
# Full-pipeline convenience functions
# ---------------------------------------------------------------------------

def run_sscl_fit(mc_events_dict, mcfit_table, data, hm, comp_cfg,
                 nominal_priors, target_ds, ch_sigma0, delta_p1, chains, draws,
                 cuts_mc=None, cuts_ndbd=None, rng=None):
    """Full sscl pipeline: transform → histogram → fit.

    Returns
    -------
    trace : InferenceData
    """
    mctable, signal_df = build_sscl_mctable(
        mc_events_dict, mcfit_table, hm, target_ds,
        ch_sigma0, delta_p1,
        cuts_mc=cuts_mc, cuts_ndbd=cuts_ndbd, rng=rng,
    )
    return run_lss_fit(mctable, signal_df, data, hm, comp_cfg,
                       nominal_priors, chains, draws)


def run_dq_fit(mc_events_dict, mcfit_table, data, hm, comp_cfg,
               nominal_priors, target_ds, q_pars_mean, q_pars_sample,
               chains, draws, cuts_mc=None, cuts_ndbd=None):
    """Full dq pipeline: transform → histogram → fit.

    Returns
    -------
    trace : InferenceData
    """
    mctable, signal_df = build_dq_mctable(
        mc_events_dict, mcfit_table, hm, target_ds,
        q_pars_mean, q_pars_sample,
        cuts_mc=cuts_mc, cuts_ndbd=cuts_ndbd,
    )
    return run_lss_fit(mctable, signal_df, data, hm, comp_cfg,
                       nominal_priors, chains, draws)
