"""
Systematic uncertainty helpers for M2 analysis.

Provides LSS parameter loading, configuration, and task-list generation
for lineshape-scaling systematic fits.
"""

import numpy as np
import pandas as pd

from qpym2.cfg.config_builders import build_input_config
from qpym2.cfg.syst_config import (
    DEFAULT_LSS_DATASETS,
    DEFAULT_LSS_SYST_SAMPLES,
)
from qpym2.utils.logging import debug, log

# Re-export with underscore prefix for backwards compatibility
_DEFAULT_LSS_DATASETS = DEFAULT_LSS_DATASETS
_DEFAULT_LSS_SYST_SAMPLES = DEFAULT_LSS_SYST_SAMPLES

# ---------------------------------------------------------------------------
# LSS systematic configuration
# ---------------------------------------------------------------------------

def get_lss_config(cfg, model_key="default"):
    """Build LSS systematic configuration from Hydra config.

    Extracts the systematic staging directory, dataset list, and
    systematic sample definitions.  Falls back to hard-coded defaults
    when the Hydra YAML does not yet carry a ``systematics.lss`` section.

    Args:
        cfg: Hydra configuration object.
        model_key: Named model key (default: 'default').

    Returns:
        tuple: (syst_stg_cfg, datasets, syst_samples_dict)
    """
    input_cfg = build_input_config(cfg, model_key=model_key)
    syst_stg_cfg = input_cfg.copy()

    if hasattr(cfg, 'systematics') and hasattr(cfg.systematics, 'lss'):
        lss_cfg = cfg.systematics.lss
        datasets = list(range(
            lss_cfg.get('ds_start', 3601),
            lss_cfg.get('ds_end', 3616),
        ))
        syst_samples_dict = dict(lss_cfg.get('samples', _DEFAULT_LSS_SYST_SAMPLES))
    else:
        datasets = list(_DEFAULT_LSS_DATASETS)
        syst_samples_dict = dict(_DEFAULT_LSS_SYST_SAMPLES)

    return syst_stg_cfg, datasets, syst_samples_dict


def get_lss_syst_defs(syst_samples_dict, datasets):
    """Generate the full list of (sname, scfg, ds) systematic parameter tuples.

    These are the input parameters for individual LSS systematic fits.
    Each tuple pairs one systematic sample definition with one dataset.

    Args:
        syst_samples_dict: dict mapping sample name to sample config dict.
        datasets: iterable of dataset identifiers.

    Returns:
        list of (sname, scfg, ds) tuples.
    """
    return [(sname, syst_samples_dict[sname], ds)
            for sname in syst_samples_dict
            for ds in datasets]


# ---------------------------------------------------------------------------
# LSS parameter loading — numpy helpers
# ---------------------------------------------------------------------------

def sigma0_arr(df_lspars, n_ds=15, n_ch=988, ds_start=3601) -> np.ndarray:
    """Flat numpy array of sigma0 values, shape (n_ds * n_ch,).

    Numpy equivalent of gintp.sigma0_arrstr — same lookup logic.
    Missing (ch, ds) entries are left as 0.

    .. deprecated::
        Prefer ``qpym2.io.qls_scaling.read_sigma0``, which returns the
        already-reshaped ``(n_ds, n_ch)`` array directly.
    """
    arr = np.zeros(n_ds * n_ch)
    for k in range(n_ds * n_ch):
        try:
            s0 = float(df_lspars.sigma.sel(
                {'ds': ds_start + k // n_ch, 'ch': 1 + k % n_ch}, method=None))
            arr[k] = s0
        except Exception:
            pass
    return arr


def delta_p1_arr(fns_sigma, covs_sigma, spt=1) -> np.ndarray:
    """Numpy array of delta_p1 values per dataset, shape (n_ds,).

    Numpy equivalent of gintp.delta_p1_arrstr, with configurable sample point.
    """
    from qpym2.io.qls_scaling import sample_pars

    def _fnpars(fn):
        return np.array([fn.GetParameter(i) for i in range(fn.GetNpar())])

    sigma_pars0 = [_fnpars(fn) for fn in fns_sigma]
    sigma_pars_spt = [
        sample_pars(fn, cov, skip_last_par=False, nsamples=1, sample_point=spt)
        for fn, cov in zip(fns_sigma, covs_sigma)
    ]
    return np.array([(c1 - c2)[1] for c1, c2 in zip(sigma_pars_spt, sigma_pars0)])


# ---------------------------------------------------------------------------
# Histogram helpers (private)
# ---------------------------------------------------------------------------

def _empty_hist(hm):
    nbins_u = round((hm.umax - hm.umin) / hm.binsize)
    nbins_v = round((hm.vmax - hm.vmin) / hm.binsize)
    return np.zeros((nbins_u, nbins_v), dtype=float)


def _hist_edges(hm):
    nbins_u = round((hm.umax - hm.umin) / hm.binsize)
    nbins_v = round((hm.vmax - hm.vmin) / hm.binsize)
    u_edges = np.linspace(hm.umin, hm.umax, nbins_u + 1)
    v_edges = np.linspace(hm.vmin, hm.vmax, nbins_v + 1)
    return u_edges, v_edges


def _hist_from_rdf(rdf, hm, filters):
    """Histogram a transformed RDF (must already have u, v, e1, e2, …)."""
    if rdf is None:
        return _empty_hist(hm)
    r = rdf
    for f in filters:
        if f:
            r = r.Filter(f)
    from qpym2.utils.hists import get_hist_from_rdf
    return get_hist_from_rdf(r, hm, rtype='numpy')[0]


def _cxx_filter_to_mask(expr, ns):
    """Evaluate a simple C++ filter expression as a numpy boolean mask.

    Handles the arithmetic/comparison subset used by QPyM2 cut strings:
    ``&&``, ``||``, ``abs()``, ``sqrt()``, and column names.
    The *ns* namespace is fully controlled (no builtins exposed).
    """
    py_expr = expr.replace('&&', ') & (').replace('||', ') | (')
    py_expr = '(' + py_expr + ')'
    return eval(py_expr, {"__builtins__": {}}, ns)  # noqa: S307


def _hist_from_np(uv_dict, hm, filters=None):
    """Histogram transformed numpy arrays with energy-window cuts.

    Parameters
    ----------
    filters : list[str] | None
        Extra C++ filter expressions (e.g. ``cuts['common']``,
        ``cuts['ndbd_only']``) to apply as numpy masks.  Energy-window
        cuts (ecuts) are always applied via direct numpy; pass only
        *additional* physics cuts here.
    """
    if uv_dict is None:
        return _empty_hist(hm)
    u, v = uv_dict['u'], uv_dict['v']
    e1 = (u + v) / np.sqrt(2)
    e2 = (u - v) / np.sqrt(2)

    fmask = (
        (u > hm.umin) &
        (u < hm.umax) &
        (v > 0) &
        (e2 > hm.threshold)
    )
    if getattr(hm, 'skip_roi', None):
        roimin, roimax = hm.skip_roi
        fmask &= ~((u > roimin / np.sqrt(2)) & (u < roimax / np.sqrt(2)))

    if filters:
        ns = {'u': u, 'v': v, 'e1': e1, 'e2': e2,
              'esum': u * np.sqrt(2), 'ediff': v * np.sqrt(2),
              'abs': np.abs, 'sqrt': np.sqrt}
        for f in filters:
            if f:
                fmask &= _cxx_filter_to_mask(f, ns)

    u_f, v_f = u[fmask], v[fmask]
    if len(u_f) == 0:
        return _empty_hist(hm)

    u_edges, v_edges = _hist_edges(hm)
    hist, _, _ = np.histogram2d(u_f, v_f, bins=[u_edges, v_edges])
    return hist


# ---------------------------------------------------------------------------
# LSS transforms over evtrees dicts
# ---------------------------------------------------------------------------

def _apply_dq_transform(evtrees, scfg, ds, lss_params, ds_start=3601):
    """Apply dq (Q-residual) transform to all component RDFs.

    Returns a new dict of RDFs with columns u, v, e1, e2, esum, ediff
    defined.  Events from the target dataset get the polynomial shift;
    all others keep u0, v0.
    """
    from qpym2.io.staging_transforms import residual_sampling_uv

    ptn = scfg['ptn']
    ids = ds - ds_start
    q_mean = np.array(lss_params['q_pars_arr'][ids])
    q_sample = np.array(lss_params[f'q_pars_sample_{ptn}'][ids])

    u_expr, v_expr = residual_sampling_uv(q_mean, q_sample, invars=('u0', 'v0'))
    u_def = f'(ds=={ds} ? ({u_expr}) : u0)'
    v_def = f'(ds=={ds} ? ({v_expr}) : v0)'

    out = {}
    for comp, rdf in evtrees.items():
        if rdf is None:
            out[comp] = None
            continue
        r = rdf.Define('u', u_def).Define('v', v_def)
        r = (r.Define('e1',   '(u+v)/sqrt(2)')
              .Define('e2',   '(u-v)/sqrt(2)')
              .Define('esum', 'u*sqrt(2)')
              .Define('ediff', 'v*sqrt(2)'))
        out[comp] = r
    return out


def _apply_sscl_transform(evtrees, scfg, ds, lss_params, ds_start=3601, rng=None):
    """Apply sscl (sigma-scaling) transform to all component RDFs.

    Per-channel sigma lookup cannot be a simple C++ expression, so each
    component's events are extracted via AsNumpy, transformed, and returned
    as ``{comp: {'u': array, 'v': array}}``  (or None).
    """
    from qpym2.io.staging_transforms import sigma_scaling_np

    if rng is None:
        rng = np.random.default_rng()

    ids = ds - ds_start
    spt = scfg.get('spt', 1)
    dp1_key = 'delta_p1_1r' if spt >= 0 else 'delta_p1_1l'
    sigma0_row = lss_params['sigma0_2d'][ids]
    delta_p1_val = lss_params[dp1_key][ids]

    out = {}
    for comp, rdf in evtrees.items():
        if rdf is None:
            out[comp] = None
            continue
        np_data = rdf.AsNumpy(columns=['u0', 'v0', 'ch1', 'ch2', 'ds'])
        u0 = np_data['u0'].astype(np.float64)
        v0 = np_data['v0'].astype(np.float64)
        u, v = u0.copy(), v0.copy()
        mask = (np_data['ds'] == ds)
        if np.any(mask):
            u[mask], v[mask] = sigma_scaling_np(
                u0[mask], v0[mask],
                np_data['ch1'][mask], np_data['ch2'][mask],
                sigma0_row, delta_p1_val, rng=rng,
            )
        out[comp] = {'u': u, 'v': v}
    return out


# ---------------------------------------------------------------------------
# mctable builder from evtrees
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Nominal mctable builder (shared by both LSS and MC syst runners)
# ---------------------------------------------------------------------------

def build_nominal_mctable(mcfit_table, input_cfg, fit_cfg, ndbd_fname='ndbd.root'):
    """Build nominal mctable with mchist filled from staging ROOT files.

    Reads each background component's histogram from its staging ROOT file
    using the energy-window and physics cuts from *fit_cfg*.  Also reads the
    ndbd signal histogram.

    This helper is shared by ``run_syst_mc`` (needs histograms for every
    systematic sample) and can be used anywhere the nominal histograms are
    needed without running a transform.

    Args:
        mcfit_table: DataFrame returned by ``read_jagsh5``; indexed by source
            name with a ``fname`` column giving the staging ROOT basename.
        input_cfg: Input configuration dict (``outpath``, ``treename``).
        fit_cfg: Dict with ``h2uv_model`` and ``cuts`` keys (from ``model_cfg``).
        ndbd_fname: Staging ROOT file basename for the ndbd signal
            (default: ``'ndbd.root'``).

    Returns:
        tuple: ``(mctable, signal_df)`` — mctable copy with ``mchist`` and
        ``integral`` columns populated; signal_df with the same columns for
        the ndbd component.
    """
    import os
    from qpym2.io import staging
    from qpym2.utils.hist_model import get_hm_cuts

    hm = fit_cfg['h2uv_model']
    cuts = fit_cfg['cuts']
    staging_dir = input_cfg['outpath']
    treename = input_cfg.get('treename', 'uvtree')

    ecuts = get_hm_cuts(hm, inpars=('u0', 'v0'))
    defs = [('u', 'u0'), ('v', 'v0')]
    base_filters = [ecuts]
    if cuts.get('common'):
        base_filters.append(cuts['common'])

    filters_mc = [*base_filters, cuts['mc_only']] if cuts.get('mc_only') else base_filters
    filters_ndbd = [*base_filters, cuts['ndbd_only']] if cuts.get('ndbd_only') else base_filters

    mctable = mcfit_table.copy()
    mchists = []
    for comp_name in mctable.index:
        fpath = os.path.join(staging_dir, mctable.at[comp_name, 'fname'])
        h = staging.read_hist(fpath, hm, treename, defs, filters_mc)[0]
        mchists.append(h)
    mctable['mchist'] = mchists
    mctable['integral'] = [float(h.sum()) for h in mchists]

    ndbd_fpath = os.path.join(staging_dir, ndbd_fname)
    ndbd_h = staging.read_hist(ndbd_fpath, hm, treename, defs, filters_ndbd)[0]
    signal_df = pd.DataFrame(
        {'mchist': [ndbd_h], 'integral': [float(ndbd_h.sum())]},
        index=['ndbd'],
    )

    return mctable, signal_df


def make_lss_mctable(evtrees, mcfit_table, fit_cfg, scfg, ds, lss_params,
                     ds_start=3601, rng=None):
    """Apply LSS transform then histogram each component.

    1. Apply the transform once across all component RDFs/arrays.
    2. Histogram transformed data with energy-window + physics cuts.
    """
    from qpym2.utils.hist_model import get_hm_cuts

    hm = fit_cfg['h2uv_model']
    cuts = fit_cfg['cuts']
    stype = scfg['type']

    if stype == 'dq':
        transformed = _apply_dq_transform(evtrees, scfg, ds, lss_params, ds_start)
    elif stype == 'sscl':
        transformed = _apply_sscl_transform(evtrees, scfg, ds, lss_params, ds_start, rng)
    else:
        raise ValueError(f"Unknown syst type: {stype}")

    ecuts = get_hm_cuts(hm, inpars=('u', 'v'))
    base_filters = [ecuts, cuts.get('common', '')]
    filters_mc = [*base_filters, cuts['mc_only']] if cuts.get('mc_only') else base_filters
    filters_ndbd = [*base_filters, cuts['ndbd_only']] if cuts.get('ndbd_only') else base_filters

    def _hist(comp, filters):
        t = transformed.get(comp)
        if stype == 'dq':
            return _hist_from_rdf(t, hm, filters)
        else:
            return _hist_from_np(t, hm)

    mctable = mcfit_table.copy()
    mchists = [_hist(comp, filters_mc) for comp in mctable.index]
    mctable['mchist'] = mchists
    mctable['integral'] = [float(h.sum()) for h in mchists]

    ndbd_h = _hist('ndbd', filters_ndbd)
    signal_df = pd.DataFrame(
        {'mchist': [ndbd_h], 'integral': [float(ndbd_h.sum())]},
        index=['ndbd'],
    )

    return mctable, signal_df


def make_lss_mctable_np(evdata, mcfit_table, fit_cfg, scfg, ds, lss_params,
                        ds_start=3601, rng=None, q_pars_sample=None):
    """Apply LSS transform (numpy-only) then histogram each component.

    Same contract as :func:`make_lss_mctable` but operates on
    pre-materialized numpy arrays (*evdata*) instead of live RDataFrame
    handles.  All component data must have been extracted once via
    ``rdf.AsNumpy`` beforehand.

    Parameters
    ----------
    evdata : dict[str, dict[str, ndarray] | None]
        ``{comp_name: {'u0': arr, 'v0': arr, 'ch1': arr, 'ch2': arr,
        'ds': arr}}``.
    q_pars_sample : list[array] | None
        Optional list of sampled Q-polynomial parameters per dataset.
        When provided, used directly for dq transforms instead of
        looking up ``q_pars_sample_{ptn}`` from *lss_params*.
    """
    from qpym2.io.staging_transforms import (
        apply_sigma_scaling_ds, apply_residual_sampling_ds,
    )

    hm = fit_cfg['h2uv_model']
    cuts = fit_cfg['cuts']
    stype = scfg['type']

    # Extra filters beyond ecuts (which _hist_from_np handles internally).
    extra_mc = [f for f in [cuts.get('common'), cuts.get('mc_only')] if f]
    extra_ndbd = [f for f in [cuts.get('common'), cuts.get('ndbd_only')] if f]

    def _transform_comp(d):
        """Return ``{'u': arr, 'v': arr}`` after LSS transform."""
        if d is None:
            return None
        if stype == 'sscl':
            spt = scfg.get('spt', 1)
            dp1_key = 'delta_p1_1r' if spt >= 0 else 'delta_p1_1l'
            u, v = apply_sigma_scaling_ds(
                d['u0'], d['v0'], d['ch1'], d['ch2'], d['ds'],
                lss_params['sigma0_2d'], lss_params[dp1_key],
                target_ds=ds, ds_start=ds_start, rng=rng,
            )
        elif stype == 'dq':
            if q_pars_sample is not None:
                qps = q_pars_sample
            else:
                if 'ptn' not in scfg:
                    raise ValueError(
                        "dq transform requires either q_pars_sample kwarg "
                        "or scfg['ptn'] key")
                ptn = scfg['ptn']
                qps = lss_params[f'q_pars_sample_{ptn}']
            u, v = apply_residual_sampling_ds(
                d['u0'], d['v0'], d['ds'],
                lss_params['q_pars_arr'],
                qps,
                target_ds=ds, ds_start=ds_start,
            )
        else:
            raise ValueError(f"Unknown syst type: {stype}")
        return {'u': u, 'v': v}

    mctable = mcfit_table.copy()
    mchists = [_hist_from_np(_transform_comp(evdata.get(comp)), hm, filters=extra_mc)
               for comp in mctable.index]
    mctable['mchist'] = mchists
    mctable['integral'] = [float(h.sum()) for h in mchists]

    ndbd_h = _hist_from_np(_transform_comp(evdata.get('ndbd')), hm, filters=extra_ndbd)
    signal_df = pd.DataFrame(
        {'mchist': [ndbd_h], 'integral': [float(ndbd_h.sum())]},
        index=['ndbd'],
    )

    return mctable, signal_df


# ---------------------------------------------------------------------------
# ndbd energy-shift systematic
# ---------------------------------------------------------------------------

def make_ndbd_shift_signal_df(evdata_ndbd, ndbd_shift, hm, filters=None):
    """Build shifted signal_df by applying an energy shift to ndbd events.

    Uses the same shift formula as ``get_ushiftp_uvtr`` in
    ``qpym2.fit.create_m2uv`` but operates on numpy arrays:

        de = (u0 + v0) / 2 * ndbd_shift / 1000 / sqrt(2)
        u  = u0 + de
        v  = v0 + de

    Parameters
    ----------
    evdata_ndbd : dict or None
        ``{'u0': arr, 'v0': arr, ...}`` — raw ndbd event arrays.
    ndbd_shift : float
        Energy shift in keV (same scale as ``co_shift``).
    hm : HistType
        Histogram model (binning, ROI, threshold).
    filters : list[str] or None
        Extra filter expressions (e.g. ``cuts['ndbd_only']``).

    Returns
    -------
    signal_df : DataFrame
        Single-row DataFrame indexed ``['ndbd']`` with ``mchist`` and
        ``integral`` columns.
    """
    if evdata_ndbd is None:
        h = _empty_hist(hm)
    else:
        u0 = evdata_ndbd['u0'].astype(np.float64)
        v0 = evdata_ndbd['v0'].astype(np.float64)
        if ndbd_shift != 0.0:
            de = (u0 + v0) / 2.0 * (ndbd_shift / 1000.0) / np.sqrt(2)
            u = u0 + de
            v = v0 + de
        else:
            u, v = u0, v0
        h = _hist_from_np({'u': u, 'v': v}, hm, filters=filters)

    return pd.DataFrame(
        {'mchist': [h], 'integral': [float(h.sum())]},
        index=['ndbd'],
    )


# ---------------------------------------------------------------------------
# Generic extra-cut mctable & data builders
# ---------------------------------------------------------------------------

def build_ecut_mctable(mcfit_table, input_cfg, fit_cfg,
                       extra_cuts=None, ndbd_fname='ndbd.root'):
    """Build mctable from staging ROOT files with additional filter cuts.

    Same as :func:`build_nominal_mctable` but appends *extra_cuts* to the
    filter chain applied when reading each component's histogram.  Useful
    for studying the effect of removing specific energy regions (e.g.
    ``abs(e2 - 511) > 5`` to veto annihilation photons).

    Parameters
    ----------
    mcfit_table : DataFrame
        From ``read_jagsh5``; indexed by source name with ``fname`` column.
    input_cfg : dict
        ``{'outpath': staging_dir, 'treename': 'uvtree'}``.
    fit_cfg : dict
        ``{'h2uv_model': HistType, 'cuts': {...}}`` from ``model_cfg``.
    extra_cuts : list[str] or None
        Additional C++ filter expressions to apply (e.g.
        ``["abs(e2 - 511) > 5"]``).
    ndbd_fname : str
        Staging ROOT basename for the ndbd signal.

    Returns
    -------
    tuple : (mctable, signal_df)
    """
    import os
    from qpym2.io import staging
    from qpym2.utils.hist_model import get_hm_cuts

    hm = fit_cfg['h2uv_model']
    cuts = fit_cfg['cuts']
    staging_dir = input_cfg['outpath']
    treename = input_cfg.get('treename', 'uvtree')

    ecuts = get_hm_cuts(hm, inpars=('u0', 'v0'))
    defs = [('u', 'u0'), ('v', 'v0')]
    base_filters = [ecuts]
    if cuts.get('common'):
        base_filters.append(cuts['common'])
    if extra_cuts:
        base_filters.extend(extra_cuts)

    filters_mc = ([*base_filters, cuts['mc_only']]
                  if cuts.get('mc_only') else base_filters)
    filters_ndbd = ([*base_filters, cuts['ndbd_only']]
                    if cuts.get('ndbd_only') else base_filters)

    mctable = mcfit_table.copy()
    mchists = []
    for comp_name in mctable.index:
        fpath = os.path.join(staging_dir, mctable.at[comp_name, 'fname'])
        h = staging.read_hist(fpath, hm, treename, defs, filters_mc)[0]
        mchists.append(h)
    mctable['mchist'] = mchists
    mctable['integral'] = [float(h.sum()) for h in mchists]

    ndbd_fpath = os.path.join(staging_dir, ndbd_fname)
    ndbd_h = staging.read_hist(ndbd_fpath, hm, treename, defs, filters_ndbd)[0]
    signal_df = pd.DataFrame(
        {'mchist': [ndbd_h], 'integral': [float(ndbd_h.sum())]},
        index=['ndbd'],
    )

    return mctable, signal_df


def build_ecut_data(input_cfg, fit_cfg, extra_cuts=None):
    """Read data event list from staging ROOT with additional filter cuts.

    Same as ``create_evlist`` in :mod:`qpym2.fit.create_m2uv` but accepts
    *extra_cuts* to append to the filter chain.

    Parameters
    ----------
    input_cfg : dict
        ``{'outpath': staging_dir, 'data_fname': ..., 'treename': ...}``.
    fit_cfg : dict
        ``{'h2uv_model': HistType, 'cuts': {...}}`` from ``model_cfg``.
    extra_cuts : list[str] or None
        Additional C++ filter expressions.

    Returns
    -------
    ndarray : shape (N, 2) array of observed ``[u, v]`` coordinates.
    """
    import os
    from qpym2.io import staging
    from qpym2.utils.hist_model import get_hm_cuts

    hm = fit_cfg['h2uv_model']
    cuts = fit_cfg['cuts']
    staging_dir = input_cfg['outpath']
    treename = input_cfg.get('treename', 'uvtree')
    data_fname = input_cfg['data_fname']

    data_fpath = os.path.join(staging_dir, data_fname)
    ecuts = get_hm_cuts(hm, inpars=('u0', 'v0'))
    defs = [('u', 'u0'), ('v', 'v0')]

    filters = [ecuts]
    if cuts.get('common'):
        filters.append(cuts['common'])
    if cuts.get('data_only'):
        filters.append(cuts['data_only'])
    if extra_cuts:
        filters.extend(extra_cuts)

    return staging.read_evlist(data_fpath, treename,
                               defs=defs, filters=filters,
                               rtype='numpy')
