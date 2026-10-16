from typing import List, Optional, Tuple

import numpy as np

from qpym2.utils.logging import debug

_s2 = np.sqrt(2)


def residual_sampling_uv(rspars_mean, rspars_sample, invars=('u0', 'v0')) -> Tuple[str, str]:
    """ Get the uv transform for the shift type and value.

    Args:
        invars (tuple): input (u, v) to transform

    Returns:
        tuple: (u, v) transform
    """
    u, v = invars
    e1,e2 = f'({u}+{v})/sqrt(2)', f'({u}-{v})/sqrt(2)'  
    a, b, c = rspars_sample - rspars_mean

    tru = f'{u} + {a}*sqrt(2) + {b}*({u}) + {c}/sqrt(2)*(({u})*({u})+({v})*({v}))'
    trv = f'{v}               + {b}*({v}) + {c}*sqrt(2)*({u})*({v})'

    return (tru, trv)


def sigma_scaling_np(u0, v0, ch1, ch2, ch_sigma0, delta_p1, rng=None):
    """Sigma scaling for a single dataset.

    Args:
        u0, v0   : (n_events,) coordinates
        ch1, ch2 : (n_events,) 1-indexed channel numbers
        sigma0   : (n_ch,) relative sigma array; index i = channel (i+1)
        delta_p1 : scalar delta for this dataset
        rng      : np.random.Generator
    """
    if rng is None:
        rng = np.random.default_rng()
    ich1 = np.asarray(ch1) - 1
    ich2 = np.asarray(ch2) - 1
    e1 = (u0 + v0) / _s2
    e2 = (u0 - v0) / _s2
    dse1 = e1 * ch_sigma0[ich1] * delta_p1
    dse2 = e2 * ch_sigma0[ich2] * delta_p1
    e1a = rng.normal(e1, np.abs(dse1))
    e2a = rng.normal(e2, np.abs(dse2))
    return (e1a + e2a) / _s2, (e1a - e2a) / _s2


def sigma_scaling_all_np(u0, v0, ch1, ch2, ds_arr, ch_sigma0, delta_p1, ds_start=3601, rng=None):
    """Vectorised sigma scaling over events from multiple datasets.

    Args:
        u0, v0   : (n_events,) coordinates
        ch1, ch2 : (n_events,) physical channel numbers [1-988, < 988]
        ds_arr   : (n_events,) dataset numbers
        ch_sigma0 : (n_ds, n_ch) reindexed so [ids, ich] is valid
        delta_p1 : (n_ds,) per-dataset delta values
        ds_start : first dataset number (default 3601)
        rng      : np.random.Generator
    """
    if rng is None:
        rng = np.random.default_rng()
    ids  = np.asarray(ds_arr) - ds_start
    ich1 = np.asarray(ch1) - 1
    ich2 = np.asarray(ch2) - 1
    e1 = (u0 + v0) / _s2
    e2 = (u0 - v0) / _s2
    dse1 = e1 * ch_sigma0[ids, ich1] * delta_p1[ids]
    dse2 = e2 * ch_sigma0[ids, ich2] * delta_p1[ids]
    e1a = rng.normal(e1, np.abs(dse1))
    e2a = rng.normal(e2, np.abs(dse2))
    return (e1a + e2a) / _s2, (e1a - e2a) / _s2


def residual_sampling_np(u0, v0, rspars_mean, rspars_sample):
    """Q-residual polynomial shift, mirrors residual_sampling_uv.

        u' = u + a*√2  + b*u + (c/√2)*(u²+v²)
        v' = v         + b*v + c*√2*u*v

    Args:
        u0, v0         : (n_events,) coordinates
        rspars_mean    : (≥3,) mean parameters [a, b, c, ...]
        rspars_sample  : (≥3,) sampled parameters
    """
    a, b, c = (rspars_sample - rspars_mean)[:3]
    u = u0 + a * _s2 + b * u0 + (c / _s2) * (u0**2 + v0**2)
    v = v0           + b * v0 +  c * _s2   * u0 * v0
    return u, v


def residual_sampling_all_np(u0, v0, ds_arr, qpars_arr, qpars_1r_arr, ds_start=3601):
    """Vectorised Q-residual shift over events from multiple datasets.

    Args:
        u0, v0        : (n_events,) coordinates
        ds_arr        : (n_events,) dataset numbers
        qpars_arr     : (n_ds, ≥3) mean parameters per dataset
        qpars_1r_arr  : (n_ds, ≥3) sampled parameters per dataset
        ds_start      : first dataset number (default 3601)
    """
    ids = np.asarray(ds_arr) - ds_start
    abc = (np.asarray(qpars_1r_arr) - np.asarray(qpars_arr))[ids, :3]  # (n_events, 3)
    a, b, c = abc[:, 0], abc[:, 1], abc[:, 2]
    u = u0 + a * _s2 + b * u0 + (c / _s2) * (u0**2 + v0**2)
    v = v0           + b * v0 +  c * _s2   * u0 * v0
    return u, v


# ---------------------------------------------------------------------------
# RDF define-list factories (used by apply_lss_syst_defs in create_evtrees)
# ---------------------------------------------------------------------------

def sigma_scaling_defs(var_suffix: str, fn_dsigma: str) -> List[Tuple[str, str]]:
    """RDF define list for the sigma-scaling (sscl) systematic transform.

    Mirrors ``sigma_scaling_np`` but expressed as a sequence of
    ``(alias, ROOT-C++ expr)`` pairs suitable for ``RDataFrame.Define``.
    The JIT functions ``{fn_dsigma}`` and ``gintp_fn_sample_gaus`` must
    already be declared via ``register_lss_jit`` before these defines are
    evaluated by RDF.

    The define sequence:
        e1 = (u0+v0)/√2    (energy of higher-energy detector)
        e2 = (u0-v0)/√2    (energy of lower-energy detector)
        de1 = fn_dsigma(ch1, ds, e1)        (smearing width for e1)
        de2 = fn_dsigma(ch2, ds, e2)        (smearing width for e2)
        e1s = Gaus(e1, |de1|)               (smeared e1)
        e2s = Gaus(e2, |de2|)               (smeared e2)
        u_<var_suffix> = (e1s+e2s)/√2
        v_<var_suffix> = (e1s-e2s)/√2

    Args:
        var_suffix: Suffix appended to all intermediate and output column
            names (e.g. ``'sscl'`` → ``u_sscl``, ``v_sscl``).
            Use a unique suffix per systematic point to avoid clashes.
        fn_dsigma: Name of the JIT function that computes the smearing
            width (e.g. ``'gintp_fn_ls_dsigma_1r'``).

    Returns:
        List of ``(alias, expr)`` pairs to pass to :func:`apply_defines`.
    """
    s = var_suffix
    return [
        (f'e1_{s}',  '(u0 + v0) / sqrt(2)'),
        (f'e2_{s}',  '(u0 - v0) / sqrt(2)'),
        (f'de1_{s}', f'{fn_dsigma}(ch1, ds, e1_{s})'),
        (f'de2_{s}', f'{fn_dsigma}(ch2, ds, e2_{s})'),
        (f'e1s_{s}', f'gintp_fn_sample_gaus(e1_{s}, abs(de1_{s}))'),
        (f'e2s_{s}', f'gintp_fn_sample_gaus(e2_{s}, abs(de2_{s}))'),
        (f'u_{s}',   f'(e1s_{s} + e2s_{s}) / sqrt(2)'),
        (f'v_{s}',   f'(e1s_{s} - e2s_{s}) / sqrt(2)'),
    ]


def residual_sampling_defs(
    rspars_mean,
    rspars_sample,
    var_suffix: str,
    invars: Tuple[str, str] = ('u0', 'v0'),
) -> List[Tuple[str, str]]:
    """RDF define list for the Q-residual (dq) systematic transform.

    Thin wrapper around :func:`residual_sampling_uv`: converts the two
    returned C++ expression strings into a ``[(alias, expr)]`` list
    compatible with :func:`apply_defines`.

    Args:
        rspars_mean:   Mean Q-polynomial parameters (≥3 elements).
        rspars_sample: Sampled Q-polynomial parameters (≥3 elements).
        var_suffix:    Suffix for the output columns (e.g. ``'dq1r_ds3601'``).
        invars:        Input ``(u, v)`` column names (default ``('u0','v0')``).

    Returns:
        ``[(f'u_{var_suffix}', tru_expr), (f'v_{var_suffix}', trv_expr)]``
    """
    tru, trv = residual_sampling_uv(rspars_mean, rspars_sample, invars)
    return [(f'u_{var_suffix}', tru), (f'v_{var_suffix}', trv)]


# ---------------------------------------------------------------------------
# Per-dataset numpy transforms (for fork-safe Pool dispatch)
# ---------------------------------------------------------------------------

def apply_sigma_scaling_ds(
    u0, v0, ch1, ch2, ds_arr,
    sigma0_2d,
    delta_p1,
    target_ds: int,
    ds_start: int = 3601,
    rng=None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply sigma-scaling transform to events from a single target dataset.

    Events whose ``ds_arr != target_ds`` are left unchanged.  This is the
    numpy equivalent of the RDF ternary ``(ds==target_ds ? u_sscl : u0)``
    plus the ``sigma_scaling_defs`` column chain.

    Args:
        u0, v0     : (n_events,) raw coordinates.
        ch1, ch2   : (n_events,) 1-indexed channel numbers.
        ds_arr     : (n_events,) dataset numbers.
        sigma0_2d  : (n_ds, n_ch) relative sigma array.
        delta_p1   : (n_ds,) per-dataset delta values.
        target_ds  : Dataset number to transform.
        ds_start   : First dataset number offset (default 3601).
        rng        : ``np.random.Generator`` (created if None).

    Returns:
        ``(u_out, v_out)`` with transformed values for ``target_ds`` events.
    """
    u0 = np.asarray(u0, dtype=float)
    v0 = np.asarray(v0, dtype=float)
    ds_arr = np.asarray(ds_arr)
    mask = ds_arr == target_ds
    if not np.any(mask):
        return u0.copy(), v0.copy()
    ids = target_ds - ds_start
    u_out, v_out = u0.copy(), v0.copy()
    u_out[mask], v_out[mask] = sigma_scaling_np(
        u0[mask], v0[mask], ch1[mask], ch2[mask],
        sigma0_2d[ids], delta_p1[ids],
        rng=rng,
    )
    return u_out, v_out


def apply_residual_sampling_ds(
    u0, v0, ds_arr,
    qpars_arr,
    qpars_sample_arr,
    target_ds: int,
    ds_start: int = 3601,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply Q-residual transform to events from a single target dataset.

    Events whose ``ds_arr != target_ds`` are left unchanged.  This is the
    numpy equivalent of the RDF ternary ``(ds==target_ds ? u_dq* : u0)``
    plus ``residual_sampling_defs``.

    Args:
        u0, v0           : (n_events,) raw coordinates.
        ds_arr           : (n_events,) dataset numbers.
        qpars_arr        : (n_ds, ≥3) mean Q-polynomial parameters.
        qpars_sample_arr : (n_ds, ≥3) sampled Q-polynomial parameters.
        target_ds        : Dataset number to transform.
        ds_start         : First dataset number offset (default 3601).

    Returns:
        ``(u_out, v_out)`` with transformed values for ``target_ds`` events.
    """
    u0 = np.asarray(u0, dtype=float)
    v0 = np.asarray(v0, dtype=float)
    ds_arr = np.asarray(ds_arr)
    mask = ds_arr == target_ds
    if not np.any(mask):
        return u0.copy(), v0.copy()
    ids = target_ds - ds_start
    qpars = np.asarray(qpars_arr)
    qpars_s = np.asarray(qpars_sample_arr)
    u_out, v_out = u0.copy(), v0.copy()
    u_out[mask], v_out[mask] = residual_sampling_np(
        u0[mask], v0[mask], qpars[ids], qpars_s[ids],
    )
    return u_out, v_out


# ---------------------------------------------------------------------------
# RDF define-string helper (kept for backward compat / RDF-based callers)
# ---------------------------------------------------------------------------

def get_syst_defs(syst_cfg, ds, type, invars=('u0', 'v0')):    
    """ Setup the axis definitions for the given systematic type and dataset.
    
    NOTE: currently only needed for lineshape scaling systematics
    TODO: read params from config
    """
    scfg = syst_cfg
    u_, v_ = invars

    if type is None:
        ud, vd = u_, v_
    elif type == 'dq':
        ptn = scfg['ptn']
        if ptn == '1l':
            ud, vd = f'u_dq1l_ds{ds}', f'v_dq1l_ds{ds}'
        if ptn == '1r':
            ud, vd = f'u_dq1r_ds{ds}', f'v_dq1r_ds{ds}'
    elif type == 'sscl':
        type = scfg['type']
        ud, vd = f'u_sscl', f'v_sscl'
    else:
        raise ValueError(f"Unknown syst type: {type}")
    
    u_ = f'(ds=={ds}? {ud}:{u_})'
    v_ = f'(ds=={ds}? {vd}:{v_})'

    return u_, v_
