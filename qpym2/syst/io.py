"""HDF5 I/O for systematic results and combined marginalised outputs.

Schema
------

Per-systematic file (``save_syst_result`` / ``load_syst_result``, format v2):

::

    <path>.h5
      /raw_post                  dataset (1-D float, optional) — combined raw posterior
      /bc_post                   dataset (1-D float)           — combined bias-corr. posterior
      /per_post                  group   (pandas, optional)    — per-component stats
      attrs (root):
        name           : str
        model_pars     : JSON-encoded dict
        stats          : JSON-encoded dict {'raw': {...}, 'bc': {...}}
                         each inner dict has keys: mode, hdi_lo, hdi_hi, p16, p84, p90
        format_version : '2'

Per-systematic file (legacy ``save_syst_delta`` / ``load_syst_delta``, format v1):
    Kept only for backward-compatible reading of old files.

Combined-result file (``save_combined_result`` / ``load_combined_result``):

::

    <path>.h5
      /nominal_bc_post           dataset (1-D float)  — nominal bias-corr. posterior
      /marg_post                 dataset (1-D float)  — smeared (unclipped)
      /rate_post                 dataset (1-D float)  — rate posterior
      /delta_table               group   (pandas)     — per-syst contributions
      /stats                     group   (pandas)     — count/rate summary
      attrs (root):
        sigma_right    : float
        sigma_left     : float
        model_pars     : JSON-encoded dict
        eff_mc         : float
        exposure       : float
        psa_mean       : float (or NaN)
        psa_sig        : float (or NaN)
        format_version : '1'
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from qpym2.syst.core import SystDelta


_FORMAT_VERSION = '1'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def cs_tag(co_shift: float) -> str:
    """Convert a co_shift float to the file-name tag used everywhere.

    Examples
    --------
    >>> cs_tag(1.0)
    'cs1p0'
    >>> cs_tag(0.3)
    'cs0p3'
    """
    return 'cs' + f'{float(co_shift):.1f}'.replace('.', 'p')


def _to_jsonable(obj):
    """Best-effort conversion for json.dumps (numpy scalars, tuples)."""
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


# ---------------------------------------------------------------------------
# Per-systematic delta I/O
# ---------------------------------------------------------------------------

def save_syst_delta(delta: SystDelta, path: str | Path, overwrite: bool = True) -> Path:
    """Write a :class:`SystDelta` to HDF5.

    Parameters
    ----------
    delta : SystDelta
    path : str or Path
        Destination ``.h5`` file.  Parent directory is created.
    overwrite : bool
        If False and the file exists, :class:`FileExistsError` is raised.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists() and not overwrite:
        raise FileExistsError(p)

    with h5py.File(p, 'w') as f:
        if delta.syst_post is not None:
            f.create_dataset('syst_post',
                             data=np.asarray(delta.syst_post, dtype=float),
                             compression='gzip', compression_opts=4)
        f.attrs['name'] = delta.name
        f.attrs['delta_right'] = float(delta.delta_right)
        f.attrs['delta_left'] = float(delta.delta_left)
        f.attrs['model_pars'] = json.dumps(_to_jsonable(delta.model_pars))
        f.attrs['stats'] = json.dumps(_to_jsonable(delta.stats))
        f.attrs['format_version'] = _FORMAT_VERSION
    return p


def load_syst_delta(path: str | Path) -> SystDelta:
    """Read a :class:`SystDelta` saved by :func:`save_syst_delta`."""
    p = Path(path)
    with h5py.File(p, 'r') as f:
        syst_post = f['syst_post'][:] if 'syst_post' in f else None
        return SystDelta(
            name=str(f.attrs['name']),
            delta_right=float(f.attrs['delta_right']),
            delta_left=float(f.attrs['delta_left']),
            model_pars=json.loads(f.attrs.get('model_pars', '{}')),
            syst_post=syst_post,
            stats=json.loads(f.attrs.get('stats', '{}')),
        )


# ---------------------------------------------------------------------------
# Combined-result I/O
# ---------------------------------------------------------------------------

def save_combined_result(
    path: str | Path,
    nominal_bc_post: np.ndarray,
    marg_post: np.ndarray,
    rate_post: np.ndarray,
    delta_table: pd.DataFrame,
    stats: pd.DataFrame,
    *,
    sigma_right: float,
    sigma_left: float,
    model_pars: dict,
    eff_mc: float,
    exposure: float,
    psa_mean: float | None = None,
    psa_sig: float | None = None,
    overwrite: bool = True,
) -> Path:
    """Save a full marginalised+rate result.

    The two ``DataFrame`` arguments are written via ``pandas.HDFStore`` so
    they can be re-read as DataFrames (round-trip-friendly).  The three array
    payloads are written via ``h5py`` for compactness.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        if not overwrite:
            raise FileExistsError(p)
        p.unlink()

    # 1) arrays + attrs via h5py
    with h5py.File(p, 'w') as f:
        for k, arr in (
            ('nominal_bc_post', nominal_bc_post),
            ('marg_post', marg_post),
            ('rate_post', rate_post),
        ):
            f.create_dataset(k, data=np.asarray(arr, dtype=float),
                             compression='gzip', compression_opts=4)
        f.attrs['sigma_right'] = float(sigma_right)
        f.attrs['sigma_left'] = float(sigma_left)
        f.attrs['model_pars'] = json.dumps(_to_jsonable(model_pars))
        f.attrs['eff_mc'] = float(eff_mc)
        f.attrs['exposure'] = float(exposure)
        f.attrs['psa_mean'] = float(psa_mean) if psa_mean is not None else np.nan
        f.attrs['psa_sig'] = float(psa_sig) if psa_sig is not None else np.nan
        f.attrs['format_version'] = _FORMAT_VERSION

    # 2) DataFrames via pandas (appends to the same file)
    delta_table.to_hdf(p, key='delta_table', mode='a', format='table')
    stats.to_hdf(p, key='stats', mode='a', format='table')

    return p


def load_combined_result(path: str | Path) -> dict:
    """Read a combined result file.

    Returns
    -------
    dict with keys:
        ``nominal_bc_post``, ``marg_post``, ``rate_post``,
        ``delta_table`` (DataFrame), ``stats`` (DataFrame),
        ``sigma_right``, ``sigma_left``, ``model_pars``, ``eff_mc``,
        ``exposure``, ``psa_mean``, ``psa_sig``, ``format_version``.
    """
    p = Path(path)
    out: dict = {}
    with h5py.File(p, 'r') as f:
        out['nominal_bc_post'] = f['nominal_bc_post'][:]
        out['marg_post'] = f['marg_post'][:]
        out['rate_post'] = f['rate_post'][:]
        out['sigma_right'] = float(f.attrs['sigma_right'])
        out['sigma_left'] = float(f.attrs['sigma_left'])
        out['model_pars'] = json.loads(f.attrs.get('model_pars', '{}'))
        out['eff_mc'] = float(f.attrs['eff_mc'])
        out['exposure'] = float(f.attrs['exposure'])
        psa_mean = float(f.attrs['psa_mean'])
        psa_sig = float(f.attrs['psa_sig'])
        out['psa_mean'] = None if np.isnan(psa_mean) else psa_mean
        out['psa_sig'] = None if np.isnan(psa_sig) else psa_sig
        out['format_version'] = str(f.attrs.get('format_version', ''))

    out['delta_table'] = pd.read_hdf(p, key='delta_table')
    out['stats'] = pd.read_hdf(p, key='stats')
    return out


# ---------------------------------------------------------------------------
# SystResult I/O (format v2 — combined posteriors + per-component stats table)
# ---------------------------------------------------------------------------

_RESULT_FORMAT_VERSION = '2'


def post_stats(arr) -> dict:
    """Compute (mode, hdi_lo, hdi_hi, p16, p84, p90) for a 1-D posterior array.

    HDI is the 1σ (68.27%) highest-density interval (via ArviZ).  ``p16/p84/p90``
    are simple percentiles.  Returned values are plain Python floats so they can
    be JSON-serialised.
    """
    from qpym2.utils.mathy import calc_intervals
    a = np.asarray(arr).flatten()
    mode, lo, hi = calc_intervals([a], method='hdi')[0]
    p16, p84, p90 = np.percentile(a, [16, 84, 90])
    return {
        'mode':   float(mode),
        'hdi_lo': float(lo),
        'hdi_hi': float(hi),
        'p16':    float(p16),
        'p84':    float(p84),
        'p90':    float(p90),
    }


_post_stats = post_stats  # backward-compat alias


def save_syst_result(
    path: str | Path,
    *,
    name: str,
    bc_post: np.ndarray | None = None,
    raw_post: np.ndarray | None = None,
    model_pars: dict | None = None,
    per_post: pd.DataFrame | None = None,
    overwrite: bool = True,
) -> Path:
    """Save a per-systematic result file (format v2).

    Stats for ``bc_post`` and ``raw_post`` (when provided) are computed via
    :func:`post_stats` and stored in the root ``stats`` JSON attr under keys
    ``'bc'`` and ``'raw'``.  Either or both may be omitted (e.g. LSS-style
    studies where only ``per_post`` is meaningful).

    Parameters
    ----------
    path : str or Path
    name : str
        Short identifier (e.g. ``'bmmc'``, ``'lss'``).
    bc_post : 1-D array, optional
        Combined bias-corrected systematic posterior.
    raw_post : 1-D array, optional
        Combined raw (pre-BC) systematic posterior.
    model_pars : dict, optional
    per_post : DataFrame, optional
        Per-component stats for non-combined posteriors (per-ds, per-shift, ...).
        Conventionally has columns
        ``[name, bc, mode, hdi_lo, hdi_hi, p16, p84, p90, <param cols...>]``.
    overwrite : bool
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        if not overwrite:
            raise FileExistsError(p)
        p.unlink()

    stats: dict = {}
    if bc_post is not None:
        stats['bc'] = post_stats(bc_post)
    if raw_post is not None:
        stats['raw'] = post_stats(raw_post)

    with h5py.File(p, 'w') as f:
        if bc_post is not None:
            f.create_dataset('bc_post', data=np.asarray(bc_post, dtype=float),
                             compression='gzip', compression_opts=4)
        if raw_post is not None:
            f.create_dataset('raw_post', data=np.asarray(raw_post, dtype=float),
                             compression='gzip', compression_opts=4)
        f.attrs['name']           = name
        f.attrs['model_pars']     = json.dumps(_to_jsonable(model_pars or {}))
        f.attrs['stats']          = json.dumps(_to_jsonable(stats))
        f.attrs['format_version'] = _RESULT_FORMAT_VERSION

    if per_post is not None:
        per_post.to_hdf(p, key='per_post', mode='a', format='table')

    return p


def load_syst_result(path: str | Path) -> dict:
    """Read a file saved by :func:`save_syst_result`.

    Returns
    -------
    dict with keys:
        ``name`` (str), ``bc_post`` (np.ndarray), ``raw_post`` (np.ndarray or None),
        ``model_pars`` (dict), ``stats`` (dict with ``'bc'`` and optional ``'raw'``
        sub-dicts), ``per_post`` (DataFrame or None), ``format_version`` (str).
    """
    p = Path(path)
    out: dict = {}
    with h5py.File(p, 'r') as f:
        out['name']           = str(f.attrs['name'])
        out['bc_post']        = f['bc_post'][:]  if 'bc_post'  in f else None
        out['raw_post']       = f['raw_post'][:] if 'raw_post' in f else None
        out['model_pars']     = json.loads(f.attrs.get('model_pars', '{}'))
        out['stats']          = json.loads(f.attrs.get('stats', '{}'))
        out['format_version'] = str(f.attrs.get('format_version', ''))

    try:
        out['per_post'] = pd.read_hdf(p, key='per_post')
    except (KeyError, ValueError):
        out['per_post'] = None

    return out
