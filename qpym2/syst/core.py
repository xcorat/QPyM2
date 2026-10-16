"""Core systematic delta computation and combination.

Provides the ``SystDelta`` container and the two main computation helpers:

- :func:`compute_delta` — combine a set of systematic traces (log-likelihood
  weighted), apply bias correction, and return the mode shift relative to a
  nominal bias-corrected posterior.
- :func:`combine_deltas` — quadrature-sum a collection of ``SystDelta`` into
  split-Gaussian widths ``(sigma_right, sigma_left)``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

import numpy as np

from qpym2.fit.bias import apply_bias_correction
from qpym2.utils.mathy import calc_intervals, get_mode
from qpym2.utils.syst_rwp import combine_posteriors

@dataclass
class SystDelta:
    """Container for one systematic contribution.

    Parameters
    ----------
    name : str
        Short identifier (e.g. ``'lss_dq'``, ``'bmmc'``, ``'ndbd_shift'``).
    delta_right : float
        Positive-side shift on the bias-corrected mode (>= 0 by convention).
    delta_left : float
        Negative-side shift on the bias-corrected mode (<= 0 by convention).
    model_pars : dict
        Reference model parameter dict (``co_cut``, ``ndbd_cutoff``,
        ``co_shift``) used to compute this delta.
    syst_post : np.ndarray, optional
        Bias-corrected systematic posterior samples (kept for plotting / I/O).
    stats : dict, optional
        Free-form metrics (e.g. ``mode_nom``, ``mode_syst``, ``p90_nom``,
        ``p90_syst``).
    """

    name: str
    delta_right: float
    delta_left: float
    model_pars: dict = field(default_factory=dict)
    syst_post: np.ndarray | None = None
    stats: dict = field(default_factory=dict)

def compute_delta_mode(
    nominal_bc_post: np.ndarray,
    syst_idatas: Iterable,
    bias_pars: tuple,
    var: str = 'ndbd',
    lp_var: str = 'pik',
    estimator: str | None = 'mode',
    n_draws: int | None = 50_000,
    rng: Any = None,
) -> tuple[float, np.ndarray]:
    # TODO: return the width scaling too (for weighing?)
    """Compute the bias-corrected mode shift for a set of systematic traces.

    Wraps the standard recipe used in ``06-bmmc_ndbd.ipynb``:

    1. Combine ``syst_idatas`` via :func:`combine_posteriors` (log-likelihood
       weighted resampling).
    2. Apply bias correction via :func:`apply_bias_correction` with ``bias_pars``.
    3. Compute ``get_mode(syst_bc) - get_mode(nominal_bc_post)``.

    Parameters
    ----------
    nominal_bc_post : np.ndarray
        Nominal bias-corrected posterior (already corrected by the caller).
    syst_idatas : iterable of arviz.InferenceData
        Systematic traces to combine.  A single ``InferenceData`` is also
        accepted (wrapped in a list).
    bias_pars : tuple
        ``(popt, pcov)`` from :func:`qpym2.fit.bias.read_bias_csv` or
        :func:`get_bias_from_fitdf`.
    var, lp_var, estimator, n_draws : see :func:`combine_posteriors`.
    rng : numpy.random.Generator, optional

    Returns
    -------
    delta : float
        ``mode(syst_bc) - mode(nominal_bc)``.
    syst_bc_post : np.ndarray
        Bias-corrected systematic posterior samples.
    """
    idatas = list(syst_idatas)
    if len(idatas) == 0:
        raise ValueError('compute_delta requires at least one syst trace')

    _rng = rng if rng is not None else np.random.default_rng()

    combined = combine_posteriors(
        idatas, var=var, lp_var=lp_var,
        estimator=estimator, n_draws=n_draws, rng=_rng,
    )
    popt, pcov = bias_pars
    syst_bc = apply_bias_correction(combined['samples'], popt)

    delta = float(get_mode(syst_bc) - get_mode(nominal_bc_post))
    return delta, syst_bc

def compute_deltas(p1, p2, method='hdi', return_extensions_only=False):
    """Compute the left and right deltas between two posteriors.
    
    if return_extensions_only is True, return the left and right extensions of the intervals,
    and return zero otherwise (p1_r > p2_r -> 0, p1_l < p2_l -> 0).

    """
    if len(p1) != len(p2):
        int1 = calc_intervals([p1], method=method)[0]
        int2 = calc_intervals([p2], method=method)[0]
    else:
        int1, int2 = calc_intervals([p1, p2], method=method)
    deltas = [int2[1] - int1[1], int2[2] - int1[2]]
    if return_extensions_only:
        deltas[0] = min(0, deltas[0])
        deltas[1] = max(0, deltas[1])
    print(f'Adjusted Deltas: {deltas}')
    return deltas

def compute_1sigma_edge_delta(
    nominal_post: np.ndarray,
    syst_posts: Iterable[np.ndarray] ,
    bias_pars: tuple | None = None
) -> tuple[float, np.ndarray]:
    """Compute the delta between the nominal mode and the 1-sigma edges of the 
        syst posteriors, return the total left/right shifts if multiple systs.

        - deltas added in quadrature

    """
    syst_posts = list(syst_posts)
    if len(syst_posts) == 0:
        raise ValueError('compute_1sigma_edge_delta requires at least one syst posterior')

    if bias_pars is not None:
        popt, pcov = bias_pars
        syst_posts = [apply_bias_correction(sp, popt) for sp in syst_posts]
        nominal_post = apply_bias_correction(nominal_post, popt)

    deltas = [compute_deltas(
        nominal_post, sp, return_extensions_only=True) for sp in syst_posts]
    print(deltas)
    return combine_deltas(deltas)

def combine_deltas(syst_deltas) -> tuple[float, float]:
    """Quadrature-sum systematic deltas into split-Gaussian widths. 

        sigma_right = sqrt( sum( d.delta_right**2 ) )
        sigma_left  = sqrt( sum( d.delta_left**2 ) )

    Each side is summed independently so that asymmetric systematics
    contribute only to the side on which they push the mode.

    Parameters
    ----------
    syst_deltas : iterable of tuples (left, right)

    Returns
    -------
    sigma_left : float
    sigma_right : float
        Both returned as positive magnitudes.
    """
    deltas = list(syst_deltas)
    print(f'Combining deltas: {deltas}')
    if len(deltas) == 0:
        return 0.0, 0.0
    sigma_left = float(np.sqrt(sum(d[0] ** 2 for d in deltas))) * -1.0
    sigma_right = float(np.sqrt(sum(d[1] ** 2 for d in deltas)))
    return sigma_left, sigma_right
