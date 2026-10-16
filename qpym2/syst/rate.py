"""Rate conversion with optional PSA-efficiency marginalisation.

Constants ``F130`` and ``N_TEO2`` mirror those in
:mod:`qpym2.fit.bias` (internal note 123D).
"""

from __future__ import annotations

import arviz as az
import numpy as np

from qpym2.utils.mathy import get_mode


F130 = 0.34167                          # ¹³⁰Te isotopic fraction (internal note 123D)
N_TEO2 = 6.022e23 * 1000 / 159.6        # TeO₂ molecules per kg


def to_rate(
    smeared_post: np.ndarray,
    eff_mc: float,
    exposure: float,
    psa_mean: float | None = None,
    psa_sig: float | None = None,
    rng_seed: int = 44,
    clip: bool = True,
) -> np.ndarray:
    """Convert a (smeared) ndbd counts posterior to a decay rate posterior.

    ``rate = n / (F130 · N_TeO2 · exposure · eff_MC · eff_PSA)``

    If both ``psa_mean`` and ``psa_sig`` are given, ``eff_PSA`` is sampled
    per draw from ``N(psa_mean, psa_sig)``.  Otherwise ``eff_PSA = psa_mean``
    if provided, else 1.

    Parameters
    ----------
    smeared_post : np.ndarray
        Marginalised (split-Gaussian smeared) ndbd counts posterior.
    eff_mc : float
        MC signal selection efficiency (from
        :func:`qpym2.fit.results.read_ndbd_eff`).
    exposure : float
        Exposure in kg·yr.
    psa_mean, psa_sig : float, optional
        PSA efficiency mean and 1σ.  Set both to enable PSA marginalisation.
    rng_seed : int, optional
        RNG seed for PSA sampling.
    clip : bool, optional
        If ``True`` (default), clip ``smeared_post`` to ``>= 0`` first.

    Returns
    -------
    rate : np.ndarray
        Decay rate samples [yr⁻¹], same length as the (clipped) input.
    """
    arr = np.asarray(smeared_post, dtype=float)
    if clip:
        arr = arr[arr >= 0]

    if psa_mean is not None and psa_sig is not None:
        rng = np.random.default_rng(rng_seed)
        eff_psa = rng.normal(psa_mean, psa_sig, len(arr))
    elif psa_mean is not None:
        eff_psa = psa_mean
    else:
        eff_psa = 1.0

    return arr / (F130 * N_TEO2 * exposure * eff_mc * eff_psa)


def rate_stats(rate: np.ndarray, hdi_prob: float = 0.68) -> dict:
    """Summary stats for a rate posterior.

    TODO: DO NOT USE! use the fns in mathy.

    Returns
    -------
    dict with keys ``mode``, ``p90``, ``hdi`` (tuple of (lo, hi)),
    ``hdi_prob`` (input echoed).
    """
    print("rate_stats is deprecated; use mathy.summary_stats instead")
    rate = np.asarray(rate)
    hdi = az.hdi(rate, hdi_prob=hdi_prob)
    return {
        'mode': float(get_mode(rate)),
        'p90': float(np.percentile(rate, 90)),
        'hdi': (float(hdi[0]), float(hdi[1])),
        'hdi_prob': float(hdi_prob),
    }
