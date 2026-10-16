"""Marginalisation of the nominal posterior with a split-Gaussian smear.

KDE-resamples the (bias-corrected) nominal posterior and adds a one-sided
Gaussian on each side using the widths from :func:`combine_deltas`.
Extracted verbatim from ``06-syst-rate-new.ipynb``.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import gaussian_kde

def smear_posterior(
    nominal_bc_post: np.ndarray,
    sigma_right: float,
    sigma_left: float,
    n_smear: int = 200_000,
    rng_seed: int = 42,
    bw_method: str | float = 'scott',   
):
    """Extend nominal posterior below cutoff and add split-Gaussian smear.
    
    Combines extend_post (KDE extension in log-space) with split-Gaussian
    smearing applied to the extended samples.
    
    Parameters
    ----------
    nominal_bc_post : np.ndarray
        Nominal bias-corrected posterior samples.
    sigma_right, sigma_left : float
        Positive magnitudes of the right / left split-Gaussian arms.
    n_smear : int, optional
        Number of KDE samples to generate for extension.
    rng_seed : int, optional
        Seed for RNG.
    bw_method : str or float, optional
        Bandwidth for the KDE (default ``'scott'``).
        
    Returns
    -------
    smeared : np.ndarray
        Extended posterior + split-Gaussian smear.
    """

    extended = extend_skew_norm(nominal_bc_post,
                                upper_perc = 45,
                                lower_perc = 2,
                                fudge_ratio = 1)

    rng = np.random.default_rng(rng_seed)


    n = len(extended)
    u = rng.uniform(size=n)
    d = np.empty(n)
    m = u > 0.5
    d[m] = np.abs(rng.normal(0.0, sigma_right, int(m.sum())))
    d[~m] = -np.abs(rng.normal(0.0, sigma_left, int((~m).sum())))

    return  extended + d

def extend_post_kde(
    nominal_bc_post: np.ndarray,
    n_smear: int = 200_000,
    rng_seed: int = 42,
    bw_method: str | float = 'scott',
) -> np.ndarray:
    """KDE-extend nominal posterior below cutoff in log-space, preserving existing bins.

    Strategy:
    1. Calculate cutoff (2nd percentile) and log(cutoff)
    2. Transform data to log-space: log(nominal_bc - cutoff + epsilon)
    3. Create KDE in transformed space
    4. Sample from KDE
    5. Split samples into below/above log(cutoff), calculate ratio
    6. Resample KDE according to ratio, inverse-transform, and concatenate with original

    Parameters
    ----------
    nominal_bc_post : np.ndarray
        Nominal bias-corrected posterior samples.
    n_smear : int, optional
        Number of KDE samples to generate for extension.
    rng_seed : int, optional
        Seed for the KDE resampling.
    bw_method : str or float, optional
        Bandwidth for the KDE (default ``'scott'``).

    Returns
    -------
    extended_post : np.ndarray
        Original samples above cutoff + KDE-extended samples below cutoff.
    """
    # Step 1: Calculate cutoff (2nd percentile) and log(cutoff)
    cutoff = np.percentile(nominal_bc_post, 5)
    log_cutoff = np.log(cutoff)
    
    # Use a small epsilon to avoid log(0) when data touches cutoff
    # Use minimum to establish a proper baseline for log-space
    epsilon = 1e-8
    baseline = nominal_bc_post.min() - epsilon
    shifted_data = nominal_bc_post - baseline

    # Step 2: Transform data to log-space
    transformed_data = np.log(shifted_data)

    # Step 3: Create KDE in transformed space
    # kde = gaussian_kde(transformed_data, bw_method=bw_method)
    kde = gaussian_kde(nominal_bc_post, bw_method=bw_method)

    # Step 4: Sample from KDE in transformed space
    kde_samples_transformed = kde.resample(n_smear, seed=rng_seed)[0]
    rng = np.random.default_rng(rng_seed + 1)

    # Adjust log_cutoff for the baseline shift
    log_cutoff_shifted = np.log(cutoff - baseline)

    # Step 5: Split and calculate ratio of below to above log(cutoff)
    below_log_cutoff = kde_samples_transformed < log_cutoff_shifted
    n_below = below_log_cutoff.sum()
    n_above = (~below_log_cutoff).sum()
    ratio = n_below / n_above if n_above > 0 else 0

    # Step 6: Resample KDE according to ratio and inverse-transform
    nom_bc_above_cutoff = nominal_bc_post[nominal_bc_post >= cutoff]
    n_kde_to_use = int(len(nom_bc_above_cutoff) * ratio)

    if n_kde_to_use > 0:
        # Randomly select KDE samples from below log_cutoff region
        kde_below_indices = np.where(kde_samples_transformed < log_cutoff_shifted)[0]
        selected_indices = rng.choice(
            kde_below_indices,
            size=min(n_kde_to_use, len(kde_below_indices)),
            replace=False,
        )
        kde_selected_transformed = kde_samples_transformed[selected_indices]

        # Inverse transform: x = exp(transformed) + baseline
        kde_selected_original = np.exp(kde_selected_transformed) + baseline
    else:
        kde_selected_original = np.array([])

    # Step 7: Concatenate original data above cutoff with KDE-extended samples
    return np.concatenate([nom_bc_above_cutoff, kde_selected_original])

def extend_skew_norm(
    nom_post_bc,
    upper_perc = 45,
    lower_perc = 2,
    rng_seed: int = 42,
    fudge_ratio = 1
):
    from scipy.stats import skewnorm
    from scipy.stats import norm

    from scipy.optimize import curve_fit

    # trunc = np.percentile(nom_post_bc, 40)
    # nom_trunc = nom_post_bc[nom_post_bc > trunc]

    tr_up = np.percentile(nom_post_bc, upper_perc)
    tr_low = np.percentile(nom_post_bc, lower_perc)
    up_trunc = nom_post_bc[nom_post_bc < tr_up]
    trunc = up_trunc[up_trunc > tr_low]
    nom_trunc = nom_post_bc[nom_post_bc > tr_low] # truncated from below, remove a little from the low edge


    histl_c, binl_e = np.histogram(up_trunc, bins=300, density=True)
    binl_c = (binl_e[:-1] + binl_e[1:]) / 2

    poptl, _ = curve_fit(lambda x, loc, scale, alpha: skewnorm.pdf(x, alpha, loc=loc, scale=scale),
                        binl_c, histl_c, p0=[nom_post_bc.mean(), 4, 0.0], maxfev=5000)


    loc_fit, scale_fit, alpha_fit = poptl
    tmp_std = nom_post_bc.std()

    xmin = nom_post_bc.mean() - 4*tmp_std
    # x_fit = np.linspace(xmin, x0.max(), 200)

    # x_fit = np.linspace(x0.min(), x0.max(), 200)
    # skspdf =  skewnorm.pdf(x_fit, alpha_fit, loc=loc_fit, scale=scale_fit) / len(nom_post_bc) * len(up_trunc)

    cdf_t = skewnorm.cdf(tr_up, alpha_fit, loc=loc_fit, scale=scale_fit)
    cdf_m = skewnorm.cdf(tr_low, alpha_fit, loc=loc_fit, scale=scale_fit)
    cdf_e = skewnorm.cdf(xmin, alpha_fit, loc=loc_fit, scale=scale_fit)

    sknowm_ratio = (cdf_m - cdf_e) / (cdf_t - cdf_e) * fudge_ratio
    n_ext = int(len(trunc) * sknowm_ratio)

    rng = np.random.RandomState(rng_seed)
    u = rng.uniform(cdf_e, cdf_m, n_ext)
    ext_samples = skewnorm.ppf(u, alpha_fit, loc=loc_fit, scale=scale_fit)
    return np.concatenate([nom_trunc, ext_samples])

def smear_posterior_old(
    nominal_bc_post: np.ndarray,
    sigma_right: float,
    sigma_left: float,
    n_smear: int = 200_000,
    rng_seed: int = 42,
    bw_method: str | float = 'scott',
) -> np.ndarray:
    """KDE-resample ``nominal_bc_post`` and add a split-Gaussian smear.

    The KDE step lets the smeared distribution extend naturally below zero
    so the smear is symmetric in samples-space rather than density-space.
    Caller is responsible for clipping to ``[0, inf)`` for physics use.

    Parameters
    ----------
    nominal_bc_post : np.ndarray
        Nominal bias-corrected posterior samples.
    sigma_right, sigma_left : float
        Positive magnitudes of the right / left split-Gaussian arms
        (from :func:`qpym2.syst.combine_deltas`).
    n_smear : int, optional
        Number of output draws.
    rng_seed : int, optional
        Seed for the smearing RNG (KDE.resample uses ``seed=rng_seed``,
        smear noise uses ``seed=rng_seed + 1``).
    bw_method : str or float, optional
        Bandwidth for the KDE (default ``'scott'``).

    Returns
    -------
    smeared : np.ndarray, shape (n_smear,)
        Unclipped smeared samples.
    """
    kde = gaussian_kde(nominal_bc_post, bw_method=bw_method)
    kde_samples = kde.resample(n_smear, seed=rng_seed)[0]

    rng = np.random.default_rng(rng_seed + 1)
    n = len(kde_samples)
    u = rng.uniform(size=n)
    d = np.empty(n)
    m = u > 0.5
    d[m] = np.abs(rng.normal(0.0, sigma_right, int(m.sum())))
    d[~m] = -np.abs(rng.normal(0.0, sigma_left, int((~m).sum())))

    return kde_samples + d
