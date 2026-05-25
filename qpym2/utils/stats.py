"""
Implement math related utility functions here.

The main functionality implemented here are the get mode -- which uses a fancy way to 
compute the mode when the numbers are too small.
"""
import numpy as np
import pandas as pd
import arviz as az

from scipy.stats import gaussian_kde
from scipy.optimize import minimize_scalar

def get_mode(arr, bw_method=None):
    """ Get the mode of the array `arr` using a gaussian kernel density estimation. 
    
    Parameters
    ----------
    arr : array-like
        The array to find the mode of.
    bw_method : str, scalar or callable, optional
        The method used to calculate the estimator bandwidth. This can be 'scott', 'silverman',
        a scalar constant or a callable. See `scipy.stats.gaussian_kde` for more details. 
        If None, the bandwidth is estimated using the scott's rule of thumb.
            `bw` = `n`**(-1./(arr.ndim+4))
    """
    # TODO: error checking and performance!!
    arr = np.array(arr)
    mean = arr.mean()
    marr = arr/mean - 1

    if not bw_method:
        bw_method = marr.size**(-1./(marr.ndim+4))
    kde = gaussian_kde(marr, bw_method=bw_method)
    xmax = minimize_scalar( lambda x: -kde.evaluate(x), bounds=(marr.min(), marr.max()), method='bounded').x.item()

    # TODO: check with other methods
    # print(mean, (xmax + 1)*mean)
    return (xmax + 1)*mean

# def get_1sig_around_mode(arr, bw_method=None, method='hdi'):
#     """Return (mode, low, high) 1σ interval around the KDE mode of `arr`.

#     Parameters
#     ----------
#     arr : array-like
#         Posterior samples.
#     bw_method : str, scalar or callable, optional
#         Bandwidth method passed to `get_mode`. See `get_mode` for details.
#     method : {'hdi', 'half_interval'}, optional
#         How to compute the interval bounds.

#         ``'hdi'`` (default)
#             68.27% highest-density interval via ``az.hdi``.  The HDI is the
#             shortest interval containing 68.27% of the mass.  The mode and
#             interval bounds are computed independently.

#         ``'half_interval'``
#             Split samples at the mode and take the 68.27th percentile of each
#             half independently.  For a symmetric (Gaussian) posterior this
#             recovers mode ± σ exactly.  Appropriate for asymmetric posteriors
#             when the mode is used as the central value:

#             - left arm  = 31.73rd percentile of samples ≤ mode
#               (= (1 − 0.6827) × 100 % of the left half)
#             - right arm = 68.27th percentile of samples ≥ mode

#     Returns
#     -------
#     mode : float
#     low : float
#         Lower 1σ bound.
#     high : float
#         Upper 1σ bound.
#     """
#     arr = np.asarray(arr)
#     mode = get_mode(arr, bw_method=bw_method)
#     if method == 'hdi':
#         low, high = az.hdi(arr, hdi_prob=0.6827)
#     elif method == 'half_interval':
#         left  = arr[arr <= mode]
#         right = arr[arr >= mode]
#         low   = np.percentile(left,  (1.0 - 0.6827) * 100)  # 31.73rd pctile of left half
#         high  = np.percentile(right,         0.6827  * 100)  # 68.27th pctile of right half
#     elif method == 'stddev':
#         # For a Gaussian, the 1σ interval is mean ± stddev.  Using the mode as the central value,
#         # we can compute the upper and lower bounds as mode ± stddev.
#         stddev = np.std(arr)
#         low = mode - stddev
#         high = mode + stddev
#     elif method == 'percentiles':
#         # This method is not centered around the mode, but simply takes the 15.87th and 84.13th percentiles,
#         # which correspond to the mean ± 1σ for a Gaussian distribution.
#         low = np.percentile(arr, 15.87)
#         high = np.percentile(arr, 84.13)
#     else:
#         raise ValueError(f"method must be 'hdi' or 'half_interval', got {method!r}")
#     return mode, low, high

# def calc_intervals(arrays, method='hdi', bw_method=None):
#     """Compute (mode, low, high) 1σ intervals for a list of equal-length posterior arrays.

#     Bounds are computed in a single batched call; mode estimation (KDE) is per-array.

#     Parameters
#     ----------
#     arrays : list of array-like
#         All arrays must have the same length.
#     method : {'hdi', 'half_interval', 'percentiles', 'stddev'}
#     bw_method : optional
#         Bandwidth for KDE mode estimation.

#     Returns
#     -------
#     list of (mode, low, high) tuples
#     """
#     arrays = [np.asarray(a).flatten() for a in arrays]
#     lengths = [len(a) for a in arrays]
#     assert len(set(lengths)) == 1, f"All arrays must have the same length, got {set(lengths)}"

#     modes = [get_mode(a, bw_method=bw_method) for a in arrays]
#     stacked = np.stack(arrays)   # (n_arrays, n_samples)

#     if method == 'hdi':
#         # az.hdi treats axis 0 as the sample axis, so pass each 1-D array separately
#         hdis = np.array([az.hdi(a, hdi_prob=0.6827) for a in arrays])   # (n_arrays, 2)
#         return list(zip(modes, hdis[:, 0].tolist(), hdis[:, 1].tolist()))
#     elif method == 'percentiles':
#         pcts = np.percentile(stacked, [15.87, 84.13], axis=1)   # (2, n_arrays)
#         return list(zip(modes, pcts[0].tolist(), pcts[1].tolist()))
#     elif method == 'stddev':
#         stds = np.std(stacked, axis=1)
#         return [(m, m - s, m + s) for m, s in zip(modes, stds.tolist())]
#     elif method == 'half_interval':
#         results = []
#         for arr, mode in zip(arrays, modes):
#             left = arr[arr <= mode]
#             right = arr[arr >= mode]
#             low   = float(np.percentile(left,  (1.0 - 0.6827) * 100))
#             high  = float(np.percentile(right,        0.6827  * 100))
#             results.append((mode, low, high))
#         return results
#     else:
#         raise ValueError(f"Unknown method {method!r}")

# def find_mode(trace, varnames=None):
#     """ Find the mode of the variables in `varnames` from the `trace` object.
    
#     Parameters
#     ----------
#     trace : arviz.InferenceData
#         The trace object containing posterior samples.
#     varnames : list of str, optional
#         Variable names to extract modes for, in the desired order.
#         If None, defaults to trace.posterior.data_vars (ALPHABETICAL order from xarray).
        
#     Returns
#     -------
#     list
#         Mode values in the same order as `varnames`.
        
#     Warning
#     -------
#     CRITICAL: xarray stores variables in ALPHABETICAL order by default.
#     If varnames is omitted, returned values will be alphabetically ordered,
#     which may NOT match your component ordering (e.g., comps.index).
    
#     Always pass explicit varnames=list(comps.index) to avoid component value swapping:
#         ✅ CORRECT: find_mode(trace, varnames=list(comps.index))
#         ❌ WRONG:   find_mode(trace)  # Returns alphabetical order!
#     """
#     if varnames is None:
#         varnames = trace.posterior.data_vars
#     return [ get_mode(trace.posterior[varname].values.flatten()) for varname in varnames ]

# def sample_weighted_mixture(posteriors, weights, n_draws=None, rng=None):
#     """Draw independent samples from a weighted mixture of posteriors.

#     Each posterior array is sampled proportionally to its weight, producing
#     a single combined sample set suitable for mode/HDI/percentile estimation.

#     Parameters
#     ----------
#     posteriors : list of array-like
#         One 1-D array of posterior samples per fit.
#     weights : array-like of float
#         Normalised per-fit weights (same length as ``posteriors``).
#     n_draws : int, optional
#         Total number of draws (default 10,000).
#     rng : numpy.random.Generator, optional
#         RNG instance.  Created with default seed if not given.

#     Returns
#     -------
#     samples : np.ndarray
#         1-D array of ``n_draws`` independent samples from the mixture.
#     """
#     if rng is None:
#         rng = np.random.default_rng()
#     weights = np.asarray(weights, dtype=float)
#     weights = weights / weights.sum()

#     arrays = [np.asarray(p).flatten() for p in posteriors]
#     if n_draws is None:
#         n_draws = min(len(a) for a in arrays)

#     # Assign each individual sample a weight = mixture_weight / posterior_size.
#     # This ensures every posterior contributes proportionally to its mixture weight
#     # regardless of how many samples it contains, so the method works correctly
#     # even when n_draws >> total posterior samples.
#     per_sample_w = np.concatenate(
#         [np.full(len(a), w / len(a)) for a, w in zip(arrays, weights)]
#     )
#     per_sample_w /= per_sample_w.sum()
#     pool = np.concatenate(arrays)

#     idx = rng.choice(len(pool), size=n_draws, p=per_sample_w, replace=True)
#     return pool[idx]

# def signed_quad_sum(deltas):
#     """Quadrature sum of deltas, split into positive and negative contributions.

#     Parameters
#     ----------
#     deltas : array-like
#         Per-source shifts (e.g. combined_mode − nominal_mode).

#     Returns
#     -------
#     sigma_plus : float
#         sqrt(sum(d² for d > 0))
#     sigma_minus : float
#         sqrt(sum(d² for d < 0))  — returned as a positive number.
#     """
#     deltas = np.asarray(deltas, dtype=float)
#     pos = deltas[deltas > 0]
#     neg = deltas[deltas < 0]
#     return np.sqrt(np.sum(pos**2)), np.sqrt(np.sum(neg**2))

# def mergeto_random(to_df, from_df):
#     """ Merge the two dataframes while picking randomly spaces from the `to_df`
#     """
#     len_total = len(to_df)
#     rnd_indices = np.random.randint(0, len_total, len(from_df)) # from_df.iloc[:,0]
#     rnd_indices.sort()
    
#     slices = []
#     last_index = 0
#     for i in range(len(from_df)):
#         # print(i, to_df.iloc[last_index])
#         slices.append(to_df.iloc[last_index: rnd_indices[i]])
#         slices.append( from_df[i:i+1] )
#         last_index = rnd_indices[i] 
    
#     slices.append(to_df.iloc[last_index: ])
        
#     return pd.concat(slices).reset_index(drop=True)

# def test_mode(arr, ax, **kwargs):
#     """ Test the mode finding method. """
    
#     from scipy.stats import gaussian_kde
#     from scipy.optimize import minimize_scalar
    
#     raise NotImplementedError()
#     mode = get_mode(arr)
#     arr = np.array(arr)
#     mean = arr.mean()
#     marr = arr#/mean - 1

#     bw_method = marr.size**(-1./(marr.ndim+4))
#     kwa_ = dict(bw_method=bw_method)
#     kwa_.update(kwargs)
    
#     kde = gaussian_kde(marr, **kwa_)
#     xmax = minimize_scalar( lambda x: -kde.evaluate(x), bounds=(marr.min(), marr.max()), method='bounded').x.item()

#     ax.hist(marr, bins=50, density=True)
#     ax.axvline(x=xmax, color='red')

#     ls = np.linspace(marr.min(), marr.max(), 100)
#     ax.plot(ls, kde.evaluate(ls))

#     return ax

# import math
# def _sant_ceil(val, pow10=5):
#     """ return ceil, but ignore the digits beyond `pow10` """
#     v = math.floor(val*math.pow(10, pow10))
#     return math.ceil(v/math.pow(10, pow10))

