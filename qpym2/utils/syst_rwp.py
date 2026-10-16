import h5netcdf
import warnings
import numpy as np
import pandas as pd
import arviz as az

from qpym2.utils.logging import debug
from qpym2.utils.mathy import get_mode, sample_weighted_mixture
from qpym2.fit.bias import apply_bias_correction, read_bias

def _validate_fit(trace):
    _DIVPERC_CUTOFF = 5
    _RHAT_CUTOFF = 0.005

    divergent = trace.sample_stats['diverging'].values
    divperc = divergent.nonzero()[0].size / len(trace) * 100

    if divperc > _DIVPERC_CUTOFF:
        debug("Number of Divergent %d" % divergent.nonzero()[0].size)
        debug("Percentage of Divergent %.1f" % divperc)
        return False
    
    rhat = np.max(abs(az.rhat(trace) - 1).max().to_array()).item()
    if rhat > _RHAT_CUTOFF:
        debug(f'max rhat: {rhat}')
        return False
        # raise ValueError(f'rhat is too high: {rhat}')

    return True

def read_idata_nowarn(fname, base_group):
    # warnings.filterwarnings('ignore')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = az.InferenceData.from_netcdf(fname, base_group=base_group)

    # validate fit rhat, and divergences
    if not _validate_fit(data):
        return None
    
    return data

def read_systmc(fname, fitname_str='systmc', n=0):

    # read the existing systematic fits,
    with h5netcdf.File(fname, 'r') as f:
        fitnames = list(gn for gn in f.groups.keys() if fitname_str in gn)

    traces = read_traces(fname, fitname_str, n)
    # def _read_mcpost(sname):
    #     mcsysts = read_idata_nowarn(fname, base_group=sname)
    #     ndbds =  mcsysts.posterior.ndbd.values.flatten()
    #     maxlogp = get_mode(mcsysts.log_likelihood.pik.values.flatten())
    #     return ndbds, maxlogp
    
    posts =  [ tr.posterior.ndbd.values.flatten() for tr in traces ]
    maxlps = [ get_mode(tr.log_likelihood.pik.values.flatten()) for tr in traces ]

    # res = [ _read_mcpost(sname) for sname in fitnames]
    # posts, maxlps =  zip(*res)
    return posts, np.array(maxlps)

def read_traces(fname, fitname_str='systmc', n=0):
    # read the existing systematic fits,
    with h5netcdf.File(fname, 'r') as f:
        fitnames = list(gn for gn in f.groups.keys() if fitname_str in gn)

    if n < 0:
        raise ValueError('n must be positive')
    # elif n == 0: is ok, all values?
    else:
        if n > len(fitnames):
            debug(f'n is larger than the number of fits, setting n to {len(fitnames)}')
            n = len(fitnames)
            
    if n > 0:
        traces = [ read_idata_nowarn(fname, base_group=sname) for sname in fitnames[:n] ]
    else:
        traces = [ read_idata_nowarn(fname, base_group=sname) for sname in fitnames ]

    traces_filt = [tr for tr in traces if tr is not None]

    debug(f'{len(traces_filt)} traces read. Discarded {len(traces) - len(traces_filt)} due to bad fits.')
    return traces_filt

def calc_weights(idatas, lp_var='pik', estimator=None, prior_log_weights=None):
    """Compute per-fit mixture weights from a list of InferenceData.

    Computes per-fit scores from the log-likelihood (``log_likelihood[lp_var]``)
    via the given estimator, applies softmax + draw-count balancing, and returns
    the normalised weights.

    Parameters
    ----------
    idatas : list of arviz.InferenceData
        One trace per fit.
    lp_var : str, optional
        Log-likelihood variable used for weighting (default ``'pik'``).
    estimator : {None, 'mean', 'mode'}, optional
        Scoring rule applied to ``log_likelihood[lp_var]``. ``None`` → equal
        weights.
    prior_log_weights : array-like of float, optional
        Per-fit log adjustments added to the scores before softmax.

    Returns
    -------
    weights : np.ndarray
        Normalised per-fit weights (length ``len(idatas)``).
    """
    if estimator is None:
        return np.full(len(idatas), 1.0 / len(idatas))
    elif estimator == 'mean':
        scores = np.array(
            [idata.log_likelihood[lp_var].values.flatten().mean() for idata in idatas]
        )
    elif estimator == 'mode':
        scores = np.array(
            [get_mode(idata.log_likelihood[lp_var].values.flatten()) for idata in idatas]
        )
    else:
        raise ValueError(f"Invalid estimator: {estimator}")

    if prior_log_weights is not None:
        scores = scores + np.asarray(prior_log_weights, dtype=float)

    weights = np.exp(scores - scores.max())
    weights[weights < 1e-10] = 0.0
    # Equalize by draw count: scale by N_avg/n_i so fits with fewer MCMC samples
    # are not under-represented. Using N_avg keeps the numbers O(1).
    n_per_fit = np.array(
        [idata.posterior.dims['draw'] for idata in idatas], dtype=float
    )
    weights = weights * (n_per_fit.mean() / n_per_fit)
    weights = weights / weights.sum()

    return weights

def combine_posteriors(idatas, var='ndbd', lp_var='pik', estimator=None, prior_log_weights=None,
                       n_draws=None, rng=None):
    """Combine posteriors from multiple fits, weighted by log-likelihood mode.

    Thin wrapper around :func:`calc_weights` + :func:`sample_weighted_mixture`.

    Parameters
    ----------
    idatas : list of arviz.InferenceData
        One trace per fit.
    var : str, optional
        Posterior variable to combine (default ``'ndbd'``).
    lp_var : str, optional
        Log-likelihood variable used for weighting (default ``'pik'``).
    estimator : {None, 'mean', 'mode'}, optional
    prior_log_weights : array-like of float, optional
        Per-fit log adjustments added to the scores before softmax.
        After softmax, weights are further scaled by ``N_avg / n_i`` (where
        ``N_avg`` is the mean draw count across fits) to equalize the
        contribution of fits with different numbers of MCMC samples.
    n_draws : int, optional
        Number of resampled draws.  If ``None``, uses the minimum posterior
        size across inputs.
    rng : numpy.random.Generator, optional
        RNG instance for reproducibility.

    Returns
    -------
    dict with keys:
        ``samples``         – 1-D resampled array of length ``n_draws``.
        ``posteriors``      – list of 1-D arrays (one per fit).
        ``weights_per_fit`` – normalised per-fit weights.
    """
    posteriors = [idata.posterior[var].values.flatten() for idata in idatas]
    weights = calc_weights(
        idatas, lp_var=lp_var,
        estimator=estimator, prior_log_weights=prior_log_weights,
    )

    if n_draws is None:
        n_draws = int(np.min([len(p) for p in posteriors]))

    samples = sample_weighted_mixture(posteriors, weights, n_draws=n_draws, rng=rng)

    return dict(
        samples=samples,
        posteriors=posteriors,
        weights_per_fit=weights,
    )

# def correct_bias(fit_dir, model, nevs):
#     popt, pcov = read_bias(fit_dir, model, xmin=10)
#     dbd_corr = apply_bias_correction(nevs, popt)
#     return dbd_corr

def read_pca_eff_TMP():
    # WARN: reading from file not implemented -- which file?
    txt_data = \
        """3601 98.13 ± 0.06 98.41 ± 0.19 12329 ± 941 16371 ± 2637 4276 ± 591 6991 ± 1765
        3602 98.61 ± 0.08 98.80 ± 0.10 9923 ± 162 9762 ± 372 2842 ± 111 2899 ± 263
        3603 98.74 ± 0.03 98.85 ± 0.06 9809 ± 114 9909 ± 271 1843 ± 180 2233 ± 189
        3604 98.84 ± 0.02 98.81 ± 0.05 10389 ± 1107 10737 ± 490 1722 ± 1237 2301 ± 323
        3605 98.59 ± 0.08 98.87 ± 0.19 11282 ± 477 13020 ± 1402 4165 ± 347 5464 ± 1050
        3606 98.88 ± 0.03 98.61 ± 0.06 9928 ± 1434 10103 ± 322 1499 ± 2029 1892 ± 236
        3607 98.98 ± 0.04 98.75 ± 0.05 9525 ± 269 9471 ± 172 1514 ± 596 1796 ± 126
        3608 98.90 ± 0.04 98.64 ± 0.05 9491 ± 266 9324 ± 160 1531 ± 599 1746 ± 108
        3609 98.65 ± 0.03 98.57 ± 0.05 9262 ± 69 9223 ± 139 1473 ± 251 1599 ± 138
        3610 99.11 ± 0.03 98.91 ± 0.04 10039 ± 1154 10053 ± 258 1583 ± 1541 1948 ± 168
        3611 99.09 ± 0.04 99.01 ± 0.05 9306 ± 74 9143 ± 141 1584 ± 223 1733 ± 99
        3612 99.21 ± 0.03 99.12 ± 0.04 9634 ± 209 9581 ± 192 1591 ± 400 1784 ± 152
        3613 98.79 ± 0.03 98.76 ± 0.05 9659 ± 249 9675 ± 227 1604 ± 470 1919 ± 166
        3614 99.02 ± 0.04 99.11 ± 0.07 9692 ± 51 9516 ± 196 2643 ± 38 2606 ± 136
        3615 99.05 ± 0.04 98.93 ± 0.05 9548 ± 215 9401 ± 164 1603 ± 451 1781 ± 1161"""

    data_ = [ ln.split() for ln in txt_data.split('\n') ]
    col_names = ['ds',
                 'e_x2', 'err_e_x2', 'e_mcmc', 'err_e_mcmc',
                 'mu_x2', 'err_mu_x2', 'mu_mcmc', 'err_mu_mcmc',
                 's_x2', 'err_s_x2', 's_mcmc', 'err_s_mcmc']
    df = pd.DataFrame(data_)[[0,1,3,4,6,7,9,10,12,13,15,16,18]] # remove '+-' parts etc

    # clean up the read values
    df = df.apply(pd.to_numeric)
    df.rename(columns=dict(zip(df.columns, col_names)), inplace=True)
    df.set_index('ds', inplace=True)

    return df


def read_ds_exposures(fname):
    df = pd.DataFrame([l.split() for l in open(fname).readlines()])

    # convert 
    df = df.apply(pd.to_numeric).set_index(0).rename({1: 'exposure'}, axis=1)

    return df.rename_axis('ds')