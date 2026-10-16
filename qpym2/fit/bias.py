from qpym2.fit.results import read_toy_df, read_toy_df_csv
def get_bias_from_fitdf(df, xmin=13, ax=None, plot_pars=None):
    df = df[df['ndbd_inj'] > xmin]  
    x = df['ndbd_inj'].values
    y = df[('ndbd', 'mode')].values
    y_uerr = df[('ndbd', 'upp34')].values - y


    # fit the plot with polynomial 1st order
    from scipy.optimize import curve_fit
    def pol1(x, a, b):
        return a * x + b

    popt, pcov = curve_fit(pol1, x, y, sigma=y_uerr) #, absolute_sigma=True)

    if ax is not None:
        name = 'bias'
        if plot_pars:
            name = plot_pars.get('name', '')

        ax.errorbar(x, y, yerr=y_uerr, fmt='none', capsize=3)
        ax.plot(x, pol1(x, *popt), 'r-')
        ax.set_xlabel('n_injected')
        ax.set_ylabel('fitted')
        ax.set_title(name)
        ax.text(0.05, 0.95, f'a = {popt[0]:.3f} +/- {pcov[0,0]**0.5:.3f}\nb = {popt[1]:.3f} +/- {pcov[1,1]**0.5:.3f}',
                transform=ax.transAxes, va='top', ha='left')

    return popt, pcov

def read_bias(fit_dir, model_pars, xmin=None):
    tdf = read_toy_df( fit_dir, model_pars, type='bias')
    popt, pcov = get_bias_from_fitdf(tdf, xmin=xmin)

    return popt, pcov

def read_bias_csv(csv_path, model_pars, xmin=None):
    """Read bias parameters from a batch CSV file.

    Like :func:`read_bias` but reads from a batch CSV instead of HDF5.
    """
    tdf = read_toy_df_csv(csv_path, model_pars, type='bias')
    popt, pcov = get_bias_from_fitdf(tdf, xmin=xmin)
    return popt, pcov

import numpy as np

def apply_bias_correction(x, popt):
    """Deterministic bias correction: corrected = (x - b) / a.

    This is the standard function to use when bias-correcting a posterior.
    Uses the central fitted parameters ``popt`` only — no sampling.

    Parameters
    ----------
    x : array-like
        Posterior samples to correct.
    popt : array-like
        Fitted bias parameters ``[a, b]`` (``fitted = a * true + b``).

    Returns
    -------
    corrected : np.ndarray
    """
    x_ = np.asarray(x, dtype=float)
    a, b = popt
    # NOTE: a ~ 1, so we don't expect to divide by zero.
    return (x_ - b) / a


def sample_corrected(x, popt, pcov, rng=None):
    """Bias correction with parameter uncertainty propagation (for syst studies).

    Draws a different ``(a, b)`` pair from ``N(popt, pcov)`` for *each*
    posterior sample, propagating bias-fit parameter uncertainty into the
    corrected posterior.  
    
    For normal bias correction use :func:`apply_bias_correction` instead.
    """
    x_ = np.array(x, dtype=float)
    ns = x_.shape
    _rng = rng if rng is not None else np.random.default_rng()
    coeff_samples = _rng.multivariate_normal(mean=popt, cov=pcov, size=x_.size)
    a = coeff_samples[:, 0].reshape(ns)
    b = coeff_samples[:, 1].reshape(ns)
    tiny = 1e-12
    a = np.where(np.abs(a) < tiny, np.sign(a) * tiny, a)
    return (x_ - b) / a

def sample_corrected_rate(x, popt, pcov, eff_MC,
                          exposure, eff_pca_gaus=None, rng=None):
    x_ = np.array(x, dtype=float)
    ns = x_.shape
    _rng = rng if rng is not None else np.random.default_rng()
    coeff_samples = _rng.multivariate_normal(mean=popt, cov=pcov, size=x_.size)
    a = coeff_samples[:, 0].reshape(ns)
    b = coeff_samples[:, 1].reshape(ns)
    tiny = 1e-12
    a = np.where(np.abs(a) < tiny, np.sign(a) * tiny, a)
    y = (x_ - b) / a

    f130 = 0.34167 # TODO: From internal note 123D (convert to nnusance param?)
    n_teO2 = 6.022e23*1000/159.6 # From internal note 123D

    efficiency = 1
    if eff_pca_gaus:
        efficiency = np.random.normal(eff_pca_gaus[0], eff_pca_gaus[1], y.shape)
        
    efficiency *= eff_MC

    rate = y/(f130*n_teO2)/(exposure*efficiency)

    return rate

def sample_rate(nevs, eff_MC, exposure, eff_pca_gaus=None):
    f130 = 0.34167 # TODO: From internal note 123D (convert to nnusance param?)
    n_teO2 = 6.022e23*1000/159.6 # From internal note 123D

    efficiency = 1
    if eff_pca_gaus:
        efficiency = np.random.normal(eff_pca_gaus[0], eff_pca_gaus[1], y.shape)
        
    efficiency *= eff_MC

    rate = nevs/(f130*n_teO2)/(exposure*efficiency)

    return rate

def get_weighted_hist(samples, maxlps, bins=100):
    """.. deprecated:: Use ``syst_rwp.combine_posteriors`` instead.
       This function averages histograms (lossy). ``combine_posteriors``
       operates on weighted samples and preserves full information.
    """
    # we should align the bins before taking the weighted average
    if isinstance( bins, int):
        # TODO: calculate the range automatically
        min, max = samples.min(), maxlps.max()
        bins = np.linspace(min, max, bins+1)
    hist_out = [np.histogram(vals, bins=bins, density=True) for vals in samples]
    hists = np.array([h[0] for h in hist_out])
    weights = np.exp((maxlps - maxlps.max()).astype(np.float64))
    weights = weights/weights.sum()
    return np.average(hists, weights=weights, axis=0)