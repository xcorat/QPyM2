""" Functions needed for plots related to the posterior distributions. """ 

import numpy as np
from matplotlib import pyplot as plt

from qpym2 import utils
from qpym2 import hists as hist_tools

def plot_posterior(trace, varname='ctest_ndbd', ax=None, **kwargs):
    """ Plot the posterior distribution for a given variable. 
    
    Parameters
    ----------
    trace : pymc3 trace
        The trace from the fit.
    varname : str
        The name of the variable to plot.
    ax : matplotlib axis
        The axis to plot on. If None, a new figure is created.
    **kwargs : dict
        Additional arguments to pass to the histogram function.

    Returns
    -------
    ax : matplotlib axis
        The axis with the plot.

    TODO: try this maybe?
        b = trace.posterior['cnoco2_no0v1740_v150_smooth_ndbd'].values.flatten()
        az.plot_dist(b,  quantiles=[.16, .5, .84])
    """


    n0v_p = trace.posterior[varname].values.flatten()
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        
    hist_kwargs = {'bins': 50, 'color': '#0504aa', 'alpha': 0.7}
    ax.hist(x=n0v_p, **hist_kwargs, **kwargs)
    ax.set_yscale('log')

    mode = utils.get_mode(n0v_p)
    p90 = np.percentile(n0v_p, 90)
    p_up1s = np.percentile(n0v_p[n0v_p>=mode], 68)
    p_low1s = np.percentile(n0v_p[n0v_p<=mode], 32)
    # mode, p90, p_up1s, p_down1s
    ax.axvline(x=mode, color='black', label=f'Mode: {mode:.2f}')
    ax.axvline(x=p90, color='red', label=f'90th perc: {p90:.2f}')
    ax.axvline(x=p_up1s, color='brown', label=f'Upper 1sig: {p_up1s:.2f}')
    ax.axvline(x=p_low1s, color='grey', label=f'Lower 1sig: {p_low1s:.2f}')
    ax.legend()

    return ax

def plot_fit_proj_numpy(table, norms, data, hm, ax=None, **kwargs):
    nbins, range = hist_tools.get_hist_settings(hm)
    h2_data, xedges, yedges = np.histogram2d(data.T[0], data.T[1], bins=nbins, range=range)
    h2_fit = hist_tools.get_sum_hist(table, norms)

    hu_data = hist_tools.get_h2proj(h2_data, axis=1)
    hu_fit = hist_tools.get_h2proj(h2_fit, axis=1)
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # return (xedges[:-1], hu_data, np.sqrt(hu_data))
    ax.errorbar(xedges[:-1], hu_data, yerr=np.sqrt(hu_data), fmt='o', label='data')
    ax.plot(xedges[:-1], hu_fit, label='fit')
    ax.set_yscale('log')
    ax.legend()

    return ax

def plot_fit_proj(table, norms, data, hm, ax=None, backend='numpy', **kwargs):
    """ TODO TODO """
    if backend == 'numpy':
        return plot_fit_proj_numpy(table, norms, data, hm, ax=ax, **kwargs)
    else:
        return None
    th2d_fit = mcroot.get_empty_hist(hm, return_numpy=False)

    cfg_name = 'ushiftp05'

    run_config =  global_vars.run_configs[run_name]
    shift_cfg = global_vars.shifts[cfg_name]

    mc_path = f'{global_vars.datapath}/MC/mcuvs_roi310' #{shift_cfg["fname_postfix"]}'
    data_fname = f'{mc_path}/data_unblinded.root'

    cuts = run_config['cuts']

    th2d_data = mcroot.create_hist_from_staged(data_fname, global_vars.hm, cuts=f"{cuts['run']}&&{cuts['data']}", return_numpy=False)

    th2d_fit = th2d_data.Clone('th2d_fit')
    th2d_fit.Reset('RSEM')
    for i in range(h2_fit.shape[0]):
        for j in range(h2_fit.shape[1]):
            th2d_fit.SetBinContent(i+1, j+1, h2_fit[i, j])