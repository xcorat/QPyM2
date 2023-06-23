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
    nbins, range_ = hist_tools.get_hist_settings(hm)
    h2_data, xedges, yedges = np.histogram2d(data.T[0], data.T[1], bins=nbins, range=range_)
    h2_fit = hist_tools.get_sum_hist(table, norms)

    hu_data = hist_tools.get_h2proj(h2_data, axis=1)
    hu_fit = hist_tools.get_h2proj(h2_fit, axis=1)
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # return (xedges[:-1], hu_data, np.sqrt(hu_data))
    ax.errorbar(xedges[:-1], hu_data, yerr=np.sqrt(hu_data), fmt='.', label='data')
    ax.plot(xedges[:-1], hu_fit, label='fit')
    ax.set_yscale('log')
    ax.legend()

    return ax

def plot_fit_proj_root(table, norms, data, hm, ax=None, **kwargs):
    """ Plot the fit projection using ROOT. "
    
    TODO: does axis=1 X axis or Y axis?
    """
    from ROOT import TCanvas, TH2D
    from qpym2.plotting.legacy import H1ComparisonPlotter

    can = ax
    axis = kwargs.get('axis', 1)

    h2_fit = hist_tools.get_sum_hist(table, norms)
    nbins, range_ = hist_tools.get_hist_settings(hm)
    h2_data, xedges, yedges = np.histogram2d(data.T[0], data.T[1], bins=nbins, range=range_)

    if can is None:
        can = TCanvas()

    # TODO: 
    # 1. maybe make a copy? (deepcopy).. 
    # 2. or change the params so we can change them individually
    # maybe implement **kwargs to be copied into the cfg dict?
    if 'hname' in kwargs:
        hname = kwargs['hname']
    else:
        hname = hm.name
    hm.name = hname + '_fit'
    th2d_fit = hist_tools.get_empty_hist(hm, return_numpy=False)
    hm.name = hname + '_data'
    th2d_data = hist_tools.get_empty_hist(hm, return_numpy=False)
    hm.name = hname

    for i in range(h2_fit.shape[0]):
        for j in range(h2_fit.shape[1]):
            th2d_fit.SetBinContent(i+1, j+1, h2_fit[i, j])
            th2d_data.SetBinContent(i+1, j+1, h2_data[i, j])

    if axis == 1:
        th1d_fit = th2d_fit.ProjectionX()
        th1d_data = th2d_data.ProjectionX()
    elif axis == 0:
        th1d_fit = th2d_fit.ProjectionY()
        th1d_data = th2d_data.ProjectionY()

    if 'rebin' in kwargs:
        th1d_fit.Rebin(kwargs['rebin'])
        th1d_data.Rebin(kwargs['rebin'])

    plot_cfg = {
        'log_yaxis_hist': True,
        'out_fpath': "",
        'plot_title': f"Data vs Fit projection ({hname})",
        'xaxis_title': "ESum/#sqrt{2} [keV]",
        'yaxis_title': "Counts/keV",
        'h1_title': "Data",
        'h2_title': "Fit",
        'fit': None,
        # 'hist_stats': 'nie'
        'out_fpath': 'test.root',
    }
    plotter = H1ComparisonPlotter(th1d_data, th1d_fit, plot_cfg, no_errorbars_on2=True)
    if 'cfg' in kwargs:
        plotter.set_plot_config(kwargs['cfg'])

    plotter.draw_comparison_plot(can)

    return plotter

def plot_fit_proj(table, norms, data, hm, ax=None, backend='numpy', **kwargs):
    """ TODO TODO """
    if backend == 'numpy':
        return plot_fit_proj_numpy(table, norms, data, hm, ax=ax, **kwargs)
    elif backend == 'root':
        return plot_fit_proj_root(table, norms, data, hm, ax=ax, **kwargs)
    else:
        raise ValueError(f'Unknown backend: {backend}')
    

from ROOT import TCanvas
from qpym2.utils import find_mode

def create_comp_plots_uvroot(comp, trace, data, hm,
                            name='test',
                            var_names = ['co60', 'u238', 'th232', 'others', 'ndbd']
                            ):
    """
    TODO: figure out what to do with `var_names` 
    """
    # debug(5, list(traces[2].posterior.data_vars), var_names)
    norms = find_mode(trace)

    CANW, CANH = 1400, 1800
    cc = TCanvas(f'can_{name}', f'can_{name}', CANW, CANH)

    cc.Divide(1, 2)
    cu = cc.GetPad(1)

    uconfig = { 'log_yaxis_hist': True,
            'out_fpath': "",
            'plot_title': f"Data vs Fit projection ({name})",
            'xaxis_title': "ESum/#sqrt{2} [keV]",
            'yaxis_title': "Counts/keV",
            'h1_title': "Data",
            'h2_title': "Fit",
            'fit': None,
            # 'hist_stats': 'nie'
            # 'out_fpath': 'test.root',
    }

    plotteru = plot_fit_proj(comp, norms, data, hm, cu, 'root',
                                            cfg=uconfig, axis=1,
                                            hname=f'{name}_u')

    # vconfig = uconfig.copy()
    # vconfig['xaxis_title'] = "#DeltaE/#sqrt{2} [keV]"
    cv = cc.GetPad(2)
    vconfig = { 'log_yaxis_hist': True,
            'out_fpath': "",
            'plot_title': f"Data vs Fit projection ({name})",
            'xaxis_title': "#DeltaE/#sqrt{2} [keV]",
            'yaxis_title': "Counts/keV",
            'h1_title': "Data",
            'h2_title': "Fit",
            'fit': None,
            # 'hist_stats': 'nie'
            # 'out_fpath': 'test.root',
    }
    plotterv = plot_fit_proj(comp, norms, data, hm, cv, 'root',
                                            cfg=vconfig, axis=0,
                                            hname=f'{name}_v', rebin=10)

    return plotteru, plotterv, cc