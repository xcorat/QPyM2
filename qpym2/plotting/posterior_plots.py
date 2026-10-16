""" Functions needed for plots related to the posterior distributions. """ 

import numpy as np
from matplotlib import pyplot as plt

from qpym2 import utils
from qpym2.utils import hists as hist_tools
from qpym2.utils.mathy import get_1sig_around_mode, get_mode

from qpym2.utils.logging import debug

# Hydra config helper for plotting backend
from qpym2.cfg import get_config
import matplotlib

def _ensure_backend(cfg=None):
    if cfg is None:
        try:
            cfg = get_config()
        except Exception:
            cfg = None
    if cfg is not None:
        try:
            backend = cfg.system.matplotlib_backend
            if backend:
                matplotlib.use(str(backend))
        except Exception:
            pass

def plot_posterior(trace, varname='ctest_ndbd', ax=None, cfg=None, **kwargs):
    """ Plot the posterior distribution for a given variable. 

    Accepts optional Hydra `cfg` (omegaconf DictConfig) to set matplotlib backend and styles.
    """

    _ensure_backend(cfg)

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


_STAT_COLORS = {
    'mode':   ('C1', '--'),
    'hdi':    ('C2', ':'),
    'mean':   ('C3', '-.'),
    'median': ('C4', '-.'),
    'p90':    ('C5', ':'),
}


def plot_1d_posterior(samples, weights=None, ax=None, bins=60, color='steelblue',
                      lw=1.5, label=None, lines=None, stats=None, xlabel=None,
                      ylabel='density', title=None, **kwargs):
    """Plot a 1-D posterior histogram with optional stat annotations and lines.

    Parameters
    ----------
    samples : array-like
        1-D posterior samples.
    weights : array-like, optional
        Per-sample weights (same length as ``samples``).
    ax : matplotlib.axes.Axes, optional
    bins : int or array, optional
        Histogram bins (default 60).
    color, lw, label : plot styling.
    lines : list of dict, optional
        Custom vertical lines.  Each dict: ``x`` (required), ``label``,
        ``color``, ``ls``, ``lw``.
    stats : list of str, optional
        Auto-computed statistic overlays: ``'mode'``, ``'hdi'``, ``'mean'``,
        ``'median'``, ``'p90'``.
    xlabel, ylabel, title : str, optional
    **kwargs : passed to ``ax.step``.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    samples = np.asarray(samples).flatten()
    if weights is not None:
        weights = np.asarray(weights).flatten()

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))

    h, xe = np.histogram(samples, bins=bins, weights=weights, density=True)
    xc = 0.5 * (xe[:-1] + xe[1:])
    ax.step(xc, h, where='mid', color=color, lw=lw, label=label, **kwargs)
    ymax = h.max()

    if stats:
        _mode, _hdi_low, _hdi_high = None, None, None
        for stat in stats:
            stat_color, stat_ls = _STAT_COLORS.get(stat, ('grey', '--'))
            if stat == 'mode':
                _mode = get_mode(samples) if _mode is None else _mode
                ax.axvline(_mode, color=stat_color, lw=lw, ls=stat_ls,
                           label=f'mode={_mode:.2f}')
            elif stat == 'hdi':
                _mode, _hdi_low, _hdi_high = get_1sig_around_mode(samples)
                ax.axvline(_hdi_low,  color=stat_color, lw=1.0, ls=stat_ls,
                           label=f'HDI=[{_hdi_low:.2f}, {_hdi_high:.2f}]')
                ax.axvline(_hdi_high, color=stat_color, lw=1.0, ls=stat_ls)
                ax.fill_betweenx([0, ymax], _hdi_low, _hdi_high,
                                 alpha=0.12, color=stat_color)
            elif stat == 'mean':
                _mean = np.average(samples, weights=weights)
                ax.axvline(_mean, color=stat_color, lw=lw, ls=stat_ls,
                           label=f'mean={_mean:.2f}')
            elif stat == 'median':
                _med = np.median(samples)
                ax.axvline(_med, color=stat_color, lw=lw, ls=stat_ls,
                           label=f'median={_med:.2f}')
            elif stat == 'p90':
                _p90 = np.percentile(samples, 90)
                ax.axvline(_p90, color=stat_color, lw=lw, ls=stat_ls,
                           label=f'p90={_p90:.2f}')

    if lines:
        for line in lines:
            ax.axvline(line['x'], color=line.get('color', 'grey'),
                       lw=line.get('lw', 1.0), ls=line.get('ls', '--'),
                       label=None) #line.get('label'))

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    # ax.legend(fontsize=8)
    return ax


def plot_fit_proj_numpy(table, norms, data, hm, ax=None, cfg=None, esum=False, **kwargs):
    """Plot the fit projection (numpy backend).

    Args:
        table: fit components DataFrame with 'pdf_hist' column
        norms: array of component normalizations (e.g. from find_mode)
        data: (N, 2) array of event coordinates in u, v space
        hm: TempConfig histogram model
        ax: matplotlib axes (created if None)
        cfg: optional Hydra config for backend/styles
        esum: if True, scale x-axis from u-space to ESum = u*√2 [keV]
    """
    _ensure_backend(cfg)

    nbins, range_ = hist_tools.get_hist_settings(hm)
    h2_data, xedges, yedges = np.histogram2d(data.T[0], data.T[1], bins=nbins, range=range_)
    h2_fit = hist_tools.get_sum_hist(table, norms)

    hu_data = hist_tools.get_h2proj(h2_data, axis=1)
    hu_fit = hist_tools.get_h2proj(h2_fit, axis=1)
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    # Scale x-axis to ESum if requested
    plot_x = xedges[:-1] * np.sqrt(2) if esum else xedges[:-1]
    xlabel = 'ESum [keV]' if esum else 'u = ESum/√2 [keV]'

    ax.errorbar(plot_x, hu_data, yerr=np.sqrt(hu_data), fmt='.', label='data')
    ax.plot(plot_x, hu_fit, label='fit')
    ax.set_xlabel(xlabel)
    ax.set_yscale('log')
    ax.legend()

    return ax

def plot_fit_proj_npnew(table, norms, data, hm, ax=None, cfg=None, esum=False, **kwargs):
    """Plot the fit projection with ROOT/HEP histogram style (numpy backend).

    Draws data as ROOT-style white histogram bars with black Poisson error caps,
    and the fit as a line with circle markers. Assumes ROOT/HEP-style rcParams
    (font, tick direction, legend frameon, etc.) are already set globally at the
    call-site (e.g. via ``mpl.rcParams.update(...)`` in the notebook).

    Args:
        table: fit components DataFrame with 'pdf_hist' column
        norms: array of component normalizations (e.g. from find_mode)
        data: (N, 2) array of event coordinates in u, v space
        hm: TempConfig histogram model
        ax: matplotlib axes (created if None)
        cfg: optional Hydra config for backend/styles
        esum: if True, scale x-axis from u-space to ESum = u*√2 [keV]
    """
    _ensure_backend(cfg)

    nbins, range_ = hist_tools.get_hist_settings(hm)
    h2_data, xedges, yedges = np.histogram2d(data.T[0], data.T[1], bins=nbins, range=range_)
    h2_fit = hist_tools.get_sum_hist(table, norms)

    hu_data = hist_tools.get_h2proj(h2_data, axis=1)
    hu_fit = hist_tools.get_h2proj(h2_fit, axis=1)

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    plot_x = xedges[:-1] * np.sqrt(2) if esum else xedges[:-1]
    xlabel = r'E$_\mathrm{sum}$ [keV]' if esum else r'u [keV/$\sqrt{2}$]'
    bin_w = float(np.median(np.diff(plot_x))) if len(plot_x) > 1 else 1.0

    # Data: histogram bars with steel-blue fill at 60% opacity + matching error caps
    bar_color = '#4e8ab5'
    yerr = np.maximum(np.sqrt(hu_data), 1)
    ax.bar(plot_x, hu_data, width=bin_w, align='edge',
           color=bar_color, linewidth=0,
           alpha=0.6, label='Data', zorder=2)
    ax.errorbar(plot_x + bin_w / 2, hu_data, yerr=yerr,
                fmt='none', color=bar_color, linewidth=1.0, capsize=2, zorder=3)

    # Fit: line with circle markers
    ax.plot(plot_x + bin_w / 2, hu_fit, 'o-', label='Fit',
            color='C3', linewidth=1.8, markersize=3.5)

    ax.set_xlabel(xlabel)
    ax.set_ylabel('Counts')
    ax.set_yscale('log')
    ax.minorticks_on()
    ax.legend()

    return ax

def plot_fit_proj_root(table, norms, data, hm, ax=None, projection=None, **kwargs):
    """ Plot the fit projection using ROOT. "
    
    TODO: does axis=1 X axis or Y axis?
    """
    from ROOT import TCanvas, TH2D
    from qpym2.plotting.legacy import H1ComparisonPlotter

    can = ax
    axis = kwargs.get('axis', 1)

    h2_fit = hist_tools.get_sum_hist(table, norms)
    nbins, range_ = hist_tools.get_hist_settings(hm)

    assert len(nbins) == 2, "Only 2D hists are supported for now."
    assert h2_fit.shape[0] == nbins[0] and h2_fit.shape[1] == nbins[1], "The hist model has different shape than the saved data"

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
        if projection is None:
            th1d_fit = th2d_fit.ProjectionX()
            th1d_data = th2d_data.ProjectionX()
        else:
            pname, pstart, pend = projection
            th1d_fit = th2d_fit.ProjectionX(pname, pstart, pend)
            print(th1d_fit.Integral(), th1d_data.Integral())
            th1d_data = th2d_data.ProjectionX(pname, pstart, pend)
            
    elif axis == 0:
        if projection is None:
            th1d_fit = th2d_fit.ProjectionY()
            th1d_data = th2d_data.ProjectionY()
        else:
            pname, pstart, pend = projection
            th1d_fit = th2d_fit.ProjectionX(pname, pstart, pend)
            th1d_data = th2d_data.ProjectionX(pname, pstart, pend)

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

    # plotter.draw_comparison_plot(can)

    return plotter


def plot_fit_proj_root_comps(table, norms, data, hm, ax=None, **kwargs):
    """ Plot the fit projection using ROOT. "
    
    TODO: does axis=1 X axis or Y axis?
    """
    from ROOT import TCanvas, TH2D
    from qpym2.plotting.legacy import H1ComparisonPlotter

    can = ax
    axis = kwargs.get('axis', 1)

    # h2_fit = hist_tools.get_sum_hist(table, norms)
    h2s = table['pdf_hist']*norms
    nbins, range_ = hist_tools.get_hist_settings(hm)

    assert len(nbins) == 2, "Only 2D hists are supported for now."
    assert h2s[0].shape[0] == nbins[0] and h2s[0].shape[1] == nbins[1], "The hist model has different shape than the saved data"

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

    th2d_fits = []
    data_done = False
    for h2 in h2s:
        th2d_fits.append( hist_tools.get_empty_hist(hm, return_numpy=False))
        for i in range(h2.shape[0]):
            for j in range(h2.shape[1]):
                h2.SetBinContent(i+1, j+1, h2[i, j])
                if not data_done:
                    th2d_data.SetBinContent(i+1, j+1, h2_data[i, j])
                    data_done = True

    hm.name = hname + '_data'
    th2d_data = hist_tools.get_empty_hist(hm, return_numpy=False)
    hm.name = hname


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

    # plotter.draw_comparison_plot(can)

    return plotter

def plot_fit_proj(table, norms, data, hm, ax=None, backend='numpy', esum=False, **kwargs):
    """Plot fit projection (data vs fit) with selectable backend.

    Args:
        esum: if True and backend='numpy', scale x-axis to ESum = u*√2 [keV]
    """
    if backend == 'numpy':
        return plot_fit_proj_numpy(table, norms, data, hm, ax=ax, esum=esum, **kwargs)
    elif backend == 'root':
        return plot_fit_proj_root(table, norms, data, hm, ax=ax, **kwargs)
    else:
        raise ValueError(f'Unknown backend: {backend}')
    

from ROOT import TCanvas
from qpym2.utils.mathy import find_mode

def create_comp_plots_uvroot(comp, trace, data, hm,
                            name='test',
                            var_names = ['co60', 'u238', 'th232', 'others', 'ndbd']
                            ):
    """
    Create component plots using ROOT framework.
    
    CRITICAL: Must use var_names parameter to ensure component ordering matches comps.index.
    Without explicit varnames, find_mode returns values in alphabetical order (xarray default),
    causing component value swapping if component names aren't alphabetically ordered.
    """
    # debug(5, list(traces[2].posterior.data_vars), var_names)
    norms = find_mode(trace, varnames=var_names)

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

def get_h2root(table, norms, data, hm, ax=None, **kwargs):
    """ Plot the fit projection using ROOT. "
    
    TODO: does axis=1 X axis or Y axis?
    """
    from ROOT import TCanvas, TH2D
    from qpym2.plotting.legacy import H1ComparisonPlotter

    can = ax
    axis = kwargs.get('axis', 1)

    h2_fit = hist_tools.get_sum_hist(table, norms)
    nbins, range_ = hist_tools.get_hist_settings(hm)

    assert len(nbins) == 2, "Only 2D hists are supported for now."
    assert h2_fit.shape[0] == nbins[0] and h2_fit.shape[1] == nbins[1], "The hist model has different shape than the saved data"

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

    return th2d_fit, th2d_data
