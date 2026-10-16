"""
Contains routines for posterior checks and diagnostics.
"""

import numpy as np
import pymc as pm
from numpy import random as rnd


from qpym2.utils.logging import debug
from qpym2.utils.hists import sample_hist
from qpym2.utils.hists import  get_sum_hist

def sample_points(trace, vars=None, stat_func=None, size=1):
    """ Sample points from the trace 
    Args: 
        trace (arviz.InferenceData): trace to sample from
        vars: list
            list of variables to sample, in the desired order.
            If None, defaults to trace.posterior.data_vars (ALPHABETICAL order from xarray).
        stat_func:
            estimator to use for sampling or None to use random.
            The default is None. `stat_funcs` string is passed on to 
            `arviz` `sample` function
        size: int
            number of points to sample #TODO make compatible with matrices

    Returns:
        Array
            array of sampled points. If vars is None, all variables are sampled
            in ALPHABETICAL ORDER (xarray default). Otherwise the order is given
            by vars.
            
    Warning:
        CRITICAL: xarray stores variables in ALPHABETICAL order by default.
        If vars is omitted, returned values will be alphabetically ordered,
        which may NOT match your component ordering (e.g., comps.index).
        
        Always pass explicit vars=list(comps.index) to avoid component value swapping:
            ✅ CORRECT: sample_points(trace, vars=list(comps.index), ...)
            ❌ WRONG:   sample_points(trace, ...)  # Returns alphabetical order!
    """
    # TODO: Move the vars ordering logic outside of this function.
    
    if vars is None:
        vars = list(trace.posterior.data_vars)
        #TODO: diagnostics. Remove `rate` if it's counted
        if 'rate' in vars:
            raise ValueError('rate is not a valid variable to sample from the trace.')

    stacked = trace.posterior[vars].stack(draws=("chain", "draw"))
    if stat_func is None:
        # TODO: there must be a direct method from arviz to do this...
        draws = np.array([ np.random.choice(stacked[var], size=size) for var in stacked]).T
    else:
        # point = stacked.apply(lambda x:stat_func(x)) #TODO: check if this can work
        # draws = np.repeat(point, size, axis=0) #TODO: check if this can work

        point = [ stat_func(stacked[var]) for var in stacked ]
        draws = np.array([point]*size)
        
    return draws

def sample_data(comps, norms, bin_edges, size=1):
    """ Sample data from the model components according to the norm.

    Args:
        comps: pd.DataFrame
            dataframe of components with columns: name, pdf_hist
        norms: list
            list of normalization factors that match the components in comps.
        bin_edges: tuple
            tuple of bin edges on component (as xedges, yedges)
        size: int
            number of datasets to sample. Default is 1.

    Returns:
        numpy.array
            array of sampled data points.
    """
    hsum = get_sum_hist(comps, norms)
    return [ sample_hist(hsum, *bin_edges) for _ in range(size) ]

from qpym2.utils.hists import get_empty_hist
from multiprocessing import Pool
from qpym2.utils import hw, mathy
from qpym2.utils.pool_monitor import run_pool_map_with_monitor

def toy_fits(comps, trace, model_cfg, nsamples=1):
    """ Sample points from the trace and fit the model to the sampled data.
    
    Args:
        comps: pd.DataFrame
            dataframe of components with columns: name, pdf_hist
        trace: arviz.InferenceData
            trace to sample from
        model_cfg: dict
            dictionary of model configuration
        nsamples: int
            number of samples to draw from the trace. Default is 1.


    Returns:
        list
            list of traces from the toy fits from the posterior.
    """
        
    _n = nsamples
    _, xedges, yedges = get_empty_hist(model_cfg['h2uv_model'])

    # Sample points?
    stat_func = mathy.get_mode

    # TODO: If we are sampling the mode we don't need to sample the same point?
    # CRITICAL: Must pass vars=list(comps.index) to ensure norms align with comps PDF order
    points = sample_points(trace, vars=list(comps.index), stat_func=stat_func, size=_n)
    
    print(f"\nDEBUG toy_fits: sampled points (in comps.index order {list(comps.index)}):")
    print(f"  points[0] = {points[0]}")
    
    # sample data returns an array. We are asking for one dataset, so take the first element
    data = [ sample_data(comps, p, (xedges, yedges))[0] for p in points ]

    toy_pars = [ [data, comps, xedges, yedges] for data in data ]

    ncores = hw.get_ncores()
    with Pool(ncores) as p:
        print('start pool...\nncores = ', ncores, ' ntasks = ', len(toy_pars))
        m2fit_traces = run_pool_map_with_monitor(p, 'toy_fits', fit_m2_formp, toy_pars)

    print('end pool...', len(m2fit_traces))

    return m2fit_traces

def sens_fits(comps, trace, model_cfg, nsamples=1, signal_comp='ndbd'):
    """ Sample points from the trace and fit the model to the sampled data.
    
    Args:
        comps: pd.DataFrame
            dataframe of components with columns: name, pdf_hist
        trace: arviz.InferenceData
            trace to sample from
        model_cfg: dict
            dictionary of model configuration
        nsamples: int
            number of samples to draw from the trace. Default is 1.
        signal_comp: str
            name of the signal component in the trace. Default is 'ndbd'
        stat_func: str
            estimator to use for sampling or None to use random.
            TODO: explain more?


    Returns:
        list
            list of traces from the toy fits from the posterior.
    """
        
    _n = nsamples
    _, xedges, yedges = get_empty_hist(model_cfg['h2uv_model'])

    # Sample points?
    stat_func = mathy.get_mode
    # CRITICAL: Must pass vars=list(comps.index) to ensure norms align with comps PDF order
    points = sample_points(trace, vars=list(comps.index), stat_func=stat_func, size=_n)

    vars = list(comps.index)
    _isig = vars.index(signal_comp)
    points[:, _isig] = 0

    # sample data returns an array. We are asking for one dataset, so take the first element
    data = [ sample_data(comps, p, (xedges, yedges))[0] for p in points ]

    toy_pars = [ [data, comps, xedges, yedges] for data in data ]

    ncores = hw.get_ncores()
    with Pool(ncores) as p:
        print('start pool...\nncores = ', ncores, ' ntasks = ', len(toy_pars))
        m2fit_traces = run_pool_map_with_monitor(p, 'sens_fits', fit_m2_formp, toy_pars)

    print('end pool...', len(m2fit_traces))

    return m2fit_traces

def bias_fits(comps, trace, model_cfg, nsamples=1,
            signal_comp='ndbd', bias_min=0, bias_max=None, ncores=None):
    """ Sample points from the trace and fit the model to the sampled data.
    
    Args:
        comps: pd.DataFrame
            dataframe of components with columns: name, pdf_hist
        trace: arviz.InferenceData
            trace to sample from
        model_cfg: dict
            dictionary of model configuration
        nsamples: int
            number of samples to draw from the trace. Default is 1.
        signal_comp: str
            name of the signal component in the trace. Default is 'ndbd'
        bias_min: float
            minimum bias to test. Default is 0.
        bias_max: float
            maximum bias to test. If not given, use the 2x the value from the fit.


    Returns:
        list
            list of traces from the toy fits from the posterior.
    """
        
    _n = nsamples
    _, xedges, yedges = get_empty_hist(model_cfg['h2uv_model'])

    # Sample points?
    # NOTE: We sample a poisson number of events from the final histogram, 
    #       so we don't need to sample the posterior.    
    stat_func = mathy.get_mode
    # CRITICAL: Must pass vars=list(comps.index) to ensure norms align with comps PDF order
    points = sample_points(trace, vars=list(comps.index), stat_func=stat_func, size=_n)


    vars = list(comps.index)
    _isig = vars.index(signal_comp)
    points[:, _isig] = rnd.uniform(bias_min, bias_max, size=_n)

    # sample data returns an array. We are asking for one dataset, so take the first element
    data = [ sample_data(comps, p, (xedges, yedges))[0] for p in points ]
    signal_norms = points[:, vars.index(signal_comp)]
    toy_pars = [ [data, comps, xedges, yedges] for data in data ]

    if ncores is None: ncores = hw.get_ncores()
    with Pool(ncores) as p:
        print('start pool...\nncores = ', ncores, ' ntasks = ', len(toy_pars))
        m2fit_traces = run_pool_map_with_monitor(p, 'bias_fits', fit_m2_formp, toy_pars)

    print('end pool...', len(m2fit_traces))

    return m2fit_traces, signal_norms

def fit_m2(data, comps, xe, ye, single_core=False):
    """ Fit data to model 
    
    Args:
        data: numpy.array
            2D array of data points to fit
        comps: pd.DataFrame
            dataframe of components with columns: index (name), pdf_hist, prior, integral.
            The integral column is used to calculate the event selection efficiency.
        xe: numpy.array
            array of bin edges in x
        ye: numpy.array
            array of bin edges in y
        single_core: bool
            Whether to use multiprocessing when sampling with pymc. False indicates
            support for multiprocessing. The default is False.
    """    
    print(f"DEBUG fit_m2: comps.index = {list(comps.index)}")
    from qpym2.models import m2

    # TODO: this parameter is a setting of the fit, so must be stored and 
    #       passed to this fit function.
    _WIDTH_NSIGMA = m2._WIDTH_NSIGMA

    # We first need to update the priors to match the number of events in the data
    def _update_prior_sigma(prior, ntotal, nsigma=_WIDTH_NSIGMA):
        # tuples cannot be modified, so we need to make a copy
        sigma = nsigma*np.sqrt(ntotal)
        if prior[0] == 'halfnormal' or prior[0] == 'normal':
            _prior = (prior[0], { 'mu': prior[1]['mu'], 'sigma': sigma })
        elif prior[0] == 'flat':
            _prior = ('flat', {'lower': 0, 'upper': ntotal})
        else:
            raise ValueError(f'Unknown prior type: {prior[0]}')
        return _prior
    
    # TODO: do we need to update the priors? 
    # comps['prior'] = comps.apply(lambda comp: _update_prior_sigma(comp['prior'], len(data)), axis=1)

    fik = m2.make_fik(comps, data, xe, ye)
    model = m2.m2_model(comps, fik)

    with model:
        if single_core:
            trace = pm.sample(4000, tune=1000, chains=4, return_inferencedata=True,
                              target_accept=0.9, progressbar=False,
                              idata_kwargs={"log_likelihood": True})
        else:
            trace = pm.sample(4000, tune=1000, chains=4, cores=1,
                               return_inferencedata=True,
                              target_accept=0.9, progressbar=False,
                              idata_kwargs={"log_likelihood": True})
        debug('means: ', trace.posterior.mean().items())
        
        return trace


def fit_m2_formp(pars):
    """ Fit data to model using multiprocessing. pars are expanded and passed to fit_m2.
    
    Returns:
        forward return from fit_m2.
    """
    data, comps, xedges, yedges = pars
    
    return fit_m2(data, comps, xedges, yedges)
