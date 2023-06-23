"""
Contains routines for posterior predictive tests. Uses Arviz traces
"""

import numpy as np
import pandas as pd
import pymc as pm

from qpym2.utils import debug

def sample_points(trace, vars=None, stat_func=None, size=1):
    """ Sample points from the trace 
    Args: 
        trace (arviz.InferenceData): trace to sample from
        vars: list
            list of variables to sample
        stat_func:
            estimator to use for sampling or None to use random.
            The default is None. `stat_funcs` string is passed on to 
            `arviz` `sample` function
        size: int
            number of points to sample #TODO make compatible with matrices
    Returns:
        Array
            array of sampled points. If vars is None, all variables are sampled
            in the order they appear in the trace. Otherwise the order is given
            by vars.
    """
    
    if vars is None:
        vars = list(trace.posterior.data_vars)

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

# return tools_old.sample_rate(hsum, xedges, yedges, fit_comps, return_trace=True)

from qpym2.hists import sample_hist
from qpym2.hists import  get_sum_hist

def sample_data(comps, norms, bin_edges, size=1):
    """ Sample data from the model

    Args:
        comps: pd.DataFrame
            dataframe of components with columns: name, pdf_hist
        norms: list
            list of normalization factors that match the components in comps.
        bin_edges: tuple
            tuple of bin edges on component (as xedges, yedges)
        size: int
            number of points to sample. Default is 1.

    Returns:
        numpy.array
            array of sampled data points.
    """
    hsum = get_sum_hist(comps, norms)
    return [ sample_hist(hsum, *bin_edges) for _ in range(size) ]

def fit_m2_old(data, comps, xedges, yedges, multiprocess=True):
    """ Fit data to model 
    
    Args:
        data: numpy.array
            2D array of data points to fit
        comps: pd.DataFrame
            dataframe of components with columns: index (name), pdf_hist, prior, integral.
            The integral column is used to calculate the event selection efficiency.
        xedges: numpy.array
            array of bin edges in x
        yedges: numpy.array
            array of bin edges in y
        multiprocess: bool
            whether to use multiprocessing when sampling with pymc3. Default is True.
        if turned off, the `cores` argument in `pm.sample` is set to 1. Turn this off
        when calling from a multiprocessing context. 
    """

    # TODO: Move model specific stuff to seperate files per model, include 
    #       the model constants in a global file...

    import arviz as az
    import pymc3 as pm

    from qm2mc.ComponentPDF import ComponentPDF

    f130 = 0.34167 # TODO: From internal note 123D (convert to nnusance param?)
    NMC_0v = 1e8 # TODO:
    n_teO2 = 6.022e23*1000/159.6 # From internal note 123D
    exposure = 1038.4
    efficiency = 1
    event_selection_efficiency = 0 

    data_toymc_binx, data_toymc_biny = np.digitize(data.T[0], bins=xedges) , np.digitize(data.T[1], bins=yedges)

    #TODO: add boundary checking and see why there are over/underflow bins
    model_name = 'toymc'

    # bpt_toymc = np.array([[comp[by-1][bx-1] for bx,by in zip(data_toymc_binx, data_toymc_biny) if (bx < 171 and by < 1871 )] for comp in comps['pdf_hist']])
    bpt_toymc = np.array([[comp[bx-1][by-1] for bx,by in zip(data_toymc_binx, data_toymc_biny) if (bx < comp.shape[0] and by < comp.shape[1] )] for comp in comps['pdf_hist']])

    with pm.Model() as m2:
        norms = [] # the RVs for norms
        comp_varnames = [] # names of the RVs for norms
        for name, comp in comps.iterrows():
            prior = comp['prior']
            varname = f'c{model_name}_{comp.name}'
            if prior['type'] == 'halfnormal':
                BoundedNormal = pm.Bound(pm.Normal, lower=0.0)
                norms.append(BoundedNormal(varname, mu=comp.prior['mu'], sigma=comp.prior['sigma']))
            elif prior['type'] == 'flat': #TODO elif
                norms.append(pm.Uniform(varname, lower=comp.prior['lower'], upper=comp.prior['upper']))
            elif prior['type'] == 'normal':
                norms.append(pm.Normal(varname, mu=comp.prior['mu'], sigma=comp.prior['sigma']))

            # TODO: work with 0v/signal differently?
            if comp.name == 'ndbd':
                event_selection_efficiency = comp.integral/NMC_0v
                var_0v = norms[-1]

            comp_varnames.append(varname)

        # fit_comps['varname'] = comp_varnames
        ev_likelihoods = ComponentPDF('ev_likelihoods', norms, testval=1, observed=bpt_toymc)

        # calculate the rate given the number of events expected
        total_efficiency = event_selection_efficiency*efficiency

        # TODO: do we need this? 
        rate = pm.Deterministic('rate', var_0v/(f130*n_teO2)/(exposure*total_efficiency))

        if multiprocess:
            trace = pm.sample(4000, tune=1000, chains=2, cores=1, return_inferencedata=True,
                              target_accept=0.9, progressbar=False)
        else:
            trace = pm.sample(4000, tune=1000, chains=2, return_inferencedata=True,
                              target_accept=0.9, progressbar=False)
        print(trace.posterior.mean())
        
        return trace
  
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
            whether to use multiprocessing when sampling with pymc. Default is False.
    """
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
    
    comps['prior'] = comps.apply(lambda comp: _update_prior_sigma(comp['prior'], len(data)), axis=1)

    fik = m2.make_fik(comps, data, xe, ye)
    model = m2.m2_model(comps, fik)

    with model:
        if single_core:
            trace = pm.sample(4000, tune=1000, chains=2, return_inferencedata=True,
                              target_accept=0.9, progressbar=False)
        else:
            trace = pm.sample(4000, tune=1000, chains=2, cores=1, return_inferencedata=True,
                              target_accept=0.9, progressbar=False)
        debug('means: ', trace.posterior.mean())
        
        return trace

def fit_m2_formp(pars):
    """ Fit data to model using nultiprocessing. pars are expanded and passed to fit_m2.
    
    Returns:
        forward return from fit_m2.
    """
    debug('in fit_m2_formp:', len(pars)) # TODO: remove
    data, comps, xedges, yedges = pars
    # print(data, comps, xedges, yedges)
    
    return fit_m2(data, comps, xedges, yedges)

def sample_sensitivity_jobs(trace, comps, bin_edges, size=1, signal_comp='ndbd',
                       remove_vars=[]):
    """ Sample sensitivity from trace 

    Args:
        trace: arviz.InferenceData
            trace to sample from
        comps: pd.DataFrame
            dataframe of components with columns: name, pdf_hist, prior
        bin_edges: tuple
            tuple of bin edges on component (as xedges, yedges)
        size: int
            number of points to sample. Default is 1.
        signal_comp: str
           name of the signal component in the trace. Default is 'ndbd'
        remove_vars: list
            list of variables to remove from the trace before sampling. Useful
            for removing deterministic variables like `rate` from the trace.

    Returns:
        numpy.array
            array of sampled sensitivities.

        Sensitivity sampling is done as follows:
        1. Sample the posterior from trace for the background only model,
        2. Sample (toy-)data from the sum histogram constructed using the norms
    from above (poisson smeared)
        3. Fit the data to the background only model and return the posterior.

    TODO:
        - Add functionality to skip more than one component (ex. Co-60).
        - Or better, conver this to return only the fit_function and list of
          parameters per fit call.
    """
    from qpym2 import tools_noqpymc as tools
    from multiprocessing import Pool

    _stat_func = tools.get_mode

    vars = [var for var in trace.posterior.data_vars if var not in remove_vars]
    bkg_vars = [var for var in vars if var != signal_comp]
    # sample_comps = comps.drop('ndbd') #TODO: make this remove by the comp_name ('ctest_ndbd')
    print(bkg_vars, signal_comp)

    points = sample_points(trace, vars=bkg_vars, stat_func=_stat_func, size=size)

    # set 0vbb to 0
    points = np.insert(points, comps.index.get_loc('ndbd'), 0, axis=1)
    toy_data = [ sample_data(comps, point, bin_edges)[0] for point in points ]

    sens_comps = comps[ ['pdf_hist', 'prior', 'integral'] ]
    sens_pars = [ (data, sens_comps, *bin_edges) for data in toy_data ]

    return fit_m2_formp, sens_pars
    # with Pool(1) as p:
    #     print('start pool...', len(toy_data))
    #     sens_traces = p.map(fit_m2_formp, sens_pars)
        
    #     print('end pool...', len(sens_traces))

    # return sens_traces


