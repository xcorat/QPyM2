""" Model related functions. """
import math
import numpy as np
import pandas as pd
import pymc as pm

from qpym2.utils.logging import log, debug
from qpym2.utils import hists as hist_tools

_WIDTH_NSIGMA = 5
def priors(comps, pars='normal', mcstds=None, ntotal=0, nsigma=_WIDTH_NSIGMA):
    """ Create and return a list of priors for the fit components.
    
    Parameters
    ----------
    comps : list or list-like
        List of integers representing the number of events in each component.
    pars : str | list(str) 
        if str: The type of prior to be used for all components.
        if list(str): A list of strings representing the prior_type.
            `len(pars)` must be equal to `len(comps)`. 
    ntotal : int
        The total number of events in the data. If 0, the integral (sum) of the
        components will be used.
    nsigma : int
        The number of sigma to use for the width of the prior.

    Returns
    -------
    priors : list
        The list of priors.

    Notes
    -----
        Supported prior types and populating parameters: (with `val` is the value from `comp`)
            - 'normal' : normal distribution, parameters: mu, sigma 
            - 'halfnormal' : halfnormal distribution, parameters: mu, sigma
            - 'flat' : flat distribution, parameters: lower, upper

            mu, sigma, lower, upper are calculated as follows:
                mu = val
                sigma = nsigma * sqrt(ntotal)
                lower = max(0, val-sigma)
                upper = min(ntotal, val+sigma)

    Examples
    --------
    >>> from qpym2.models import m2
    >>> priors = fit_comps['integral'].copy()
    >>> priors['ndbd'] = 0
    >>> p = m2.priors(priors, ['normal']*4 + ['flat'])
    >>> p
    [('normal', {'mu': 0, 'sigma': 5}), ('normal', {'mu': 0, 'sigma': 5}), ('normal', {'mu': 0, 'sigma': 5}), ('normal', {'mu': 0, 'sigma': 5}), ('flat', {'lower': 0, 'upper': 5})]
    
    TODO: 
        - I don't like the `get` in the name.
    """

    if ntotal == 0: ntotal = comps.sum()
    if mcstds is not None and isinstance(mcstds, list):
        assert len(mcstds) == len(comps)
        sigma = [ math.sqrt(nev + mcstd**2) for nev, mcstd in zip(comps, mcstds) ]
    else:
        sigma = [nsigma * np.sqrt(ntotal)]*len(comps)

    outpars = []
    if isinstance(pars, str): pars = [pars] * len(comps)
    if isinstance(pars, list): 
        assert len(pars) == len(comps)
        for i in range(len(comps)):
            val = comps.iloc[i] if hasattr(comps, 'iloc') else comps[i]
            p, s = pars[i], sigma[i]
            if isinstance(p, str):
                if p == 'normal': outpars.append((p, {'mu': val, 'sigma': s}))
                elif p == 'halfnormal': outpars.append((p, {'mu': val, 'sigma': s}))
                elif p == 'flat':
                    # TODO: maybe it is better to have a flat prior from 0 to ntotal
                    min = 0 #np.max([0, val-sigma])
                    max = ntotal #np.min([ntotal, val+sigma])
                    outpars.append((p, {'lower': min, 'upper': max}))
                else:
                    raise ValueError(f'Unsupported prior type {p}')
                
    debug(outpars)
    return outpars


def make_bkgcomps_table(comp_cfgs, mctable, smooth=None):
    """ Create the fit input components table, including the merged histograms """
    # TODO: DEBUG: we probably don't want to add a column to the mctable here...
    #       We can move this block below the component table creation and track 
    #       the norm separately instead?
    if not 'norm' in mctable:
        if 'mode' in mctable:
            mctable['norm'] = mctable['mode']
            print(f"using 'mode' as the norm...")
        else:
            mctable['norm'] =  mctable['mean']
            print(f"using 'mean' as the norm...")
        
    def _get_component_table(comp):
        if 'filter_regex' in comp:
            return mctable.filter(regex=comp['filter_regex'], axis=0)
        elif 'filter_regex_inverse' in comp:
            filtered = mctable.filter(regex=comp['filter_regex_inverse'], axis=0)
            return mctable.drop(filtered.index)
        
    component_tables = [ _get_component_table(comp) for comp in comp_cfgs ]
    hists = [ (ct['norm'] * ct['mchist']).sum()  for ct in component_tables ]

    # TODO: remove the negative bins
    for h in hists:
        if np.any(h<0): debug('WARNING: negative bins in the histogram, setting to 0.')
        h[h<0] = 0

    if smooth:
        hists = [hist_tools.smooth_nph2(h2np, smooth=smooth) for h2np in hists ]
        
    integrals = [ np.sum(hist) for hist in hists ]
    
    data = {
        'name': [comp['name'] for comp in comp_cfgs],
        'mchist': hists,
        'integral': integrals,
        'prior': None,
    }
        
    table = pd.DataFrame(data=data)

    debug( table['integral'] )
    table['pdf_hist'] = table['mchist']/table['integral']
    table.set_index('name', inplace=True)
    return table[['prior', 'integral', 'pdf_hist']]

def get_comps(mctable, signal_df, comp_cfg, smooth=None, eff_blinding=0, ntotal=0):
    """ Create the components as a pandas.DataFrame by merging the mctable and signal.
     
    TODO:
        - Reading the mctable and adding the signal to the comps table
             may need to be separated into two functions to make it easier
             to work with blinded/injected data.
    """
    fit_comps = make_bkgcomps_table(comp_cfg, mctable, smooth=smooth)
    sig_hist = signal_df['mchist'].iloc[0]
    if smooth:
        sig_hist = hist_tools.smooth_nph2(sig_hist, smooth=smooth)

    signal_df['pdf_hist'] = [ sig_hist/signal_df['integral'].iloc[0] ]
    # TODO: what is this?
    # signal_df[ 'integral' ] = DATA_0v_INJECTED*signal_df[ 'integral' ]/MC_0v_NEV
    # NOTE: :above:, I have no idea. :below: DF.append was removed in pandas 2.0
    # return fit_comps.append(signal_df[[ 'integral', 'pdf_hist']])
    comps = pd.concat([fit_comps, signal_df], join='inner')  
    prior_vals = comps['integral'].copy()

    print(comp_cfg)
    # expected mean value for the signal must be set to 0 for unblinded data,
    # or to the expected number of events for blinded data
    prior_vals['ndbd'] = eff_blinding*signal_df['integral']['ndbd']
    
    ncomps_bkg = len(comp_cfg)
    sigmas = [ comp['std'] for comp in comp_cfg ]
    mcstds = [*sigmas, 0]
    prs = priors(prior_vals, ['halfnormal']*ncomps_bkg + ['flat'], mcstds=mcstds)
    # prs = priors(prior_vals, ['normal']*ncomps_bkg + ['flat'], ntotal=ntotal, nsigma=5)
                

    #TODO: TODO:
    # comps['prior']['th232'] = ('normal', {'mu': 622.2579596900074, 'sigma': 97})
    # prs[3][1]['sigma'] = 10
    # prs[2] = ('normal', {'mu': 622.2579596900074, 'sigma': 27})
    comps['prior'] = prs
    # print(comps['prior']['th232'])
    print('new priors: ', comps['prior'])
    comps['varname'] = comps.index     
    return comps

def make_fik(comps, data, xedges, yedges):
    """ Create the extended likelihood matrix for each event i, in each component k. 
    
    TODO: Add proper documentation.
    """
    data_binx = np.clip(np.digitize(data.T[0], bins=xedges), 1, len(xedges) - 1)
    data_biny = np.clip(np.digitize(data.T[1], bins=yedges), 1, len(yedges) - 1)

    fik = [[comp_hist[bx-1][by-1] for bx,by in zip(data_binx, data_biny)] for comp_hist in comps['pdf_hist']]
    return np.array(fik)

def m2_model(comps, fik, rate_pars={'ndbd_name':'ndbd'}):
    """ Create and return a pymc model for the M2 fit.

    Parameters
    ----------
    comps : pandas.DataFrame
        The components table with `name` and `prior` columns. The `prior` column
        must be a list of tuples with the prior type and parameters.
    fik : np.ndarray (n_events, n_components)
        The extended likelihood matrix for each event i, in each component k.

    Returns
    -------
    model : pymc.Model

    Notes
    -----
        halfnormal prior uses pm.Truncated (PyMC5 replacement for deprecated pm.Bound)
    """

    import pytensor.tensor as tt
    def _logp(fik, *norms):
        """ Extended log likelihood function for fitting the components to the data. """
        # print(norms)
        norms_list = list(norms) if not isinstance(norms, (list, tuple)) else norms
        lambda_ = tt.sum(tt.stack(norms_list))  # PyMC5: Use tensor ops instead of Python sum()
        ext_prob = tt.dot(norms_list, fik)
        ext_logp = tt.sum(tt.log( ext_prob + 1e-10) )
        return ext_logp - lambda_
    
    def _random(*norms):
        """ TODO: Not implemented """
        raise NotImplementedError()

    with pm.Model() as m2:

        def _get_param(prior, name):
            if prior[0] == 'halfnormal':
                # PyMC5: Use Truncated instead of deprecated Bound
                return pm.Truncated(name,
                                pm.Normal.dist(mu=prior[1]['mu'], sigma=prior[1]['sigma']), lower=0.0)
            elif prior[0] == 'flat':
                return pm.Uniform(name, lower=prior[1]['lower'], upper=prior[1]['upper'])
            elif prior[0] == 'normal':
                return pm.Normal(name, mu=prior[1]['mu'], sigma=prior[1]['sigma'])
            else:
                raise ValueError(f'Unsupported prior type {prior[1]["type"]}')
            
        priors = comps['prior']
        names = comps.index
        
        ck = [_get_param(prior, name) for name, prior in zip(names, priors)]
        
        # This calls the likelihood function for the observable. The Observed variable here (fik)
        # is the not a list of events, but the list of probabilities (per component) for events.
        ev_probs = pm.CustomDist('pik', ck, logp=_logp, random=_random, observed=fik)

        # TODO: TODO
        # if rate_params:
        #     f130 = 0.34167 # TODO: From internal note 123D (convert to nnusance param?)
        #     NMC_0v = 1e8 # TODO:
        #     n_teO2 = 6.022e23*1000/159.6 # From internal note 123D
        #     exposure = 1122.9 #1038.4


        #     var_ndbd = filter(lambda param: param.name==ndbd_name, ck)

        return m2

def m2_binned(comps, fik):
    """ Create and return a pymc model (binned likelihood) for the M2 fit.

    TODO: this implementation doesn't seem to work. 
    """

    raise NotImplementedError('m2_binned is not implemented yet.')

    import pytensor.tensor as tt
    def _logp(fik, norms):
        """ Extended log likelihood function for fitting the components to the data. """
        lambda_ = tt.sum(tt.stack(norms))  # PyMC5: Use tensor ops instead of Python sum()
        ext_prob = tt.dot(norms, fik)
        ext_logp = tt.sum(tt.log( ext_prob + 1e-10) )
        return ext_logp - lambda_
    
    def _random(norms):
        """ TODO: Not implemented """
        raise NotImplementedError()

    with pm.Model() as m2:

        m2priors = comps['prior']
        names = comps.index
        hists = np.array([comps['pdf_hist'][name] for name in names])


    with pm.Model() as m2:

        def _get_param(prior, name):
            if prior[0] == 'halfnormal':
                # PyMC5: Use Truncated instead of deprecated Bound
                return pm.Truncated(name, pm.Normal.dist(mu=prior[1]['mu'], sigma=prior[1]['sigma']), lower=0.0)
            elif prior[0] == 'flat':
                return pm.Uniform(name, lower=prior[1]['lower'], upper=prior[1]['upper'])
            elif prior[0] == 'normal':
                return pm.Normal(name, mu=prior[1]['mu'], sigma=prior[1]['sigma'])
            else:
                raise ValueError(f'Unsupported prior type {prior[1]["type"]}')

        
        ck = [_get_param(prior, name) for name, prior in zip(names, m2priors)]

        # data_var = pm.Data('data', h2)
        # mu = np.array([c*hist for c, hist in zip(ck, hists)])
        mu = pm.math.dot(hists.T, ck).T
        obs = pm.Poisson('obs', mu=mu, observed=h2)
        return m2
    
def sample_rate(hist, xedges, yedges, comps, return_trace=False, sample_pca=None):
    """

    Args:
        hist (_type_): _description_
        xedges (_type_): _description_
        yedges (_type_): _description_
        comps (_type_): _description_

    Returns:
        _type_: _description_
    """

    import arviz as az

    f130 = 0.34167 # TODO: From internal note 123D (convert to nnusance param?)
    NMC_0v = 1e8 # TODO:
    n_teO2 = 6.022e23*1000/159.6 # From internal note 123D
    exposure = 1122.9 #1038.4

    if sample_pca:
        efficiency = sample_pca[0]
        eff_sigma = sample_pca[1]
    else:
        efficiency = 1
    event_selection_efficiency = 0 

    data_toymc = sample_hist(hist, xedges, yedges)
    data_toymc_binx, data_toymc_biny = np.digitize(data_toymc.T[0], bins=xedges) , np.digitize(data_toymc.T[1], bins=yedges)

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
                # PyMC5: Use Truncated instead of deprecated Bound
                norms.append(pm.Truncated(varname, pm.Normal.dist(mu=comp.prior['mu'], sigma=comp.prior['sigma']), lower=0.0))
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
        rate = pm.Deterministic('rate', var_0v/(f130*n_teO2)/(exposure*total_efficiency))

        trace = pm.sample(4000, tune=1000, chains=2, return_inferencedata=True, target_accept=0.9, progressbar=False)
        print(trace.posterior.mean())
        # az.plot_trace(trace)

    if return_trace:
        return trace
    else:
        rate = trace.posterior['rate'].values #.mean().item()
        return rate
