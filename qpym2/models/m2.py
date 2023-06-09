""" Model related functions. """
import numpy as np
import pandas as pd
import pymc as pm

from qpym2.utils import log, debug
from qpym2 import hists as hist_tools

_WIDTH_NSIGMA = 5
def priors(comps, pars='normal', ntotal=0, nsigma=_WIDTH_NSIGMA):
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
    sigma = nsigma * np.sqrt(ntotal)

    outpars = []
    if isinstance(pars, str): pars = [pars] * len(comps)
    if isinstance(pars, list): 
        assert len(pars) == len(comps)
        for val, p in zip(comps, pars):
            if isinstance(p, str):
                if p == 'normal': outpars.append((p, {'mu': val, 'sigma': sigma}))
                elif p == 'halfnormal': outpars.append((p, {'mu': val, 'sigma': sigma}))
                elif p == 'flat':
                    min = np.max([0, val-sigma])
                    max = np.min([ntotal, val+sigma])
                    outpars.append((p, {'lower': min, 'upper': max}))
                else:
                    raise ValueError(f'Unsupported prior type {p}')
                
    debug(outpars)
    return outpars

comp_cfg = [{ 'name': 'co60', 'filter_regex': '60Co'},
            { 'name': 'u238', 'filter_regex': '238U|210Pb'},
            { 'name': 'th232', 'filter_regex': '232Th|208Pb'},
            { 'name': 'others', 'filter_regex_inverse': '60Co|238U|210Pb|232Th|208Pb'}]

def make_bkgcomps_table(comp_cfgs, mctable, smooth=None):
    """ Create the fit input components table, including the merged histograms """
    # TODO: DEBUG: we probably don't want to add a column to the mctable here...
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
    if smooth:
        hists = [hist_tools.smooth_nph2(h2np) for h2np in hists ]
        
    integrals = [ np.sum(hist) for hist in hists ]
    
    data = {
        'name': [comp['name'] for comp in comp_cfgs],
        'mchist': hists,
        'integral': integrals,
        'prior': None,
    }
        
    table = pd.DataFrame(data=data)

    table['pdf_hist'] = table['mchist']/table['integral']
    table.set_index('name', inplace=True)
    return table[['prior', 'integral', 'pdf_hist']]

def get_comps(mctable, signal_df, comp_cfg, smooth=None, eff_blinding=0, ntotal=0):
    """ Create the components as a pandas.DataFrame by merging the mctable and signal.
     
    TODO:
    """
    fit_comps = make_bkgcomps_table(comp_cfg, mctable, smooth=smooth)
    sig_hist = signal_df['mchist'][0]
    if smooth:
        sig_hist = hist_tools.smooth_nph2(sig_hist, smooth=smooth)

    signal_df['pdf_hist'] = [ sig_hist/signal_df['integral'][0] ]
    # TODO: what is this?
    # signal_df[ 'integral' ] = DATA_0v_INJECTED*signal_df[ 'integral' ]/MC_0v_NEV
    # NOTE: :above:, I have no idea. :below: DF.append was removed in pandas 2.0
    # return fit_comps.append(signal_df[[ 'integral', 'pdf_hist']])
    comps = pd.concat([fit_comps, signal_df], join='inner')  
    prior_vals = comps['integral'].copy()
    # expected mean value for the signal must be set to 0 for unblinded data,
    # or to the expected number of events for blinded data
    prior_vals['ndbd'] = eff_blinding*signal_df['integral']['ndbd']
    priors = m2.priors(prior_vals, ['normal']*4 + ['flat'], ntotal=ntotal, nsigma=5)

    comps['prior'] = priors
    comps['varname'] = comps.index
     
    return comps

def make_fik(comps, data, xedges, yedges):
    """ Create the extended likelihood matrix for each event i, in each component k. 
    
    TODO: Add proper documentation.
    """
    data_binx, data_biny = np.digitize(data.T[0], bins=xedges) , np.digitize(data.T[1], bins=yedges)
    fik = [[comp_hist[bx-1][by-1] for bx,by in zip(data_binx, data_biny)] for comp_hist in comps['pdf_hist']]
    return np.array(fik)

def m2_model(comps, fik, eff_blinding=0, smooth=None):
    """ Create and return a pymc3 model for the M2 fit.

    Parameters
    ----------
    comps : pandas.DataFrame
        The components table.
    fik : np.ndarray (n_events, n_components)
        The extended likelihood matrix for each event i, in each component k.
    eff_blinding : float
        The blinding efficiency (n_injected/n_mc) to use for the signal component.
    smooth : (int, )
        The smoothing factor to use for the components.

    Returns
    -------
    model : pymc.Model

    Notes
    -----
        half_normal is a pymc.Bound(pm.Normal, lower=0.0)
    """

    import pytensor.tensor as tt
    def _logp(fik, norms):
        """ Extended log likelihood function for fitting the components to the data. """
        lambda_ = sum(norms)
        ext_prob = tt.dot(norms, fik)
        ext_logp = tt.sum(tt.log( ext_prob + 1e-10) )
        return ext_logp - lambda_
    
    def _random(norms):
        """ TODO: Not implemented """
        raise NotImplementedError()

    with pm.Model() as m2:

        def _get_param(prior, name):
            if prior[0] == 'halfnormal':
                _BoundedNormal = pm.Bound(pm.Normal, lower=0.0)
                return _BoundedNormal(name, mu=prior[1]['mu'], sigma=prior[1]['sigma'])
            elif prior[0] == 'flat':
                return pm.Uniform(name, lower=prior[1]['lower'], upper=prior[1]['upper'])
            elif prior[0] == 'normal':
                return pm.Normal(name, mu=prior[1]['mu'], sigma=prior[1]['sigma'])
            else:
                raise ValueError(f'Unsupported prior type {prior[1]["type"]}')
        
        ck = [_get_param(comp['prior'], name) for name, comp in comps.iterrows()]
        # This calls the likelihood function for the observable. The Observed variable here (fik)
        # is the not a list of events, but the list of probabilities (per component) for events.
        ev_probs = pm.CustomDist('pik', ck, logp=_logp, random=_random, observed=fik)

        return m2

    
    
