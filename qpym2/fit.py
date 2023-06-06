import math
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
from qpym2 import hists as hist_tools
from qpym2.utils import debug

comp_cfg = [{ 'name': 'co60', 'filter_regex': '60Co'},
            { 'name': 'u238', 'filter_regex': '238U|210Pb'},
            { 'name': 'th232', 'filter_regex': '232Th|208Pb'},
            { 'name': 'others', 'filter_regex_inverse': '60Co|238U|210Pb|232Th|208Pb'}]

def create_bkgcomps_table(comp_cfgs, mctable, smooth=None):
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

def setup_fit(model_name, mctable, data,  signal_df, hm, smooth=None,
              comp_cfg=comp_cfg, run_parallel=True, unblinded=True):    # Create configuration for components, 
    debug('comp_cfg', comp_cfg)
    fit_comps = create_bkgcomps_table(comp_cfg, mctable, smooth=smooth)
    debug('comp_cfg1', fit_comps)
    sig_hist = signal_df['mchist'][0]
    if smooth:
        sig_hist = hist_tools.smooth_nph2(sig_hist, smooth=smooth)

    signal_df['pdf_hist'] = [ sig_hist/signal_df['integral'][0] ]
    if not unblinded:
        efficiency_blinding = 175.0/2574857.
    else:
        efficiency_blinding = 0
    
    # signal_df[ 'integral' ] = DATA_0v_INJECTED*signal_df[ 'integral' ]/MC_0v_NEV
    # NOTE: :above:, I have no idea. :below: DF.append was removed in pandas 2.0
    # fit_comps.loc[len(fit_comps)] = signal_df[[ 'integral', 'pdf_hist']]
    fit_comps = pd.concat([fit_comps, signal_df], join='inner')


    total_events = len(data)
    nsigma = 5
    prior_sig = nsigma * math.sqrt(total_events)

    priors = fit_comps.apply(lambda comp: {'type': 'normal', 'mu': comp.integral, 'sigma': prior_sig}, axis=1)
    # priors[-1] = dict({'type': 'halfnormal', 'mu': efficiency_blinding*signal_df['integral'][0], 'sigma': math.sqrt(total_events)})
    
    n_0v_expected = efficiency_blinding*signal_df['integral'][0]
    lower_0vprior = max(n_0v_expected - prior_sig, 0)
    # if lower_0vprior < 0: lower_0vprior = 0
    upper_0vprior = min(n_0v_expected + prior_sig, total_events)
    # if upper_0vprior > total_events:
    #     upper_0vprior = total_events

    priors[-1] = dict({'type': 'flat', 'lower': lower_0vprior, 'upper': upper_0vprior})

    fit_comps['prior'] = priors

    _, xedges, yedges = hist_tools.get_empty_hist(hm, return_numpy=True)
    data_binx, data_biny = np.digitize(data.T[0], bins=xedges) , np.digitize(data.T[1], bins=yedges)

    #TODO: add boundary checking and see why there are over/underflow bins
    bpt = np.array([[comp[bx-1][by-1] for bx,by in zip(data_binx, data_biny) if (bx < comp.shape[0] and by < comp.shape[1] )] for comp in fit_comps['pdf_hist']])
    # bpt = np.array([[comp[by-1][bx-1] for bx,by in zip(data_binx, data_biny) if (bx < 171 and by < 1871 )] for comp in fit_comps['pdf_hist']])

    f130 = 0.34167 # TODO: From internal note 123D (convert to nusance param?)
    NMC_0v = 1e8 # TODO:
    n_teO2 = 6.022e23*1000/159.6 # From internal note 123D
    exposure = 1038.4
    efficiency = 1 # we are not using any other efficiencies
    event_selection_efficiency = 0 
    var_0v = None

    import pytensor.tensor as tt
    def _logp_test(point, norms):
        """ Test logp function """
        l_nev = sum(norms)
        l_exp_per_ev = tt.dot(norms, point)
        l_ex = tt.sum(tt.log( l_exp_per_ev + 1e-10) )
        l_i = -1*l_nev + l_ex
        return l_i
    
    def _random_test(norms):
        """ Test random function """
        return 0
        

    with pm.Model() as m2:
        norms = [] # the RVs for norms
        comp_varnames = [] # names of the RVs for norms
        for name, comp in fit_comps.iterrows():
            prior = comp['prior']
            varname = f'c{model_name}_{comp.name}'
            if prior['type'] == 'halfnormal':
                BoundedNormal = pm.Bound(pm.Normal, lower=0.0)
                norms.append(BoundedNormal(varname, mu=comp.prior['mu'], sigma=comp.prior['sigma']))
            if prior['type'] == 'flat': #TODO elif
                norms.append(pm.Uniform(varname, lower=comp.prior['lower'], upper=comp.prior['upper']))
            elif prior['type'] == 'normal':
                norms.append(pm.Normal(varname, mu=comp.prior['mu'], sigma=comp.prior['sigma']))

            # TODO: work with 0v/signal differently?
            if comp.name == 'ndbd':
                event_selection_efficiency = comp.integral/NMC_0v
                var_0v = norms[-1]

            comp_varnames.append(varname)

        fit_comps['varname'] = comp_varnames
        evl = pm.CustomDist('evl', norms, logp=_logp_test, random=_random_test, observed=bpt)


        # debug(norms, bpt)
        # ev_likelihoods = ComponentPDF('ev_likelihoods', norms, testval=1, observed=bpt)

        # calculate the rate given the number of events expected
        # total_efficiency = event_selection_efficiency*efficiency
        # rate = pm.Deterministic('rate', var_0v/(f130*n_teO2)/(exposure*total_efficiency))

        if run_parallel:
            trace = pm.sample(5000, cores=1, chains=4, tune=1000, return_inferencedata=True, target_accept=0.9)
        else:
            trace = pm.sample(5000, chains=4, tune=1000, return_inferencedata=True, target_accept=0.9)
        az.plot_trace(trace)
        
    return m2, trace, fit_comps, xedges, yedges

