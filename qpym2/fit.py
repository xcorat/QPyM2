import math
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
from qpym2 import hists as hist_tools
from qpym2.utils import debug
from qpym2.models import m2

comp_cfg = [{ 'name': 'co60', 'filter_regex': '60Co'},
            { 'name': 'u238', 'filter_regex': '238U|210Pb'},
            { 'name': 'th232', 'filter_regex': '232Th|208Pb'},
            { 'name': 'others', 'filter_regex_inverse': '60Co|238U|210Pb|232Th|208Pb'}]

def setup_fit(mctable, data,  signal_df, hm, smooth=None,
              comp_cfg=comp_cfg,  eff_blinding=0):    # Create configuration for components,
    #efficiency_blinding = 175.0/2574857. 

    comps = m2.get_comps(mctable, signal_df, comp_cfg, smooth=smooth)

    prior_vals = comps['integral'].copy()
    # expected mean value for the signal must be set to 0 for unblinded data,
    # or to the expected number of events for blinded data
    prior_vals['ndbd'] = eff_blinding*signal_df['integral']['ndbd']
    priors = m2.priors(prior_vals, ['normal']*4 + ['flat'], ntotal=len(data), nsigma=5)

    comps['prior'] = priors

    _, xe, ye = hist_tools.get_empty_hist(hm, return_numpy=True)
    fik = m2.make_fik(comps, data, xe, ye)

    model = m2.m2_model(comps, fik, eff_blinding, smooth)

    return model, comps

def run_fit(model, single_core=False):

    with model:
        if single_core:
            trace = pm.sample(5000, cores=1, chains=4, tune=1000, return_inferencedata=True, target_accept=0.9)
        else:
            trace = pm.sample(5000, chains=4, tune=1000, return_inferencedata=True, target_accept=0.9)
        az.plot_trace(trace)
        
    return trace

