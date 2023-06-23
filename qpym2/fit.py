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

def setup_fit(mctable, data, signal_df, hm, smooth=None,
              comp_cfg=comp_cfg,  eff_blinding=0, model_name='m2'):    # Create configuration for components,
    #efficiency_blinding = 175.0/2574857. 

    comps = m2.get_comps(mctable, signal_df, comp_cfg, smooth=smooth, eff_blinding=eff_blinding)

    _, xe, ye = hist_tools.get_empty_hist(hm, return_numpy=True)
    fik = m2.make_fik(comps, data, xe, ye)

    if model_name == 'm2':
        model = m2.m2_model(comps, fik)
    elif model_name == 'm2_binned':
        model = m2.m2_binned(comps, fik)

    return model, comps

def run_fit(model, single_core=False, **kwargs):
    args_ = dict(draws=5000, chains=4, tune=1000, return_inferencedata=True, target_accept=0.9,
                                idata_kwargs={"log_likelihood": True})
    args_.update(kwargs)
    if single_core:
        args_['cores'] = 1
        
    with model:
        trace = pm.sample(**args_)    
        az.plot_trace(trace)
        
    return trace

