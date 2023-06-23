import site
import numpy as np



def run_sensitivity(fullname, hm, nsamples=1, signal_comp='ndbd'):
    """ Run sensitivity for a given fit.
    
    """
    # site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/qpym2')
    from qpym2 import posterior_predictive as pp
    from qpym2.io import mcroot
    from qpym2.io import fit as iofit
    
    _, xedges, yedges = mcroot.get_empty_hist(hm)
    # TODO: change to read_fit(fpath, fname)
    # TODO: add hist model to the fit file
    comps, trace = iofit.read_fit(fullname=fullname)

    return pp.sample_sensitivity_jobs(trace, comps, (xedges, yedges), size=nsamples, 
                                        signal_comp=signal_comp, remove_vars=['rate'])

def setup_toy_draws(trace, nsamples_per_test=1, vars=None, signal_comp='ndbd',
                        bias_min=0, bias_max=None):
    """ Setup toy posterior samples for sensitivity and bias studies.
    For bias, the signal component is set to a random value between bias_min and bias_max.
    And for sensitivity, the signal component is set to 0. Each of these is done nsamples_per_test times.

    Args:
        trace: arviz.InferenceData
            trace to sample from
        nsamples_per_test: int  
            number of samples to draw from the posterior per test.
        vars: list
            list of variables to sample
        signal_comp: str
            name of the signal component in the trace. Default is 'ndbd'
        bias_min: float
            minimum bias to test. Default is 0.
        bias_max: float
            maximum bias to test. If not given, use the 2x the value from the fit.
    """
    from numpy import random as rnd
    from qpym2 import utils
    from qpym2 import posterior_predictive as pp
    
    _n = nsamples_per_test
    _stat_func = utils.get_mode

    if vars is None:
        vars = [var for var in trace.posterior.data_vars if var != 'rate']
    _isig = vars.index(signal_comp)
    
    # Create toy data samples for sensitivity. Set the signal component to 0.
    points_sens = pp.sample_points(trace, vars=vars, stat_func=_stat_func, size=_n)
    points_sens[:, _isig] = 0
    # Create toy data samples for bias. Set the signal component to min, max given.
    if bias_max is None:
        bias_max = 2 * trace.posterior[signal_comp].mean().item() # Do I need the .item()?
    points_bias = pp.sample_points(trace, vars=vars, stat_func=_stat_func, size=_n)
    points_bias[:, _isig] = rnd.uniform(bias_min, bias_max, size=_n)

    return np.concatenate((points_sens, points_bias), axis=0)
    
def run_sample_post(comps, trace, hm, nsamples_per_test=1, signal_comp='ndbd',
                        bias_min=0, bias_max=None):
    """ Run posterior sampling for bias and sensitivity studies.
    
    Args:
        fullname: str
            full path to the fit file.
        hm: tuple
            tuple of (xedges, yedges) for the histogram.
        nsamples_per_test: int
            number of samples to draw from the posterior per test.
        signal_comp: str
            name of the signal component in the trace. Default is 'ndbd'
        bias_min: float
            minimum bias to test. Default is 0.
        bias_max: float
            maximum bias to test. If not given, use the 2x the value from the fit.

    Returns:
        m2fit_traces: list
            list of traces from the posterior sampling.
        signal_norms: list
            list of signal norms used for each posterior sampling.

    TODO: 
        - Store traces in a native format that doesn't use json-stringify.
            we need this to get rid of the crashing when # runs is too large (>500, 500->1.4GB)

    """
    import time
    from qpym2 import hists as hist_utils
    from qpym2 import utils
    from qpym2 import posterior_predictive as pp
    from qpym2.models import m2
    
    _n = nsamples_per_test

    _, xedges, yedges = hist_utils.get_empty_hist(hm)
    # TODO: change to read_fit(fpath, fname)
    # TODO: add hist model to the fit file
    # comps, trace = iofit.read_fit(fullname=fullname)
    
    vars = list(comps.varname)
    
    points = setup_toy_draws(trace, nsamples_per_test=_n, vars=vars, 
                             signal_comp=signal_comp, bias_min=bias_min, bias_max=bias_max)
    signal_norms = points[:, vars.index(signal_comp)]
    fit_comps = comps[ ['pdf_hist', 'prior', 'integral'] ]
    # Sample data from the points setup.
    toy_data = [ pp.sample_data(fit_comps, point, (xedges, yedges))[0] for point in points ]

    # fit_comps['prior'] = m2.piors(, pars=['normal']*4+['flat'], ntotal=0, nsigma=5)
    # NOTE: maybe check the time taken for this...
    # toy_data = np.array(toy_data)
    sens_pars = [ [data, comps, xedges, yedges] for data in toy_data ]
    ncores = utils.get_ncores()
    with Pool(ncores) as p:
        print('start pool...\nncores = ', len(sens_pars))
        m2fit_traces = p.map(pp.fit_m2_formp, sens_pars)

        timeout = 60*30 # 30 minutes
        start = time.time()
        while time.time()-start < timeout:
            if p._taskqueue.empty():
                break
            else:
                time.sleep(1)
        else:
            p.terminate()
            p.join()
            # raise RuntimeError('Timeout in run_sample_post')
            print('Timeout in run_sample_post')
            p.close()

    print('end pool...', len(m2fit_traces))

    return m2fit_traces, signal_norms


import pandas as pd

def run_single_call(cfg):
    from qpym2.io import fit as fit_io
    outdir = cfg.fit_pars.fit_dir
    fname = f'fit_{cfg.fit_pars.name}'
    fullname = f'{outdir}/{fname}'
    comps, trace = fit_io.read_fit(fullname=fullname)
    traces, signal_norms = run_sample_post(comps, trace, cfg.hm, nsamples_per_test=500, signal_comp='ndbd',
                            bias_min=0, bias_max=200)
    
    
    sens_out = f'{outdir}/toy_fits_{fname}.h5'
    with pd.HDFStore(sens_out) as store:
        store['m2main_sens'] = pd.DataFrame({'trace': traces, 'signal_norm': signal_norms})
    
    print(signal_norms)


if __name__ == '__main__':
    import site
    import sys
    from multiprocessing import Pool
    
    site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/')
    site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/cfg/nat21/roi310_main/')
    from qpym2.cfg.common_vars import cfg_m2nat21_may23 as cfg

    run_single_call(cfg)

""" Move below to a new function """
    # run_name, shift_name = sys.argv[1:] #'noco2_no0v1740_v150_smooth', 'ushiftp05'
    # hm = cfg.hm

    # cols = pd.MultiIndex.from_product(
    #     [
    #         [run_name],
    #         [shift_name],
    #         ['trace', 'signal_norm'] 
    #     ],
    #     names=('run_config', 'shift', 'fit_out') 
    # )
    # ndbd_stats = pd.DataFrame(columns=cols)

    # fullname = f'{outdir}/fit_{run_name}_{shift_name}'

    # traces, signal_norms = run_sample_post(fullname, hm, nsamples_per_test=100, signal_comp='ndbd',
    #                             bias_min=0, bias_max=200)
    
    # print(signal_norms)
    # # traces = run_sensitivity(fullname, hm, nsamples=1, signal_comp='ndbd')
    # # fit_func, pars = run_sensitivity(fullname, hm, nsamples=100, signal_comp='ctest_ndbd')

    # # ncores = check_node()
    # # with Pool(ncores) as p:
    # #     m2fit_traces = p.map(fit_func, pars)

    # # m2fit_traces = list(m2fit_traces)

    # ndbd_stats[run_name, shift_name, 'trace'] = traces
    # ndbd_stats[run_name, shift_name, 'signal_norm'] = signal_norms
    # sens_out = f'{outdir}/toy_fits/toy_fits_{run_name}_{shift_name}.h5'
    # ndbd_stats.to_hdf(sens_out, 'm2main_sens')
