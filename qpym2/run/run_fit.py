from qpym2.io import hdf
from qpym2 import fit

# TODO: use default names from hists.py instead of hardcoding.

def read_input(input_fname, dataname='data', signal_name='ndbd'):
    print(input_fname, dataname)
    mctable = hdf.read_bkg_model(input_fname)
    data = hdf.read_data(input_fname, dataname=dataname)
    signal_df = hdf.read_signal(input_fname, signal_name)

    return mctable, signal_df, data

def run_fit(cfg):
    """ Run the fit. """
    input_fname = f'{cfg.fit_pars.fit_dir}/input_{cfg.fit_pars.name}.h5'
    mctable, signal_df, data = read_input(input_fname, dataname='data', signal_name='ndbd')
    m2model, fit_comps = fit.setup_fit(mctable, data, signal_df, cfg.hm, smooth=cfg.hm.smooth) #, eff_blinding=175.0/2574857)

    trace = fit.run_fit(m2model, single_core=False)
    # m2model, trace, fit_comps, _, _ = fit.setup_fit(cfg.fit_pars.name, mctable, data, signal_df,
    #                                                 cfg.hm, smooth=cfg.hm.smooth,
    #                                                 run_parallel=True, unblinded=True)
    
    return m2model, trace, fit_comps
