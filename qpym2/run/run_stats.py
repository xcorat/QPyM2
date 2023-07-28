from qpym2.run.create_fit_input import create_fit_input as create_input
from qpym2 import fit

def run_fit_syst(cfg, 
                 syst_pars=None,
                 chains=4, draws=5000):
    """ Run the fit. """
    # input_fname = f'{cfg.fit_pars.fit_dir}/input_{cfg.fit_pars.name}.h5'

    shift = cfg.fit_pars.shift
    mctable, signal_df, data = create_input(cfg, signal_name='ndbd', syst_test_ds=syst_pars, shift=shift, write=False)
    m2model, fit_comps = fit.setup_fit(mctable, data, signal_df, cfg.hm, smooth=cfg.hm.smooth) #, eff_blinding=175.0/2574857)
    trace = fit.run_fit(m2model, single_core=False, chains=chains, draws=draws)
    
    return m2model, trace, fit_comps

if __name__ == '__main__':
    { 'syst_ds': None, 'syst_ds_dq': 0, 'syst_dsigma': 0, 'shift': 0, }
