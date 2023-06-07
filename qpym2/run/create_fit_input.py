""" Create fit input related files.

TODO: we shouldn't need to add the qpym2 path here if the module is
installed properly. But for now, we need to add it to the path.
"""
import os
import site
import numpy as np
import pandas as pd


site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/')
from qpym2.io import bkg_model, staging, hdf
from qpym2.utils import debug

def get_uv_transform(shift, fname):
    """ Get the uv transform for the shift type and value.

    Args:
        shift (qpym2.cfg.fit_pars.Shift): shift
        fname (str): file name

    Returns:
        tuple: (u, v) transform
    """
    tr_uv = ( 'u0', 'v0')

    if shift.filter_str and shift.filter_str not in fname:
        return tr_uv

    if shift.type == 'none':
        return tr_uv
    if shift.type == 'ushiftp':
        de1 = f'(u0+v0)/2*{shift.val/1000}/sqrt(2)'
        tr_uv = ( f'u0 + {de1}', f'v0 + {de1}')
    else:
        raise NotImplementedError(f"shift type {shift.type} not implemented.")
    return tr_uv

def create_fit_input(cfg, write=True):
    """ Create the fit input files.

    Args:
        cfg (dict): configuration dictionary that holds all the configuration settings.

    TODO: implement multiprocessing?
    """
    mcfit_table, mkchain = bkg_model.read_jagsh5(cfg.jags_fpath)
    staging_dir = cfg.staging_config['outpath']
    
    fpaths = mcfit_table.fname.apply(lambda fname: f"{staging_dir}/{fname}") #mcfit_table.fname[f"{staging_dir}/{fname}" for fname in ]
    treename = cfg.staging_config['treename']

    shift = cfg.fit_pars.shift

    filters = [ cfg.fit_pars.cuts.common ]
    filters_mc = [ *filters, cfg.fit_pars.cuts.mc_only ] if cfg.fit_pars.cuts.mc_only else filters
    filters_ndbd = [ *filters, cfg.fit_pars.cuts.ndbd_only ] if cfg.fit_pars.cuts.ndbd_only else filters
    filters_data = [ *filters, cfg.fit_pars.cuts.data_only ] if cfg.fit_pars.cuts.data_only else filters

    defs0 =  [('u', 'u0'), ('v', 'v0')]
    def _get_hist_mc(fpath):
        if shift.filter_str in fpath:
            # the def above doesn't get changed right? make sure..
            fname = os.path.basename(fpath)
            tr_uv = get_uv_transform(shift, fname)
            defs = [('u', tr_uv[0]), ('v', tr_uv[1]) ]
        else:
            defs = defs0
        return staging.read_hist(fpath, cfg.hm, treename, defs, filters_mc, rtype='numpy')[0]
    
    mchists = pd.Series(fpaths).apply(_get_hist_mc)
    mchists.name = 'mchist'
    mctable = pd.concat([mcfit_table, mchists], axis=1, join='inner')

    ndbd_mcpath = staging_dir + '/0vbb.root'
    ndbd_name = 'ndbd'
    defs_ndbd = defs0
    ndbd_hist, _, __ = staging.read_hist(ndbd_mcpath, cfg.hm, treename, defs_ndbd, filters_ndbd, rtype='numpy')
    
    signal_df = pd.DataFrame({
        'mean': 0,
        'std': 0,
        'mode': 0,
        'integral': np.sum(ndbd_hist),
        'mchist': [ndbd_hist]
    }, index=[ndbd_name] )

    data_fpath = f"{cfg.data_config['out_fpath']}"
    defs_data = defs0
    data = staging.read_evlist(data_fpath, treename, defs=defs_data, filters=filters_data, rtype='numpy')

    fit_outfname = f'{cfg.fit_pars.fit_dir}/input_{cfg.fit_pars.name}.h5'

    if write:
        hdf.write_input(fit_outfname, mctable, signal_df, data,
                        signal_name=ndbd_name, mkchain=None, overwrite=False)
    
    return mctable, signal_df, data
    
if __name__ == '__main__':
    from qpym2.cfg.common_vars import cfg_m2nat21_may23 as cfg
    create_fit_input(cfg)