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

def get_uv_transform(shift, fname, invars=('u0', 'v0')):
    """ Get the uv transform for the shift type and value.

    Args:
        shift (qpym2.cfg.fit_pars.Shift): shift
        fname (str): file name
        invars (tuple): input (u, v) to transform
    Returns:
        tuple: (u, v) transform
    """

    if shift.filter_str and shift.filter_str not in fname:
        return invars
    if shift.type == 'none':
        return invars
    
    u, v = invars   
    if shift.type == 'ushiftp':
        de1 = f'({u}+{v})/2*{shift.val/1000}/sqrt(2)'
        tr_uv = ( f'{u} + {de1}', f'{v} + {de1}')
    else:
        raise NotImplementedError(f"shift type {shift.type} not implemented.")
    return tr_uv

def get_uv_transform_shift(shift, invars=('u0', 'v0')):
    """ Get the uv transform for the shift type and value.

    Args:
        shift (qpym2.cfg.fit_pars.Shift): shift
        invars (tuple): input (u, v) to transform

    Returns:
        tuple: (u, v) transform
    """
    u, v = invars    

    if shift.type == 'none':
        return invars
    
    if shift.type == 'ushiftp':
        de1 = f'({u}+{v})/2*{shift.val/1000}/sqrt(2)'
        tr_uv = ( f'{u} + {de1}', f'{v} + {de1}')
        return tr_uv
    else:
        raise NotImplementedError(f"shift type {shift.type} not implemented.")
    
def get_uv_transform_dq(dq_pars, invars=('u0', 'v0')):
    """ Get the uv transform for the shift type and value.

    Args:
        **dq_pars (tuple): dq_pars
        invars (tuple): input (u, v) to transform

    Returns:
        tuple: (u, v) transform
    """
    u, v = invars
    e1,e2 = f'({u}+{v})/sqrt(2)', f'({u}-{v})/sqrt(2)'  
    a, b, c = dq_pars

    tru = f'{u} + {a}*sqrt(2) + {b}*({u}) + {c}/sqrt(2)*(({u})*({u})+({v})*({v}))'
    trv = f'{v}               + {b}*({v}) + {c}/sqrt(2)*({u})*({v})'

    return (tru, trv)

    # if lss_defs.type == 'none':
    #     return invars
    
    # if shift.type == 'ushiftp':
    #     de1 = f'({u}+{v})/2*{shift.val/1000}/sqrt(2)'
    #     tr_uv = ( f'{u} + {de1}', f'{v} + {de1}')
    #     return tr_uv
    
    # if shift.type == 'lss_residual':
    #     de1 = f'({u}+{v})/2*{shift.val/1000}/sqrt(2)'
    #     tr_uv = ( f'{u} + {de1}', f'{v} + {de1}')
    #     return tr_uv
    
    # else:
    #     raise NotImplementedError(f"shift type {shift.type} not implemented.")

def read_ls_syst(cfg):
    """
    """
    def _fnpars(fn, skip_last=True):
        p = [fn.GetParameter(i) for i in range(fn.GetNpar()) ]
        if skip_last: return p[:-1]
        else: return p

    dpath = '/global/homes/x/xcorat/Software/cuore-nersc-modern_v4.0.0/cuoremc/ares/data/test'
    from qpym2.io.qls_scaling import read_scaling_pars, read_lspars, sample_pars
    datasets = range(3601, 3616, 1)
    fns_Q, covs_Q, fns_sigma, covs_sigma = read_scaling_pars(
        [f'{dpath}/residual_and_width_vs_energy_ds{ds}.root' for ds in datasets]
    )

    fname_ls = f'{dpath}/PeakShape_ds*.root'
    dfps = read_lspars(fname_ls)

    qpars0 = [_fnpars(fn) for fn in fns_Q]
    qpars_1l = [sample_pars(fn, cov, skip_last_par=True, nsamples=1, sample_point=1) for fn, cov in zip(fns_Q, covs_Q)]
    qpars_1r = [sample_pars(fn, cov, skip_last_par=True, nsamples=1, sample_point=-1) for fn, cov in zip(fns_Q, covs_Q)]

    sigma_pars0 = [_fnpars(fn, skip_last=False) for fn in fns_sigma]
    sigma_pars_1l = [sample_pars(fn, cov, skip_last_par=False, nsamples=1, sample_point=1) for fn, cov in zip(fns_sigma, covs_sigma)]
    sigma_pars_1r = [sample_pars(fn, cov, skip_last_par=False, nsamples=-1, sample_point=1) for fn, cov in zip(fns_sigma, covs_sigma)]

    import pandas as pd
    qpars_df = pd.DataFrame({'q_pars': qpars0,
                             'q_pars_sample_1l': qpars_1l,
                             'q_pars_sample_1r': qpars_1r,
                              'siqma_pars': sigma_pars0,
                              'sigma_pars_sample_1l': sigma_pars_1l,
                              'sigma_pars_sample_1r': sigma_pars_1r,
                            },
                            index=datasets).rename_axis('ds')
    
    return qpars_df, dfps

from qpym2.io.staging_transforms import residual_sampling_uv as rsuv_tansform

def create_fit_input(cfg, signal_name='ndbd', syst_test_ds=None, shift=None, write=True):
    mcfit_table, mkchain = bkg_model.read_jagsh5(cfg.jags_fpath)
    staging_dir = cfg.staging_config['outpath']

    fpaths = mcfit_table.fname.apply(lambda fname: f"{staging_dir}/{fname}") #mcfit_table.fname[f"{staging_dir}/{fname}" for fname in ]
    treename = cfg.staging_config['treename']

    # don't use the default, remove this.
    if shift is None:
        shift = cfg.fit_pars.shift

    filters = [ cfg.fit_pars.cuts.common ]
    filters_mc = [ *filters, cfg.fit_pars.cuts.mc_only ] if cfg.fit_pars.cuts.mc_only else filters
    filters_ndbd = [ *filters, cfg.fit_pars.cuts.ndbd_only ] if cfg.fit_pars.cuts.ndbd_only else filters
    filters_data = [ *filters, cfg.fit_pars.cuts.data_only ] if cfg.fit_pars.cuts.data_only else filters

    _U, _V = 'u0', 'v0'
    _LSS_DEFS = ('ulss', 'vlss')

    # add the definitions depending on the parameter set for the run...
    def _get_defs(fpath, syst_test_ds, shift, inpars=('u0', 'v0')):
        _u, _v = inpars
        defs = []
        if syst_test_ds is not None:
            ds = syst_test_ds['ds']
            # defs = [*defs, (_LSS_DEFS[0], f'ds=={syst_test_ds}? {qpars_df["udef"][syst_test_ds["ds"]]}:{_u}'),
            #                (_LSS_DEFS[1], f'ds=={syst_test_ds}? {qpars_df["vdef"][syst_test_ds["ds"]]}:{_v}')]
            defs = [(_LSS_DEFS[0], syst_test_ds["udef"]),
                    (_LSS_DEFS[1], syst_test_ds["vdef"])],
            
            _u, _v = _LSS_DEFS

        if shift.filter_str in fpath:
            # the def above doesn't get changed right? make sure..
            # fname = os.path.basename(fpath)
            tr_uv = get_uv_transform_shift(shift, invars=(_u, _v))
            print(defs)
            defs = [*defs, ('u', tr_uv[0]), ('v', tr_uv[1]) ]
            print(defs)
        else:
            defs = [*defs, ('u', _u), ('v', _v)]

        return defs
    
    comp_defs = pd.Series(fpaths).apply(
        lambda fpath: _get_defs(fpath, syst_test_ds, shift, inpars=(_U, _V)))
    comp_defs.name = 'defs'
    # return comp_defs
    
    for comp_def in comp_defs:
        for def_ in comp_def:
            pass
            # print('--', def_)
    mchists = [staging.read_hist(fpath, cfg.hm, treename, defs, filters_mc, rtype='numpy')[0]
            for fpath, defs in zip(fpaths, comp_defs)]
    

    mctable = pd.concat([ mcfit_table,
                        comp_defs,
                        pd.Series(mchists, name='mchist', index=mcfit_table.index)
                        ], axis=1, join='inner')
    
    # Read the signal info -----------------------------------------------------
    ndbd_mcpath = f'{staging_dir}/{signal_name}.root'
    defs_ndbd = _get_defs(ndbd_mcpath, syst_test_ds, shift, inpars=(_U, _V))
    ndbd_hist, _, __ = staging.read_hist(ndbd_mcpath, cfg.hm, treename, defs_ndbd, filters_ndbd, rtype='numpy')

    signal_df = pd.DataFrame({
        'mean': 0,
        'std': 0,
        'mode': 0,
        'integral': np.sum(ndbd_hist),
        'mchist': [ndbd_hist]
    }, index=[signal_name] )

    # Read the data -----------------------------------------------------------
    data_fpath = f"{cfg.data_config['out_fpath']}"
    defs_data = (('u', _U), ('v', _V))
    data = staging.read_evlist(data_fpath, treename, defs=defs_data, filters=filters_data, rtype='numpy')

    return mctable, signal_df, data

def create_fit_input_old(cfg, write=True):
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

    qpars0, qpars, dfps = read_ls_syst(cfg.lspars_config)
    # dq_pars = qpars - qpars0

    # --------

    def _get_hist_mc(fpath, default_axes=('u0', 'v0')):
        # TODO: TODO: untangle this. The order of transformations matter, and we need to keep 
        #           track of what transformations are to be applied and what their output
        #           variables (u,v) are.
        defs = []


        if shift.filter_str in fpath:
            # the def above doesn't get changed right? make sure..
            fname = os.path.basename(fpath)
            tr_uv = get_uv_transform(shift, fname, invars=defaults)
            defs = [*defs, ('u', tr_uv[0]), ('v', tr_uv[1]) ]
        else:
            defs = [*defs, ('u', defaults[0]), ('v', defaults[1])]

        for def_ in defs:
            print(def_)

        return staging.read_hist(fpath, cfg.hm, treename, defs, filters_mc, rtype='numpy')[0]
    

    _LSS_DEFS = ('ulss', 'vlss')
    if 'lss' in cfg.fit_pars.__dir__() and cfg.fit_pars.lss is not None:
        truv_lss = rsuv_tansform(qpars, qpars0, invars=defaults)
        defs = [(_LSS_DEFS[0], truv_lss[0]), (_LSS_DEFS[1], truv_lss[1])]  
        # IMPORTANT: change the current u,v values, and use them for the next transformation.
        defaults = _LSS_DEFS

    mchists = pd.Series(fpaths).apply(_get_hist_mc)
    mchists.name = 'mchist'
    mctable = pd.concat([mcfit_table, mchists], axis=1, join='inner')

    # ndbd_mcpath = staging_dir + '/0vbb.root'
    # ndbd_name = 'ndbd'
    # defs_ndbd = defaults
    # ndbd_hist, _, __ = staging.read_hist(ndbd_mcpath, cfg.hm, treename, defs_ndbd, filters_ndbd, rtype='numpy')
    
    # signal_df = pd.DataFrame({
    #     'mean': 0,
    #     'std': 0,
    #     'mode': 0,
    #     'integral': np.sum(ndbd_hist),
    #     'mchist': [ndbd_hist]
    # }, index=[ndbd_name] )

    # data_fpath = f"{cfg.data_config['out_fpath']}"
    # defs_data = defaults
    # data = staging.read_evlist(data_fpath, treename, defs=defs_data, filters=filters_data, rtype='numpy')

    # fit_outfname = f'{cfg.fit_pars.fit_dir}/input_{cfg.fit_pars.name}.h5'

    # if write:
    #     hdf.write_input(fit_outfname, mctable, signal_df, data,
    #                     signal_name=ndbd_name, mkchain=None, overwrite=False)
    
    # return mctable, signal_df, data
    
if __name__ == '__main__':
    from qpym2.cfg.common_vars import cfg_m2nat21_may23 as cfg
    create_fit_input(cfg)