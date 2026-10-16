""" Create fit input related files for M2 u-v axis related inputs.



"""
# TODO: we shouldn't need to add the qpym2 path here if the module is
# installed properly. But for now, we need to add it to the path.
# site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/')

import os
import site
import numpy as np
import pandas as pd


from qpym2.io import bkg_model, staging, hdf
from qpym2.utils.logging import debug

from qpym2.utils import TempConfig

# TODO: how to keep global config variables?
"""
input filename config?
"""
_INPUT_NAME_ = 'input'
_CONAME = 'co60'
_NDBDNAME = 'ndbd'

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
    
    if shift.type == 'ushiftp':
        tr_uv = get_ushiftp_uvtr(shift.val, invars)
    else:
        raise NotImplementedError(f"shift type {shift.type} not implemented.")
    return tr_uv

def get_ushiftp_uvtr(shift: float, invars=('u0', 'v0')):
    """ Get the uv transform for the shift type and value.

    Args:
        shift (float): shift
        invars (tuple): input (u, v) to transform
    Returns:
        tuple: (u, v) transform
    """
    u, v = invars   
    de1 = f'({u}+{v})/2*{shift/1000}/sqrt(2)'
    tr_uv = ( f'{u} + {de1}', f'{v} + {de1}')
        
    return tr_uv

def get_uv_transform_shift(shift, invars=('u0', 'v0')):
    """ Get the uv transform for the shift type and value.

    Args:
        shift (qpym2.cfg.fit_pars.Shift): shift
        invars (tuple): input (u, v) to transform

    Returns:
        tuple: (u, v) transform

    TODO: this should just take `shift_type`, `shift_values` instead of a cfg obj. 
    """
    u, v = invars    

        
    if shift.type == 'ushiftp':
        de1 = f'({u}+{v})/2*{shift.val/1000}/sqrt(2)'
        tr_uv = ( f'{u} + {de1}', f'{v} + {de1}')
        return tr_uv
    elif shift.type == 'none':
        return invars
    
    else:
        raise NotImplementedError(f"shift type {shift.type} not implemented.")
    
from qpym2.io import staging
from qpym2.utils.hist_model import get_hm_cuts

def _cache_key(fpath, defs, filters):
    """ Build a hashable cache key from (fpath, defs, filters). """
    return (fpath,
            tuple(tuple(d) for d in defs),
            tuple(f for f in filters))

def create_evlist(fit_cfg, staging_cfg, inpars=('u0', 'v0'), hist_cache=None):
    """ Create the data event list for the fit. """
    data_fpath = f"{staging_cfg['outpath']}/{staging_cfg['data_fname']}"
    debug('datapath:', data_fpath)
    _u, _v = inpars
    defs_data = (('u', _u), ('v', _v))
    cuts = fit_cfg['cuts']

    ecuts = get_hm_cuts(fit_cfg['h2uv_model'])

    filters = [ ecuts, cuts['common'] ]
    if cuts.get('data_only'):
        filters = [ *filters, fit_cfg['cuts']['data_only'] ]

    key = _cache_key(data_fpath, defs_data, filters) if hist_cache is not None else None
    if key is not None and key in hist_cache:
        debug('cache hit (evlist):', data_fpath)
        return hist_cache[key]

    data = staging.read_evlist(data_fpath, staging_cfg['treename'],
                               defs=defs_data,
                               filters=filters,
                               rtype='numpy')

    if key is not None:
        hist_cache[key] = data

    return data

# TODO: merge? or ?
def get_shift(val, type='ushiftp', filter_str=None):
    return TempConfig(val=val, type=type, filter_str=filter_str)

def get_shiftname(shift):
    return f"{shift.type}{shift.val*10:02.0f}"

def get_systname(syst_pars):
    name = ''
    if 'ndbd_shift' in syst_pars:
        name = _NDBDNAME + get_shiftname(get_shift(syst_pars['ndbd_shift']))

    return name

def get_uv_transform(shift, fname, invars=('u0', 'v0'), outvars=('u', 'v')):
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

def get_fit_transformers(fit_cfg, fit_pars={}, syst_pars={}, invars=('u0', 'v0'), outvars=('u', 'v')):
    """ Return a function that creates the transformations for each fit setting
        according to the fit configurations.
        
    Args:
        fit_cfg (dict): fit configuration
        fit_pars (dict): fit parameters
        syst_pars (dict): systematic parameters
        invars (tuple): input (u, v) to transform
        outvars (tuple): output (u, v) to transform

    Returns:
        function: function that creates the transformations for each fit setting

    TODO: move fit_pars and syst_cfg to fit_cfg? why not there again?
    """

    # NOTE: the only transforms that are implemented now are the shifts for mc and ndbd

    _u, _v = outvars
    if not fit_pars: fit_pars = {}
    if not syst_pars: syst_pars = {}

    default_defs = ( (_u, invars[0]), (_v, invars[1]) )

    if 'co_shift' in fit_pars:
        shift = get_shift(val=fit_pars['co_shift'], filter_str=_CONAME)
        mc_defs = list(zip( (_u, _v), get_ushiftp_uvtr(shift.val, invars) ))
    else:
        shift = None
    
    if getattr(syst_pars, 'ndbd_shift', None):
        ndbd_shift = get_shift(val=syst_pars.get('ndbd_shift'), filter_str=_NDBDNAME)
        ndbd_defs = list(zip( (_u, _v), get_ushiftp_uvtr(ndbd_shift.val, invars) ))
    else:
        ndbd_shift = None

    # if no shift is specified, return the identity transform
    if shift is None and ndbd_shift is None:
        return lambda fpath: default_defs
    
    # TODO optimize!
    # NOTE: this is the logic of all the systematic transforms, so we
    #       should add more here when needed.
    def _get_transforms(fpath):
        if ndbd_shift != None and \
            ndbd_shift.val != 0 and \
            _NDBDNAME in fpath:
            return ndbd_defs
        
        elif shift != None and \
             shift.val != 0 and \
             _CONAME in fpath:
            return mc_defs
        
        else:
            return default_defs
        
    return _get_transforms

def _cached_read_hist(fpath, hm, treename, defs, filters, hist_cache=None):
    """ Read a hist with optional caching.

    Returns:
        tuple: (hist, edges) where edges is a tuple of edge arrays (xedges, yedges, ...).
    """
    key = _cache_key(fpath, defs, filters) if hist_cache is not None else None
    if key is not None and key in hist_cache:
        debug('cache hit (hist):', fpath)
        return hist_cache[key]

    read_result = staging.read_hist(fpath, hm, treename, defs, filters, rtype='numpy')
    result = (read_result[0], tuple(read_result[1:]))

    if key is not None:
        hist_cache[key] = result

    return result

def get_mctable(fit_cfg, staging_cfg, jags_fpath, 
                fit_pars={}, syst_pars={}, invars=('u0', 'v0'), hist_cache=None):
    """ TODO: This need to be refactored """
    ndbd_table = {  'fname': 'ndbd.root',
                    'mean': 0,
                    'std': 0,
                    'mode': 0,}

    mcfit_table, _ = bkg_model.read_jagsh5(jags_fpath)

    mctable = pd.concat([mcfit_table, pd.DataFrame([ndbd_table], index=['ndbd'])])
    
    fpaths = mctable.fname.apply(lambda fname: f"{staging_cfg['outpath']}/{fname}")

    # TODO: NOTE: does the ecuts applied somewhere else? or should it be?
    filters = [get_hm_cuts(fit_cfg['h2uv_model'], invars)]
    
    cuts = fit_cfg['cuts'] 
    filters = [*filters, cuts['common']]

    # TODO: Add all filters and remove the empty ones later? Or is this the right way?
    if cuts.get('mc_only'):
        filters_mc = [ *filters, fit_cfg['cuts']['mc_only'] ]
    else:
        filters_mc = filters

    if cuts.get('ndbd_only'):
        filters_ndbd = [ *filters, fit_cfg['cuts']['ndbd_only'] ]
    else:
        filters_ndbd = filters

    # add the definitions depending on the parameter set for the run...
    comp_defs = pd.Series(fpaths).apply(
        get_fit_transformers(fit_cfg, fit_pars, syst_pars, invars) )
    comp_defs.name = 'defs'

    # MC components only
    hist_edges_pairs = [_cached_read_hist(fpath, fit_cfg['h2uv_model'], staging_cfg['treename'],
                                          defs, filters_mc, hist_cache)
                        for fpath, defs in zip(fpaths[:-1], comp_defs[:-1])]

    # NDBD only
    hist_edges_pairs.append(_cached_read_hist(fpaths.iloc[-1], fit_cfg['h2uv_model'],
                                              staging_cfg['treename'], comp_defs.iloc[-1],
                                              filters_ndbd, hist_cache))

    mchists = [h for h, _ in hist_edges_pairs]
    edges_list = [e for _, e in hist_edges_pairs]

    mctable['defs'] = comp_defs
    mctable['integral'] = [h.sum() for h in mchists]
    mctable['mchist'] = mchists
    mctable['edges'] = edges_list

    return mctable

def create_input(fit_cfg, input_cfg, jags_fpath,
                 fit_name=None, write_cfg=None, invars=('u0', 'v0'),
                 hist_cache=None):
    fit_pars = fit_cfg['fit_pars']

    out_fname = _INPUT_NAME_ 
    if fit_cfg['h2uv_model'].skip_roi:
        out_fname = _INPUT_NAME_ + '_noroi'

    syst_pars = fit_cfg.get('syst_pars', None)

    staging_comps = get_mctable(fit_cfg, input_cfg, jags_fpath, 
                    fit_pars, syst_pars, invars=invars,
                    hist_cache=hist_cache)

    # Validate all components share identical edges before extracting
    ref_edges = staging_comps['edges'].iloc[0]
    for idx, e in staging_comps['edges'].items():
        for dim, (ea, eb) in enumerate(zip(e, ref_edges)):
            if not np.array_equal(ea, eb):
                raise ValueError(f"Edge mismatch for component '{idx}', dim {dim}")
    edges = list(ref_edges)

    mctable = staging_comps.loc[staging_comps.index != 'ndbd']
    signal_df = staging_comps.loc[staging_comps.index == 'ndbd']
    data = create_evlist(fit_cfg, input_cfg, hist_cache=hist_cache)

    bfn = fit_cfg.get('blinded', None)
    if bfn:
        # Add blinded data to the fit
        if not isinstance(bfn, str): 
            bfn = input_cfg.get('data_blinded_name', 'data_blinded')
        
        bs_cfg = input_cfg.copy()
        bs_cfg['data_fname'] = f'{bfn}.root'
        b_data = create_evlist(fit_cfg, bs_cfg, hist_cache=hist_cache)

    else: b_data = None

    if fit_name is None:
        fit_name = fit_cfg.get('dir_name', fit_cfg['name'])
    #     fitname = get_shiftname(get_shift(fit_pars['co_shift']))
    # if syst_pars is not None:
    #     fitname += '_' + get_systname(syst_pars)

    shift_suffix = fit_cfg.get('shift_suffix', '')

    if write_cfg:
        # outdir = f"{write_cfg['outdir']}/input_{fit_cfg['name']}"
        outdir = f"{write_cfg['outdir']}/{fit_name}"
        os.makedirs(outdir, exist_ok=True)
        fit_outfname = f"{outdir}/{out_fname}{shift_suffix}.h5"

        hdf.write_input(fit_outfname, mctable, signal_df, data, blinded_data=b_data,
                        signal_name='ndbd', mkchain=None, edges=edges,
                        overwrite=write_cfg.get('overwrite', False))
        
    return mctable, signal_df, data, b_data

def create_input_ndbd(fit_cfg, staging_cfg, syst_pars=None,
                 fit_name=None, write_cfg=None, invars=('u0', 'v0')):
    
    fit_pars = fit_cfg['fit_pars']

    out_fname = _INPUT_NAME_ 
    if fit_cfg['h2uv_model'].skip_roi:
        out_fname = _INPUT_NAME_ + '_noroi'
    

    """ TODO: This need to be refactored """
    fname = 'ndbd.root'
    ndbd_table = {  'fname': fname,
                    'mean': 0,
                    'std': 0,
                    'mode': 0,}
    fpath = f"{staging_cfg['outpath']}/{fname}"

    # TODO: NOTE: does the ecuts applied somewhere else? or should it be?
    filters = [get_hm_cuts(fit_cfg['h2uv_model'], invars)]
    
    cuts = fit_cfg['cuts'] 
    filters = [*filters, cuts['common']]

    if cuts.get('ndbd_only'):
        filters_ndbd = [ *filters, fit_cfg['cuts']['ndbd_only'] ]
    else:
        filters_ndbd = filters

    ndbd_defs = get_fit_transformers(fit_cfg, fit_pars, syst_pars, invars)(fpath)
    
    # NDBD only
    hist, *edge_arrs = staging.read_hist(fpath, fit_cfg['h2uv_model'],
                                         staging_cfg['treename'], ndbd_defs,
                                         filters_ndbd, rtype='numpy')

    signal_dict = {
        **ndbd_table,
        'defs': ndbd_defs,
        'integral': hist.sum(),
        'mchist': hist,
        'edges': tuple(edge_arrs),
    }

    signal_df = pd.DataFrame([signal_dict], index=['ndbd'])
    
    # staging_comps.loc[staging_comps.index == 'ndbd']


    if fit_name is None:
        fit_name = fit_cfg.get('dir_name', fit_cfg['name'])
    #     fitname = get_shiftname(get_shift(fit_pars['co_shift']))
    # if syst_pars is not None:
    #     fitname += '_' + get_systname(syst_pars)

    shift_suffix = fit_cfg.get('shift_suffix', '')

    if write_cfg:
        # outdir = f"{write_cfg['outdir']}/input_{fit_cfg['name']}"
        ow = write_cfg.get('overwrite', False)
        outdir = f"{write_cfg['outdir']}/{fit_name}"
        os.makedirs(outdir, exist_ok=True)
        fit_outfname = f"{outdir}/{out_fname}{shift_suffix}.h5"

        hdf._write_signal(signal_df, fit_outfname, signal_name='ndbd', overwrite=ow)

    # dummy return
    return signal_df
