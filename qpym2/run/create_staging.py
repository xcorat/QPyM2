""" Create the staging files. 

TODO:
    + Check segfault errors while running multiprocessing, seem to happen with generic errors.
    + Refactor individual calls to be made from multiprocessing code?
"""

import site
site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/')
# site.addsitedir('/global/homes/x/xcorat/Software//QPyMC/')

from qpym2.utils import debug

from qpym2.io import bkg_model
from qpym2.io import staging
from qpym2.cfg.common_vars import cfg_m2nat21_may23 as cfg

_NDBD_STAGING_NAME = 'ndbd'

def create_staging(cfg):
    """ Create the staging files. 
        
    All the configuration is in the cfg object.    
    """


    table, _ = bkg_model.read_jagsh5(cfg.jags_fpath)
    read_config = cfg.ares_config
    mcpath = read_config['inpath']

    for fname in table.fname:
        staging.create_rdf_m2mc(fname, mcpath,  
                                read_config=read_config,
                                write_config=cfg.staging_config, 
                                filters=[cfg.staging_cuts])

def _create_rdf_partial(**kwargs):
    """ Create a partial function for create_rdf_m2mc to be run with multiprocessing. 
    
    TODO: Use mp.Value and mp.Array to share common varaibles between processes.
    """
    from functools import partial
    return partial(staging.create_rdf_m2mc, **kwargs)

def create_staging_data(cfg):
    """ Create the staging files. 
        
    Args:
        cfg (dict): configuration dictionary that holds all the configuration settings.

    TODO: not tested.
    """
    from multiprocessing import Pool
    from qpym2.utils import get_ncores

    read_config_data = cfg.data_config
    fname_data = read_config_data['fname']

    staging.create_rdf_m2mc(fname_data, 
                            mcpath=read_config_data['inpath'],
                            read_config=read_config_data,
                            write_config=cfg.staging_config, 
                            filters=[cfg.staging_cuts],
                            out_fname=read_config_data['out_fpath'])


def create_staging_mp(cfg, multiproc=True, ncores=None):
    """ Create the staging files. 
        
    Args:
        cfg (dict): configuration dictionary that holds all the configuration settings.
        multiproc (bool): if True, use multiprocessing.
        ncores (int): number of cores to use. If None, use all available cores.
    """
    from multiprocessing import Pool
    from qpym2.utils import get_ncores

    table, _ = bkg_model.read_jagsh5(cfg.jags_fpath)
    read_config = cfg.ares_config
    mcpath = read_config['inpath']

    _create_rdf = _create_rdf_partial(mcpath, read_config, cfg.staging_config, [cfg.staging_cuts])

    if multiproc:
        debug('starting pool...')
        if ncores is None or ncores == 1:
            ncores = get_ncores()
        with Pool() as pool:
            ret = pool.map(_create_rdf, table.fname)
        debug('end pool...', len(ret))
    else:
        ret = [_create_rdf(fname) for fname in table.fname]

def create_staging_all_mp(cfg, multiproc=True, ncores=None):
    """ Create MC, 0v and data staging files in parallel. 
    
    Args:
        cfg (dict): configuration dictionary that holds all the configuration settings.
        multiproc (bool): if True, use multiprocessing.
        ncores (int): number of cores to use. If None, use all available cores.

    TODO: this actually does not run in parallel, we need to replace `map` method with 
        an async method and read the output later.
    """
    from multiprocessing import Pool
    from qpym2.utils import get_ncores
    from os.path import split

    table, _ = bkg_model.read_jagsh5(cfg.jags_fpath)
    read_config = cfg.ares_config
    read_config_ndbd = cfg.ares_config_0v
    read_config_data = cfg.data_config
    
    mcpath = read_config['inpath']
    out_fname_ndbd = f'{cfg.staging_config["outpath"]}/{_NDBD_STAGING_NAME}.root'
    fname_ndbd = read_config_ndbd['fname']

    fname_data = read_config_data['fname']

    cuts_mc = [cfg.staging_cuts]
    cuts_ndbd = [cfg.staging_cuts]
    cuts_data = [cfg.staging_cuts]
    if 'extra_cuts' in read_config:
        cuts_mc = [cfg.staging_cuts, read_config['extra_cuts']]
    if 'extra_cuts' in read_config_ndbd:
        cuts_ndbd = [cfg.staging_cuts, read_config_ndbd['extra_cuts']]
    if 'extra_cuts' in read_config_data:
        cuts_data = [cfg.staging_cuts, read_config_data['extra_cuts']]


    # TODO: can we create the partial directly here? or do we need the top level function 
    #        to work with pickling in multiprocessing?
    _create_rdf = _create_rdf_partial(mcpath=mcpath, read_config=read_config,
                                      write_config=cfg.staging_config, 
                                      filters=cuts_mc)
    # NOTE: we need to pass these as kwargs for the partial to work
    _create_rdf_0v = _create_rdf_partial(mcpath=read_config_ndbd['inpath'],
                                        read_config=read_config_ndbd,
                                        write_config=cfg.staging_config,
                                        filters=cuts_ndbd,
                                        out_fname=out_fname_ndbd)
    _create_rdf_data = _create_rdf_partial(mcpath=read_config_data['inpath'],
                                        read_config=read_config_data,
                                        write_config=cfg.staging_config, 
                                        filters=cuts_data,
                                        out_fname=read_config_data['out_fpath'])

    if multiproc:
        debug('starting pool...')
        if ncores is None or ncores < 1:
            ncores = get_ncores()
        with Pool() as pool:
            # TODO: use an async method instead of map.
            ret_data = pool.map(_create_rdf_data, [fname_data])
            ret_0v = pool.map(_create_rdf_0v, [fname_ndbd])
            ret_mc = pool.map(_create_rdf, table.fname)
        debug('end pool...', len(ret_mc)+len(ret_0v)+len(ret_data))
    else:
        debug('creating staging files without parallelization...')
        ret_mc = [_create_rdf(fname) for fname in table.fname]
        ret_0v = [_create_rdf_0v(fname_ndbd)]
        ret_data = [_create_rdf_data(fname_data)]

def create_staging_0v(cfg):
    """ Create the staging files. 
        
    Args:
        cfg (dict): configuration dictionary that holds all the configuration settings.
    """
    from multiprocessing import Pool
    from qpym2.utils import get_ncores

    read_config_0v = cfg.ares_config_0v
    out_fname = cfg.staging_config['outpath'] + '/0vbb.root'

    # NOTE: we need to pass these as kwargs for the partial to work
    _create_rdf_0v = _create_rdf_partial(mcpath=read_config_0v['inpath'], read_config=read_config_0v,
                                         write_config=cfg.staging_config, filters=[cfg.staging_cuts],
                                         out_fname=out_fname)
    fname_0v = read_config_0v['fname']

    ret = [_create_rdf_0v(fname) for fname in [fname_0v]]

if __name__ == '__main__':
    """ Create the staging files. """
    create_staging_all_mp(cfg)