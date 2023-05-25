import site
site.addsitedir('/global/homes/x/xcorat/Software//QPyM2/')
# site.addsitedir('/global/homes/x/xcorat/Software//QPyMC/')

from qpym2.io import bkg_model
from qpym2.io import staging
from qpym2.cfg.common_vars import cfg_m2nat21_may23 as cfg

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

def _create_rdf_partial(mcpath, read_config, write_config, filters):
    """ Create a partial function for create_rdf_m2mc to be run with multiprocessing. """
    from functools import partial
    return partial(staging.create_rdf_m2mc, mcpath=mcpath, read_config=read_config,
                   write_config=write_config , filters=filters)

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
        print('starting pool...')
        if ncores is None or ncores == 1:
            ncores = get_ncores()
        with Pool() as pool:
            ret = pool.map(_create_rdf, table.fname)
        print('end pool...', len(ret))
    else:
            ret = [_create_rdf(fname) for fname in table.fname]
            
if __name__ == '__main__':
    """ Create the staging files. """
    create_staging_mp(cfg)