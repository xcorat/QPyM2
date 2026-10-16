""" Routines for saving and reading fit runs for systematics studies """

import pandas as pd
from qpym2.io.hdf import _write_df

SYST_TYPES = ['bm_stat', 'lss', 'ndbd_shift']

def write_syst(df, outdir, type, fname='syst'):
    if type not in SYST_TYPES:
        raise ValueError(f'Invalid type {type}. Must be one of {SYST_TYPES}')
    
    outpath = outdir + f'/{fname}.h5'
    _write_df(df, outpath, key=type, overwrite=True)


def read_syst(outdir, type, fname='syst'):
    """ read ststematics run results from a file. """
    if type not in SYST_TYPES:
        raise ValueError(f'Invalid type {type}. Must be one of {SYST_TYPES}')
    
    inpath = outdir + f'/{fname}.h5'
    return pd.read_hdf(inpath, key=type)