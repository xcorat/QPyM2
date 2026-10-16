import os
import arviz as az
import numpy as np
import pandas as pd

from qpym2.io import hdf

"""
TODO: Add the model configuration to the output files, read/write.
"""


def read_input(input_fname, dataname='data', signal_name='ndbd'):
    """Read fit input HDF5 file.

    Convenience wrapper that reads the background model table, signal MC,
    and data histogram from a single HDF5 input file.

    Args:
        input_fname: Path to input HDF5 file.
        dataname: HDF5 dataset key for the data histogram ('data' or
            'data_blinded' for blind analysis).
        signal_name: HDF5 dataset key for the signal MC histogram.

    Returns:
        Tuple of (mctable, signal_df, data, edges).
    """
    print(f"Reading fit input from {input_fname}")
    mctable, edges = hdf.read_bkg_model(input_fname)
    data = hdf.read_data(input_fname, dataname=dataname)
    signal_df = hdf.read_signal(input_fname, signal_name)

    return mctable, signal_df, data, edges

def write_fit(table, trace, outdir, data=None, fname='fit', group=None):
    """
    TODO: implementation of group name is weird...
    """
    gn = group+"/" if group else ""
    os.makedirs(outdir, exist_ok=True)
    # TODO: change to `out_fpath`
    tf = az.to_netcdf(trace, f'{outdir}/{gn}{fname}.nc')

    # TODO: Is it ok to store the hists in the dataframe too? just for 5 rows?
    cf = table.to_hdf(f'{outdir}/{fname}.h5', key=f'{gn}components')

    if data is not None:
        s = pd.Series([p for p in data])
        s.to_hdf(f'{outdir}/{fname}.h5', key=f'{gn}data')

    return cf, tf

def read_fit(outdir=None, fname='fit', fullname=None, group=None):
    """ Read fit from disk (trace and the components table)

    Args:
        outdir: str
            output directory
        fname: str  
            output file name
        fullname: str
            full path to the output file (without the extension).
            If None, it is constructed from outdir and fname.

    Returns:
        tuple: (component_table, trace)
    """
    gn = group+"/" if group else ""
    if fullname is None:
        if outdir is None:
            raise ValueError('Either outdir or fullname must be provided')
        fullname = f'{outdir}/{gn}{fname}'
    
    trace = az.from_netcdf(fullname + '.nc')
    table = pd.read_hdf(fullname + '.h5', key=f'{gn}components')
    try:
        data = pd.read_hdf(fullname + '.h5', key=f'{gn}data')
    except:
        data = None

    return table, trace, data


