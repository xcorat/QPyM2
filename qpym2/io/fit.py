import arviz as az
import pandas as pd

def write_fit(table, trace, outdir, data=None, fname='fit', group=None):
    """
    TODO: implementation of group name is weird...
    """
    gn = group+"/" if group else ""
    
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