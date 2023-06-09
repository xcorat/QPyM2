import arviz as az
import pandas as pd

def write_fit(table, trace, outdir, fname='fit'):
    # TODO: change to `out_fpath`
    az.to_netcdf(trace, f'{outdir}/{fname}.nc')

    # TODO: Is it ok to store the hists in the dataframe too? just for 5 rows?
    table.to_hdf(f'{outdir}/{fname}.h5', key='components')

def read_fit(outdir=None, fname='fit', fullname=None):
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
    if fullname is None:
        if outdir is None:
            raise ValueError('Either outdir or fullname must be provided')
        fullname = f'{outdir}/{fname}'
    
    trace = az.from_netcdf(fullname + '.nc')
    table = pd.read_hdf(fullname + '.h5', key='components')

    return table, trace