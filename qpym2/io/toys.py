""" Routines for saving and reading data related to toy fits 

    NOTE: need to be finalized

    The toy runs are saved as a dataframe in a HDF5 file. 
    The dataframe has the following columns:

    idx (int): index of the toy run.
    type (str): type of the toy run. Can be 'fit', 'bias' or 'sens'.
                fit: toy runs with data created using the fit results.
                bias: toy runs with injected 0v when creating data
                sens: toy runs without 0v component when creating data.
    ndbd_inj (float): number of 0vbb events injected (for bias runs).
    logp_mode (float): mode of the logp posterior for the toy fit.
    ...[(fit results) : columns with the fit results.
"""

import pandas as pd

TOY_TYPES = ['fit', 'bias', 'sens']
_STAT_SUFFIXES = ['mode', 'p16', 'p90', 'upp34']
# Metadata columns in batch CSVs (trailing _ stripped on read)
_META_COLS = ['logp_mode_', 'ndbd_inj_', 'type_', 'task_id',
              'co_cut', 'ndbd_cutoff', 'co_shift']


def _unflatten_columns(df):
    """Convert flat CSV columns (e.g. ``ndbd_mode``) back to MultiIndex ``(var, stat)``.

    Metadata columns (``logp_mode_``, ``ndbd_inj_``, etc.) are kept as
    plain single-level columns with the trailing underscore stripped.

    Returns:
        pd.DataFrame with MultiIndex stat columns + plain metadata columns.
    """
    meta = {}
    stat_data = {}
    for col in df.columns:
        if col in _META_COLS:
            clean = col.rstrip('_')
            meta[clean] = df[col].values
        else:
            # Try to split "var_stat" on known suffixes
            matched = False
            for suf in _STAT_SUFFIXES:
                if col.endswith(f'_{suf}'):
                    var = col[: -(len(suf) + 1)]
                    stat_data[(var, suf)] = df[col].values
                    matched = True
                    break
            if not matched:
                meta[col] = df[col].values

    out = pd.DataFrame(stat_data)
    out.columns = pd.MultiIndex.from_tuples(out.columns)
    for k, v in meta.items():
        out[k] = v
    return out


def read_toys_csv(csv_path, type=None):
    """Read toy results from a batch CSV file.

    Args:
        csv_path (str): Path to the CSV file.
        type (str | None): Filter by toy type ('fit', 'bias', 'sens').
            If None, all rows are returned.

    Returns:
        pd.DataFrame with MultiIndex stat columns, matching the format
        returned by :func:`read_toys`.
    """
    df = pd.read_csv(csv_path)
    if type is not None:
        df = df[df['type_'] == type].reset_index(drop=True)
    return _unflatten_columns(df)


def write_toys(toy_df, outdir, type, fname='toys'):
    """ write toy runs to a file. 

    type can be 'fit', 'bias' or 'sens'. 
    
    """
    from qpym2.io.hdf import _write_df
    if type in TOY_TYPES:
        outpath = outdir + f'/{fname}.h5'
        _write_df(toy_df, outpath, key=type, overwrite=True)
    else:
        raise ValueError(f'Unknown type \'{type}\'. Type must be one of \'fit\', \'bias\' or \'sens\'.')
    

def read_toys(outdir, type, fname='toys'):
    """ read toy runs from a file. """

    if type in TOY_TYPES:
        inpath = outdir + f'/{fname}.h5'
        key = type
        return pd.read_hdf(inpath, key=key)
    else:
        raise ValueError(f'Unknown type \'{type}\'. Type must be one of \'fit\', \'bias\' or \'sens\'.')    

