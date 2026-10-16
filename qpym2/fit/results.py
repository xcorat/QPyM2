import pandas as pd
import numpy as np
import itertools

from qpym2.cfg.model_utils import build_model_name, build_model_dir, get_shift_suffix
from qpym2.io import fit as fit_io
from qpym2.utils.mathy import get_mode
from qpym2.utils.logging import debug

_VARS_COL_NAME_ = 'var'

def _get_model_pars(model_coords):
    """ Get the list of model parameters from the model coordinates."""
    mvals = itertools.product(*model_coords.values()) # model values
    return [dict(zip(model_coords.keys(), v)) for v in mvals] # model parameters


def _get_stats(x):
    return {
        'mode': get_mode(x),
        'mean': float(np.mean(x)),
        'p16': np.percentile(x, 16),
        'p90': np.percentile(x, 90),
        'upp34':  np.percentile(x[x>=get_mode(x)], 34)
    }

def read_fit(parent_dir, model_pars, hm=None):
    """ Read a fit from the fit directory.
    
    Args:
        parent_dir: Base fit directory
        model_pars: Dict with co_cut, ndbd_cutoff, co_shift, ndbd_shift, smooth
        hm: Histogram model (deprecated, kept for backward compatibility)
    """
    model_dir = build_model_dir(
        co_cut=model_pars.get('co_cut'),
        ndbd_cutoff=model_pars.get('ndbd_cutoff'),
        smooth=model_pars.get('smooth', True)
    )
    shift_suffix = get_shift_suffix(
        co_shift=model_pars.get('co_shift', 0.0),
        ndbd_shift=model_pars.get('ndbd_shift', 0.0),
    )
    fdir = f'{parent_dir}/{model_dir}'
    fname = f'fit{shift_suffix}'

    return fit_io.read_fit(fdir, fname=fname)

def read_input_data(parent_dir, model_pars):
    """ Read input data from the fit input HDF5 file.
    
    Mirrors read_fit() API: constructs input file path from model_pars internally.
    
    Args:
        parent_dir: Base fit directory
        model_pars: Dict with co_cut, ndbd_cutoff, co_shift, ndbd_shift, smooth
    
    Returns:
        Tuple of (mctable, signal_df, data) where:
            - mctable: DataFrame with background model info
            - signal_df: DataFrame with signal MC info
            - data: (N, 2) array of event coordinates [u, v]
    """
    model_dir = build_model_dir(
        co_cut=model_pars.get('co_cut'),
        ndbd_cutoff=model_pars.get('ndbd_cutoff'),
        smooth=model_pars.get('smooth', True)
    )
    input_path = f'{parent_dir}/{model_dir}/input.h5'
    
    return fit_io.read_input(input_path)

def read_ndbd_eff(parent_dir, model_pars, ndbd_gen=1e8):
    comps, _, _ = read_fit(parent_dir, model_pars, hm=None)
    return comps.integral.ndbd/ndbd_gen

def read_posterior(parent_dir, model_pars, var='ndbd'):
    """ Read the posterior of a fit from the fit directory for the given variable"""

    _, t, __ = read_fit(parent_dir, model_pars)
    return t.posterior[var].values.flatten() 

def read_posterior_stats(parent_dir, model_pars, hm):
    """  Read the posterior stats from a fit directory.

    Args:
        parent_dir (str): The fit directory where each model directory is.
        model_pars (dict): The model parameters.
        include_trace (str): The name of the param to include in the trace of.

    Returns:
        (pd.DataFrame): The posterior stats.
        the DataFrame has the following columns:
            - mode: The mode of the posterior.
            - p16: The 16th percentile of the posterior.
            - p90: The 90th percentile of the posterior.
            - upp34: The 34th percentile of the upper half of the posterior.
            - trace: The trace of the param if include_trace is not None.
    """
    # We only care about logp for now.

    mctable, trace, __ = read_fit(parent_dir, model_pars, hm)
    
    stats = {}
    for var in trace.posterior.data_vars:
        stats[var] = _get_stats(trace.posterior[var].values.flatten())
        stats[var]['prior_norm'] = mctable['integral'][var]
    
    stats['ndbd']['prior_norm'] = 0
    
    stats['logp'] = _get_stats(trace.log_likelihood.pik.values.flatten())
    # TODO: we probably should save the prior_norms with the input,
    #   but for now, manually set the 0vbb to 0

    df = pd.DataFrame(stats).T
    df.index.name = _VARS_COL_NAME_

    debug(stats)
    return df

def read_posterior_stats_wrong(parent_dir, model_pars, hm):
    """  Read the posterior stats from a fit directory.

    Args:
        parent_dir (str): The fit directory where each model directory is.
        model_pars (dict): The model parameters.
        include_trace (str): The name of the param to include in the trace of.

    Returns:
        (pd.DataFrame): The posterior stats.
        the DataFrame has the following columns:
            - mode: The mode of the posterior.
            - p16: The 16th percentile of the posterior.
            - p90: The 90th percentile of the posterior.
            - upp34: The 34th percentile of the upper half of the posterior.
            - trace: The trace of the param if include_trace is not None.
    """
    # We only care about logp for now.

    mctable, trace, __ = read_fit(parent_dir, model_pars, hm)
    
    stats = {}
    for var in trace.posterior.data_vars:
        stats[var] = _get_stats(trace.posterior[var].values.flatten())
        stats[var]['prior_norm'] = mctable['integral'][var]
    
    stats['ndbd']['prior_norm'] = 0
    
    stats['logp'] = _get_stats(trace.log_likelihood.pik.values.flatten())
    # TODO: we probably should save the prior_norms with the input,
    #   but for now, manually set the 0vbb to 0

    df = pd.DataFrame(stats).T
    df.index.name = _VARS_COL_NAME_

    debug(stats)
    return df

def read_stats_dataset(parent_dir, model_coords, hm):
    """ Read the stats dataset from a fit directory.
    
    Args:
        parent_dir (str): The fit directory where each model directory is.
        model_coords (dict): The model coordinates where the key is the 
            coordinate name and the value is a list of coordinate values.
        include_trace (str): The name of the param to include in the trace of.

    Returns:
        (xarray.Dataset): The stats dataset
    
    """
    mpars = _get_model_pars(model_coords)

    # Read all the fit models and create their stats dataframes.
    stat_dfs = [ read_posterior_stats(parent_dir, m, hm) 
                    for m in mpars ]

    for df, m in zip(stat_dfs, mpars):
        for k,v in m.items():
            df[k] = v

    sss = pd.concat(stat_dfs).reset_index().set_index(['ndbd_cutoff',	'co_cut',	'co_shift', _VARS_COL_NAME_])
    # debug(sss.rename(index={np.NaN:'None'}))

    # stat_xr = pd.concat(stat_dfs).reset_index() \
    #             .set_index(['ndbd_cutoff',	'co_cut',	'co_shift', _VARS_COL_NAME_])\
    #             .to_xarray()

    stat_xr = sss.rename(index={np.nan:'None'}).to_xarray()
    return stat_xr

from qpym2.io.toys import read_toys, read_toys_csv, TOY_TYPES
def read_toy_df(fit_dir, model_pars, hm=None, type=None, cols=None):
    """ Read the toy fits from the fit directory.
    
    Args:
        fit_dir (str): The fit directory where each model directory is.
        model_pars (dict): The model parameters (co_cut, ndbd_cutoff, co_shift, etc.).
        hm: Histogram model (deprecated, kept for backward compatibility)
        type (str): The type of toy to read. If None, all toy types are read.
        cols (list): The columns to include in the DataFrame. If None, all columns are included.

    Returns:
        (pd.DataFrame): The toy fits.

    """
    model_dir = build_model_dir(
        co_cut=model_pars.get('co_cut'),
        ndbd_cutoff=model_pars.get('ndbd_cutoff'),
        smooth=model_pars.get('smooth', True)
    )
    shift_suffix = get_shift_suffix(
        co_shift=model_pars.get('co_shift', 0.0),
        ndbd_shift=model_pars.get('ndbd_shift', 0.0),
    )
    fdir = f'{fit_dir}/{model_dir}'
    if type:
        toyfits = read_toys(fdir, type=type)
    else:
        toyfits = pd.concat([read_toys(fdir, type=t) for t in TOY_TYPES], ignore_index=True)

    if cols:
        toyfits = toyfits[cols]
        
    return toyfits


def read_toy_df_csv(csv_path, model_pars, type=None, cols=None):
    """Read toy fits from a batch CSV file, filtering by model parameters.

    Args:
        csv_path (str): Path to the batch CSV (e.g. ``toys_batch_bias.csv``).
        model_pars (dict): Model parameters to filter on
            (``co_cut``, ``ndbd_cutoff``, ``co_shift``).
            ``None`` values in the dict match NaN entries in the CSV.
        type (str | None): Filter by toy type ('fit', 'bias', 'sens').
        cols (list | None): Columns to keep. If None, all columns returned.

    Returns:
        pd.DataFrame: Toy fits with MultiIndex stat columns, same format
        as :func:`read_toy_df`.
    """
    df = read_toys_csv(csv_path, type=type)

    mask = np.ones(len(df), dtype=bool)
    for key in ('co_cut', 'ndbd_cutoff', 'co_shift'):
        val = model_pars.get(key)
        if key not in df.columns:
            continue
        if val is None:
            mask &= df[key].isna()
        elif isinstance(val, float):
            mask &= np.isclose(df[key].astype(float), val)
        else:
            mask &= df[key] == val
    df = df[mask].reset_index(drop=True)

    if cols:
        df = df[cols]

    return df