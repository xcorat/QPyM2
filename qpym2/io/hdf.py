""" Routines needed for creating fit input (spectral histograms and data) and other HDF5 files. 

## Input HDF5 files

Used to make sure all the input data needed for analysis is self contained.

+ input[_postfix].h5 -- data, and MC event dataframes as groups.
    bkg_model: bkg_model params
    ndbd: a row similar to bkg_model table
    mchists: a group with np arrays representing the histograms as datasets. Name refers to index in table
    data: data to fit to
    mkchain: bkg_model markov chain (Not necessary)

+ fit[_postfix].h5 : fit output and component histograms. (Everything needed to replicate a fit?)
+ fit[_postfix].nc : arviz inference dataor now are {u, v, Channel, Dataset}

## BM Markov chain HDF5 files

"""

import os
import h5py
import pandas as pd
import numpy as np

from qpym2.utils.logging import log

""" Define the columns output in the bkg_model table. """
_BKG_MODEL_TABLE_COLS = [ 'mean' , 'mode']
""" Name of the group containing the mchists. """
_MCHISTS_GNAME = 'mchists'
""" Name of the group containing the signal. 
TODO: Not implenetd yet.
"""
_SIGNAL_KEY_NAME = 'signal'
""" Name of the group containing the data. """
_DATA_NAME = 'data'
""" Name of the group containing the data. """
_BLINDED_DATA_NAME = 'data_blinded'
""" Name of the group containing the bkg_model. """
_BM_KEY_NAME = 'bkg_model'
""" Name of the group containing the bkg_model markov chain. """
_MKCHAIN_KEY_NAME = 'mkchain'
""" Name of the group containing the histogram bin edges. """
_HIST_EDGES_KEY = 'hist_edges'

def _write_df(df, outpath, key, overwrite=False, **kwargs):
    """ write a dataframe to a file. If overwrite is set, overwrite data.
    
    Args:
        df (pd.DataFrame): dataframe to write.
        outpath (str): path to the output file. 
        key (str): name of the group to write to.
        overwrite (bool): overwrite the file if it exists. Default is False.
        **kwargs: additional arguments to be passed to `Pandas.to_hdf`.

    TODO: implement `mode='a'` in kwargs is needed.
    """
    with pd.HDFStore(outpath, 'a') as h5f:
        if key in h5f:
            if overwrite:
                log(2, f"Recreating group '{key}'. Previous data might be lost.")
                del h5f[key]
            else:
                raise ValueError(f'Unable to create group \'{key}\' (name already exists).'
                                    'You can overwrite data by passing the argument `overwrite=True`')

        h5f.put(key, df, **kwargs)

def _write_hists(hists, outpath, group_name, names=None, overwrite=False, append_hist=False):
    """ write the list of mchists to a file. The hists are stored as is (numpy arrays) with names given. 

    Args:
        hists (list or pd.Series of np.array): list of histograms to write.
        outpath (str): path to the output file.
        group_name (str): name of the group to write to.
        names (list of str): names of the histograms. If None, use the index of the list.
        overwrite (bool): Indicate whether to overwrite all histograms. If false and the `group_name` exists,
            throw an error.
        append_hist (bool): Indicate whether to append the histograms to the existing group. If false and the
            `group_name` exists, throw an error.
        
    """
    if not names:
        if 'name' in hists:
            names  = [ hist.name for hist in hists ]
        else:
            names = hists.index

    with h5py.File(outpath, 'a') as h5f:
        if group_name in h5f:
            if overwrite:
                log(2, f"Recreating group '{group_name}'. Previous data might be lost.")
                del h5f[group_name]
            elif append_hist:
                log(2, f"Appending histograms to group '{group_name}'. Previous data might be lost.")
            else:
                raise ValueError(f'Unable to create group \'{group_name}\' (name already exists).'
                                    'You can overwrite data by passing the argument `overwrite=True`')

        for name, hist in zip(names, hists):
            h5f.create_dataset(f'{group_name}/{name}', data=hist)

def _write_edges(edges, outpath, key=_HIST_EDGES_KEY, overwrite=False):
    """ Write histogram bin edges to an HDF5 file.

    Args:
        edges (list of np.ndarray): edge arrays, one per axis. Stored as
            integer-keyed datasets (``{key}/0``, ``{key}/1``, ...).
        outpath (str): path to the output file.
        key (str): HDF5 group name. Default is ``_HIST_EDGES_KEY``.
        overwrite (bool): overwrite if the group already exists.
    """
    with h5py.File(outpath, 'a') as h5f:
        if key in h5f:
            if overwrite:
                log(2, f"Recreating group '{key}'. Previous data might be lost.")
                del h5f[key]
            else:
                raise ValueError(f"Unable to create group '{key}' (name already exists). "
                                 "Pass overwrite=True to overwrite.")
        for i, edge in enumerate(edges):
            h5f.create_dataset(f'{key}/{i}', data=edge)


def _write_data(data, outpath, dataname, overwrite=False):
    """ write data (np.array) to h5 file. If overwrite is set, overwrite data.
    
    Args:
        data (np.array): data to write.
        outpath (str): path to the output file.
        dataname (str): name of the group to write to.
        overwrite (bool): overwrite the file if it exists. Default is False.
    """

    with h5py.File(outpath, 'a') as h5f:
        if dataname in h5f:
            if overwrite:
                log(2, f"Recreating group '{dataname}'. Previous data might be lost.")
                del h5f[dataname]
            else:
                raise ValueError(f'Unable to create group \'{dataname}\' (name already exists).'
                                    'You can overwrite data by passing the argument `overwrite=True`')

        h5f.create_dataset(dataname, data=data)

def _write_signal(signal_df, outpath, signal_name, overwrite=False, **kwargs):
    """ Write the signal component to the file. 

    Args:
        signal_df (pd.DataFrame): a single row table representing the signal.
        outpath (str): path to the output file.
        signal_name (str): name of the signal component. If None, use the index of the dataframe.
        overwrite (bool): overwrite the file if it exists. Default is False.
        **kwargs: additional arguments to be passed to `Pandas.to_hdf`.
    """
    if not overwrite and os.path.exists(outpath):
        log(2, f'WARNING: The file {outpath} exists. '
            'Any new keys will be added, but existing keys will not be overwritten.'
            'To overwrite the keys, pass the argument `overwrite=True`')
    if not signal_name:
        # TODO: check if this is a pd.DataFrame. if not use _SIGNAL_KEY_NAME?
        signal_name = signal_df.index[0]
    _write_df(signal_df[_BKG_MODEL_TABLE_COLS], outpath, key=signal_name, overwrite=overwrite, **kwargs)
    # We are adding signal to the existing group, so do not overwrite here and turn on append.
    _write_hists(signal_df['mchist'], outpath, group_name=_MCHISTS_GNAME, names=[signal_name],
                 overwrite=False, append_hist=True, **kwargs)

def write_input(outpath, mctable, signal_df, data, blinded_data=None,
                signal_name=None, mkchain=None, edges=None, overwrite=False, **kwargs):
    """ Write the input (data, hists, bm table) needed for fitting to a file.
    
    Args:
        outpath (str): path to the output file.
        mctable (pd.DataFrame): bkg_model table.
        signal_df (pd.DataFrame): a single row table representing the signal.
        data (np.array): data to fit to.
        blinded_data (np.array): blinded data to fit to. Default is None.
        signal_name (str): name of the signal component. Default is None.
        mkchain (pd.DataFrame): bkg_model markov chain. Default is None.
        overwrite (bool): overwrite the file if it exists. Default is False.
        edges (list of np.ndarray): bin edge arrays, one per histogram axis.
            Stored once in ``hist_edges/``. Default is None.
        **kwargs: additional arguments to be passed to `Pandas.to_hdf`.

    TODO: NOTE: seperate functions to write table, data, signal, etc might be better 
        if we are trying to improve performance using multiprocessing. By seperating
        the write functions, we don't need to wait till all the reading is done.

    TODO: We probably want custom 
    """
    if not overwrite and os.path.exists(outpath):
        log(2, f'WARNING: The file {outpath} exists. '
            'Any new keys will be added, but existing keys will not be overwritten.'
            'To overwrite the keys, pass the argument `overwrite=True`')
        
    _write_df(mctable[_BKG_MODEL_TABLE_COLS], outpath, key=_BM_KEY_NAME, overwrite=overwrite, **kwargs)
    _write_hists(mctable['mchist'], outpath, group_name=_MCHISTS_GNAME, overwrite=overwrite, **kwargs)

    _write_signal(signal_df, outpath, signal_name, overwrite=overwrite)
    # if not signal_name:
    #     # TODO: check if this is a pd.DataFrame. if not use _SIGNAL_KEY_NAME?
    #     signal_name = signal_df.index[0]
    # _write_df(signal_df[_BKG_MODEL_TABLE_COLS], outpath, key=signal_name, overwrite=overwrite, **kwargs)
    # # We are adding signal to the existing group, so do not overwrite here and turn on append.
    # _write_hists(signal_df['mchist'], outpath, group_name=_MCHISTS_GNAME, names=[signal_name],
    #              overwrite=False, append_hist=True, **kwargs)

    if data is not None:
        _write_data(data, outpath, dataname=_DATA_NAME, overwrite=overwrite, **kwargs)

    if blinded_data is not None:
        _write_data(blinded_data, outpath, dataname=_BLINDED_DATA_NAME,
                    overwrite=overwrite, **kwargs)

    if mkchain is not None:
        _write_df(mkchain, outpath, key=_MKCHAIN_KEY_NAME, overwrite=overwrite, **kwargs)

    if edges is not None:
        _write_edges(edges, outpath, overwrite=overwrite)
    
# TODO: redo these with proper docs
def read_mkchain(h5path):
    return pd.read_hdf(h5path, key='mkchain')

def read_bkg_model_fit(h5path):
    return pd.read_hdf(h5path, key='bkg_model')

def read_hist_edges(h5path, key=_HIST_EDGES_KEY):
    """ Read histogram bin edges from an HDF5 file.

    Args:
        h5path (str): path to the HDF5 file.
        key (str): HDF5 group name. Default is ``_HIST_EDGES_KEY``.

    Returns:
        list of np.ndarray: edge arrays sorted by integer key (axis 0, 1, ...).
    """
    with h5py.File(h5path, 'r') as h5f:
        group = h5f[key]
        return [group[str(i)][:] for i in sorted(int(k) for k in group.keys())]


def read_bkg_model(h5path):
    """ Read the background model table from input.

    Returns:
        tuple: (mctable, edges) where mctable is a pd.DataFrame with columns
            _BKG_MODEL_TABLE_COLS_ + 'mchist', and edges is a list of np.ndarray
            bin edge arrays (one per histogram axis).
    """
    mcfit_table = read_bkg_model_fit(h5path)
    mchists = read_mchists(h5path, _MCHISTS_GNAME)
    mchists.name = 'mchist'
    mctable = pd.concat([mcfit_table, mchists], axis=1, join='inner')
    try:
        edges = read_hist_edges(h5path)
    except KeyError:
        edges = None
    return mctable, edges

def read_data(h5path, dataname='data'):
    with h5py.File(h5path, 'r') as h5f:
        # TODO: Implement error handling to improve error messages
        data = np.array(h5f[dataname])
    
    return data

def read_signal(h5path, name):
    """ read the table row related to the signal component named `name`"""
    # TODO: we can do this more efficiently
    sig_table = pd.read_hdf(h5path, key=name)
    hist = read_single_mchist(h5path, mcname=name, group_name=_MCHISTS_GNAME)
    sig_table['mchist'] = [ hist ]
    sig_table['integral'] = [ np.sum(hist) ]
    return sig_table


def read_ndbd(h5path):
    return read_signal(h5path, name='ndbd')

def read_single_mchist(h5path, mcname, group_name=_MCHISTS_GNAME):
    """ read the list of mchists from a file. The hists are returned as a pd series indexed 
    by the name in the dataset
    """
    with h5py.File(h5path, 'a') as h5f:
        mchist_group = h5f[group_name]
        return mchist_group[mcname][:]

def read_mchists(h5path, group_name=_MCHISTS_GNAME):
    """ read the list of mchists from a file. The hists are returned as a pd series indexed 
    by the name in the dataset
    """
    with h5py.File(h5path, 'a') as h5f:
        mchist_group = h5f[group_name]
        hists = {name: mchist_group[name][:] for name in mchist_group.keys()}
        hists_series = pd.Series(data=hists)

        return  hists_series
         