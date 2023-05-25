"""
Provide functionality to read background model related input (JAGS output and MCs)

"""

import os

import pandas as pd

from qpym2.utils import get_mode
def read_mkchain(fpath: str) -> pd.DataFrame:
    """ Read the background model fit markov chain from the JAGS h5 output (`fpath`).
    Returns -- markov chain as a table with columns indicating the background component.
    """
    mkchain_df = pd.read_hdf(fpath, 'MCMC_Chains')
    mc_names = [col for col in mkchain_df if not col.startswith('expCounts')]
    mkchain_df = mkchain_df[mc_names]

    return mkchain_df

def read_jagsh5(fpath:str) -> tuple:
    """ Read JAGS output hdf file from `fpath` and create a table.

    Args:
        mctable: A pd.DataFrame with:
            index -- is the background model component name,
            fname -- basename of the MC file path
            mean -- mean of the full mkchain parameter
            std -- 
            median --
        
        mkchain: the jags markov chain

    Returns:
        tuple: (mctable, mkchain)
    """

    mkchain = read_mkchain(fpath)

    jags_data = pd.read_hdf(fpath, 'MC_Data')

    # TODO: Fix - maybe a nicer way to code this?
    # jags_data.set_index('MCName', inplace=True)
    # fnames = [os.path.basename(jags_data.loc[cname]['FileName']) for cname in mkchain]
    # fname actually isn't used by the mc_staged files anyway
    
    fnames = [ os.path.basename(jags_data[jags_data['MCName'] == mcname].iloc[0]['FileName'] )
                    for mcname in mkchain.columns]

    means = [mkchain[cname].mean() for cname in mkchain]
    stds = [mkchain[cname].std() for cname in mkchain]
    modes = [get_mode(mkchain[cname]) for cname in mkchain]

    table = pd.DataFrame({'fname': fnames, 'mean': means, 'std': stds, 'mode': modes})
    table.set_index(mkchain.columns, inplace=True)
    table.index.name = 'name'

    return table, mkchain