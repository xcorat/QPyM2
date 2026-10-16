"""
Provide functionality to read background model related input (JAGS output and MCs)

"""

import os

import numpy as np
import pandas as pd

from qpym2.utils.logging import debug
from qpym2.utils.mathy import get_mode


def read_bmsyst_corrmat(fpath: str = None) -> pd.DataFrame:
    """Read the background systematic covariance matrix from CSV file.

    Note: despite the name, this file contains a full covariance matrix
    (diagonal != 1), not a correlation matrix.

    Args:
        fpath: Path to covariance matrix CSV file.

    Returns:
        pd.DataFrame: Covariance matrix.
    """
    if fpath is None:
        raise ValueError("corrmat_fpath is required. Add to config: cfg.infrastructure.data_files.bkg_syst_corrmat")
    syst_corr = pd.read_csv(fpath, index_col=0)
    print(f"✓ Loaded correlation matrix from {fpath}: shape={syst_corr.shape}")
    return syst_corr


def read_mkchain(fpath: str) -> pd.DataFrame:
    """Read the raw JAGS MCMC chain from the h5 output file.

    Args:
        fpath: Path to the JAGS h5 file.

    Returns:
        pd.DataFrame: Shape (n_iterations, n_sources), columns = source names.
            Use ``.sample(n=N).T`` to get N random samples ready for
            ``make_priors()``.
    """
    mkchain_df = pd.read_hdf(fpath, 'MCMC_Chains')
    mc_names = [col for col in mkchain_df if not col.startswith('expCounts')]
    return mkchain_df[mc_names].reset_index(drop=True)


def get_mcsyst_samples(chain: pd.DataFrame, corrmat_fpath: str = None, n: int = 100) -> pd.DataFrame:
    """Generate BM-stat systematic samples via multivariate-normal perturbation.

    Selects *n* random points from the JAGS chain and, for each, draws one
    sample from a multivariate normal centred at that point with the
    background-model correlation matrix as covariance.

    This is used **exclusively** by the MC systematic runner
    (``run_syst_mc``).  Normal fits and LSS systematics should use the raw
    chain directly (``read_mkchain(...).sample(n).T``).

    Args:
        chain: Raw JAGS chain DataFrame, shape (n_iterations, n_sources).
            Typically the return value of ``read_mkchain()``.
        corrmat_fpath: Path to the BM-stat correlation matrix CSV.
        n: Number of perturbed samples to generate.

    Returns:
        pd.DataFrame: Shape (n, n_sources), columns = source names.
    """
    selected = chain.sample(n=n, replace=False)
    if corrmat_fpath is None:
        print(f"Generating {n} MC-stat samples (raw chain, no covariance perturbation)")
        return selected.reset_index(drop=True)
    syst_corr = read_bmsyst_corrmat(fpath=corrmat_fpath)
    print(f"Generating {n} MC-systematic samples with covariance matrix")
    syst_samples = [
        np.random.multivariate_normal(mean=row.values, cov=syst_corr)
        for _, row in selected.iterrows()
    ]
    return pd.DataFrame(syst_samples, columns=chain.columns)


def read_jagsh5(fpath: str) -> tuple:
    """Read JAGS output h5 file and build a source summary table.

    Args:
        fpath: Path to the JAGS h5 file.

    Returns:
        tuple: (table, mkchain) where

        * ``table`` — DataFrame indexed by source name with columns
          ``{fname, mean, std, mode}`` summarising the full chain.
        * ``mkchain`` — raw MCMC chain DataFrame (n_iterations, n_sources).
    """
    mkchain = read_mkchain(fpath)
    jags_data = pd.read_hdf(fpath, 'MC_Data')

    fnames = [
        os.path.basename(jags_data[jags_data['MCName'] == mcname].iloc[0]['FileName'])
        for mcname in mkchain.columns
    ]
    means = [mkchain[cname].mean() for cname in mkchain]
    stds  = [mkchain[cname].std()  for cname in mkchain]
    modes = [get_mode(mkchain[cname]) for cname in mkchain]

    table = pd.DataFrame({'fname': fnames, 'mean': means, 'std': stds, 'mode': modes})
    table.set_index(mkchain.columns, inplace=True)
    table.index.name = 'name'

    return table, mkchain