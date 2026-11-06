"""qls_scaling.py

This module contains the input/output functions needed for CUORE
lineshape and lineshape scaling root files.

"""

import glob
import pandas as pd
import numpy as np

from ROOT import TFile
from ROOT import RDataFrame as RDF


def read_scaling_pars(
    fnames,
    fnQ_name='fit_function_Q',
    covQ_name='covariance_matrix_Q',
    fn_sigma_name='fit_function_sigma',
    cov_sigma_name='covariance_matrix_sigma',
):
    """Read residual and resolution scaling fit out put from the given list of files

    Parameters
    ----------
    fnames : list of str
        List of root file names to read. Each file should contain the fit functions
        and covariance matrices for the residual and resolution scaling.
    fnQ_name : str (default: 'fit_function_Q')
        Name of the residual scaling fit function
    covQ_name : str (default: 'covariance_matrix_Q')
        Name of the residual scaling covariance matrix
    fn_sigma_name : str (default: 'fit_function_sigma')
        Name of the resolution scaling fit function
    cov_sigma_name : str (default: 'covariance_matrix_sigma')
        Name of the resolution scaling covariance matrix

    Returns
    -------
    fnQ : list of ROOT.TF1
        List of residual scaling fit functions
    covQ : list of ROOT.TMatrixDSym
        List of residual scaling covariance matrices
    fn_sigma : list of ROOT.TF1
        List of resolution scaling fit functions
    cov_sigma : list of ROOT.TMatrixDSym
        List of resolution scaling covariance matrices
    """
    fnQ = []
    covQ = []
    fn_sigma = []
    cov_sigma = []
    for fname in fnames:
        f = TFile(fname)
        fnQ.append(f.Get(fnQ_name).Clone())
        covQ.append(f.Get(covQ_name).Clone())
        fn_sigma.append(f.Get(fn_sigma_name).Clone())
        cov_sigma.append(f.Get(cov_sigma_name).Clone())
        f.Close()

    return fnQ, covQ, fn_sigma, cov_sigma


def read_lspars(fname, tree_name="fitParamTree", rename_cols={}):
    """Read the lineshape parameters from the given root file. Usually this contains
    the lineshape parameters for each dataset/channel, and the exposures.

    Parameters
    ----------
    fname : str
    tree_name : str (default: 'fitParamTree')
        Name of the tree containing the lineshape parameters
    rename_cols : dict (default: {})
        Dictionary mapping the names of the columns in the root file to the names
        of the columns in the returned xarray dataset

    Returns
    -------
    xarray.Dataset
        Dataset containing the lineshape parameters. The default column names are:
        'ch', 'ds', 'exposure', 'q0', 'sigma', 'pl_ratio', 'pr_ratio', 'el_ratio', 'er_ratio'
    """

    _NEW_COLS = {
        'Channel': 'ch',
        'Dataset': 'ds',
        'Exposure': 'exposure',
        # 'FWHM', 'fwhm',
        'Q_value': 'q0',
        'Sigma': 'sigma',
        'SubPeakEnergyRatio': 'pl_ratio',
        'SubPeakEnergyRatio2': 'pr_ratio',
        'SubPeakRatio': 'el_ratio',
        'SubPeakRatio2': 'er_ratio',
    }
    cols = {**_NEW_COLS, **rename_cols}

    rdf = RDF(tree_name, fname)
    npdf = rdf.AsNumpy(columns=list(str(col) for col in rdf.GetColumnNames()))
    # TODO: can we do this without pandas? or two different conversions?
    df = pd.DataFrame(npdf).rename(columns=cols).set_index(['ch', 'ds'])

    return df.to_xarray()


def read_peakshape_exposures(cfg, datasets=None, tree_name="fitParamTree", lss_data_path=None):
    """Read per-dataset exposures from PeakShape ROOT files configured in Hydra."""
    from qpym2.cfg.syst_config import (
        get_lss_data_path,
        build_lss_file_paths,
        DEFAULT_LSS_DATASETS,
    )

    if datasets is None:
        datasets = list(DEFAULT_LSS_DATASETS)

    dpath = lss_data_path or get_lss_data_path(cfg)
    _, peakshape_glob = build_lss_file_paths(dpath, datasets, cfg)

    files = sorted(glob.glob(peakshape_glob))
    if not files:
        raise FileNotFoundError(f"No PeakShape files found for glob: {peakshape_glob}")

    rows = []
    for fname in files:
        df = read_lspars(fname, tree_name=tree_name)
        try:
            exposures = df['exposure'].to_series()
        except Exception:
            exposures = df.to_dataframe()['exposure']
        exposure_by_ds = exposures.groupby('ds').sum()
        for ds, exposure in exposure_by_ds.items():
            rows.append((int(ds), float(exposure)))

    if not rows:
        raise ValueError(f"No exposures extracted from PeakShape files under {dpath}")

    exposures_df = (
        pd.DataFrame(rows, columns=['ds', 'exposure_kg_yr'])
        .drop_duplicates(subset='ds', keep='last')
        .set_index('ds')
        .sort_index()
    )
    return exposures_df


def sample_pars(fn, cov, skip_last_par=True, skip_last_cov=False, nsamples=1, sample_point=None):
    """Sample parameters from a function and covariance matrix"""
    # TODO: TODO
    # cov.Print()
    means = np.array([fn.GetParameter(i) for i in range(fn.GetNpar())])
    means = means[:-1] if skip_last_par else means

    n = cov.GetNrows()
    if skip_last_cov:
        n -= 1
    cv = np.array([[cov(i, j) for j in range(n)] for i in range(n)])

    if len(means) != cv.shape[0]:
        # print(means, cv)
        raise ValueError("Means and covariance matrix have different dimensions")

    if sample_point is None:
        sampled = np.random.multivariate_normal(means, cv, size=nsamples)
    elif sample_point == 0:
        sampled = means
    else:
        psigma = float(sample_point)

        # print(cv)
        # l = np.linalg.cholesky(cv)
        # errs = np.ones(shape=means.shape)*psigma
        # sampled = means + np.dot(l, errs)
        errs = np.sqrt(np.diag(cv)) * psigma
        sampled = means + errs
    return sampled


def read_ls_syst(cfg=None, lss_data_path=None):
    """Read LSS systematic parameters from ROOT files."""
    from qpym2.cfg.syst_config import (
        get_lss_data_path,
        build_lss_file_paths,
        DEFAULT_LSS_DATASETS,
    )

    def _fnpars(fn, skip_last=True):
        p = [fn.GetParameter(i) for i in range(fn.GetNpar())]
        if skip_last:
            return p[:-1]
        else:
            return p

    dpath = lss_data_path or (get_lss_data_path(cfg) if cfg is not None else None)
    datasets = list(DEFAULT_LSS_DATASETS)
    fnames_root, fname_ls = build_lss_file_paths(dpath, datasets, cfg)

    fns_Q, covs_Q, fns_sigma, covs_sigma = read_scaling_pars(fnames_root)
    dfps = read_lspars(fname_ls)

    qpars0 = [_fnpars(fn) for fn in fns_Q]
    qpars_1l = [
        sample_pars(fn, cov, skip_last_par=True, nsamples=1, sample_point=-1)
        for fn, cov in zip(fns_Q, covs_Q)
    ]
    qpars_1r = [
        sample_pars(fn, cov, skip_last_par=True, nsamples=1, sample_point=1)
        for fn, cov in zip(fns_Q, covs_Q)
    ]

    sigma_pars0 = [_fnpars(fn, skip_last=False) for fn in fns_sigma]
    sigma_pars_1l = [
        sample_pars(fn, cov, skip_last_par=False, nsamples=1, sample_point=-1)
        for fn, cov in zip(fns_sigma, covs_sigma)
    ]
    sigma_pars_1r = [
        sample_pars(fn, cov, skip_last_par=False, nsamples=1, sample_point=1)
        for fn, cov in zip(fns_sigma, covs_sigma)
    ]

    qpars_df = pd.DataFrame(
        {
            'q_pars': qpars0,
            'q_pars_sample_1l': qpars_1l,
            'q_pars_sample_1r': qpars_1r,
            'sigma_pars': sigma_pars0,
            'sigma_pars_sample_1l': sigma_pars_1l,  # useless?
            'sigma_pars_sample_1r': sigma_pars_1r,
        },
        index=datasets,
    ).rename_axis('ds')

    return qpars_df, dfps


def read_sigma0(df_lspars, n_ds=15, n_ch=988, ds_start=3601) -> np.ndarray:
    """2D array of sigma0 values, shape (n_ds, n_ch).

    Combines the flat lookup from ``fit.syst.sigma0_arr`` with the
    subsequent ``reshape(n_ds, n_ch)`` that every call site needs.
    Missing (ch, ds) entries are left as 0.
    """
    arr = np.zeros(n_ds * n_ch)
    for k in range(n_ds * n_ch):
        try:
            s0 = float(
                df_lspars.sigma.sel({'ds': ds_start + k // n_ch, 'ch': 1 + k % n_ch}, method=None)
            )
            arr[k] = s0
        except Exception:
            pass
    return arr.reshape(n_ds, n_ch)


# ---------------------------------------------------------------------------
# NEW: Wrapper to get lss_params dict (added 2026-03-29)
# ---------------------------------------------------------------------------


def load_lss_params(cfg=None, n_ds=15, n_ch=988, lss_data_path=None):
    """Load all LSS lineshape parameters from ROOT files via config.

    Wrapper around read_ls_syst that returns the dict format expected by
    make_lss_mctable_np and other LSS systematic functions.

    Args:
        cfg: Hydra configuration object.
        n_ds: Number of datasets (default: 15).
        n_ch: Number of channels (default: 988).
        lss_data_path: Optional explicit path to LSS ROOT files.

    Returns:
        dict with keys:
            sigma0_2d        (n_ds, n_ch)
            delta_p1_1r      (n_ds,)
            delta_p1_1l      (n_ds,)
            q_pars_arr       list[array] per ds
            q_pars_sample_1r list[array] per ds
            q_pars_sample_1l list[array] per ds
    """
    from qpym2.cfg.syst_config import DEFAULT_LSS_DATASETS

    qpars_df, dfps = read_ls_syst(cfg=cfg, lss_data_path=lss_data_path)

    datasets = list(DEFAULT_LSS_DATASETS)
    ds_start = datasets[0]

    # Extract lists from DataFrame
    q_pars_arr = qpars_df['q_pars'].tolist()
    q_pars_sample_1r = qpars_df['q_pars_sample_1r'].tolist()
    q_pars_sample_1l = qpars_df['q_pars_sample_1l'].tolist()

    # Compute delta_p1 from sigma parameter differences
    sigma_pars0 = qpars_df['sigma_pars'].tolist()
    sigma_pars_1r = qpars_df['sigma_pars_sample_1r'].tolist()
    sigma_pars_1l = qpars_df['sigma_pars_sample_1l'].tolist()

    delta_p1_1r = np.array([(s1 - s0)[1] for s1, s0 in zip(sigma_pars_1r, sigma_pars0)])
    delta_p1_1l = np.array([(s1 - s0)[1] for s1, s0 in zip(sigma_pars_1l, sigma_pars0)])

    # Build sigma0_2d
    sigma0_2d = read_sigma0(dfps, n_ds=n_ds, n_ch=n_ch, ds_start=ds_start)

    return {
        'sigma0_2d': sigma0_2d,
        'delta_p1_1r': delta_p1_1r,
        'delta_p1_1l': delta_p1_1l,
        'q_pars_arr': q_pars_arr,
        'q_pars_sample_1r': q_pars_sample_1r,
        'q_pars_sample_1l': q_pars_sample_1l,
    }


# NEW
def sample_dq_random(cfg=None, n_samples=100, seed=None, lss_data_path=None):
    """Draw N random Q-residual parameter sets from each dataset's covariance matrix.

    For each dataset, calls ``sample_pars(fn, cov, sample_point=None)``
    to draw from the multivariate normal defined by the fit function
    means and covariance matrix.

    Args:
        cfg: Hydra configuration object.
        n_samples: Number of random parameter draws per dataset.
        seed: Optional RNG seed for reproducibility.
        lss_data_path: Optional explicit path to LSS ROOT files.

    Returns:
        list of length n_samples, where each element is a list of arrays
        (one per dataset).  So ``result[i][ids]`` is the parameter
        vector for sample *i*, dataset index *ids*.
    """
    from qpym2.cfg.syst_config import get_lss_data_path, build_lss_file_paths, DEFAULT_LSS_DATASETS

    rng = np.random.default_rng(seed)

    dpath = lss_data_path or (get_lss_data_path(cfg) if cfg is not None else None)
    datasets = list(DEFAULT_LSS_DATASETS)
    fnames_root, _ = build_lss_file_paths(dpath, datasets, cfg)

    fns_Q, covs_Q, _, _ = read_scaling_pars(fnames_root)

    # For each dataset, draw n_samples parameter sets
    # samples_per_ds[ids] has shape (n_samples, n_pars)
    samples_per_ds = []
    for fn, cov in zip(fns_Q, covs_Q):
        means = np.array([fn.GetParameter(i) for i in range(fn.GetNpar())])[:-1]
        n = cov.GetNrows()
        cv = np.array([[cov(i, j) for j in range(n)] for i in range(n)])
        drawn = rng.multivariate_normal(means, cv, size=n_samples)
        samples_per_ds.append(drawn)

    # Reorganize: list of n_samples, each a list of n_ds arrays
    n_ds = len(datasets)
    result = []
    for si in range(n_samples):
        result.append([samples_per_ds[ids][si] for ids in range(n_ds)])

    return result
