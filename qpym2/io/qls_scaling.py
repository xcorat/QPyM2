""" qls_scaling.py

    This module contains the input/output functions needed for CUORE 
    lineshape and lineshape scaling root files.

"""
from random import sample
import pandas as pd
import numpy as np

from ROOT import TFile, TF1
from ROOT import RDataFrame as RDF

def read_scaling_pars(fnames, fnQ_name='fit_function_Q', covQ_name='covariance_matrix_Q',
                      fn_sigma_name='fit_function_sigma', cov_sigma_name='covariance_matrix_sigma'):
    """ Read residual and resolution scaling fit out put from the given list of files 
    
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
        fnQ.append(f.Get(fnQ_name))
        covQ.append(f.Get(covQ_name))
        fn_sigma.append(f.Get(fn_sigma_name))
        cov_sigma.append(f.Get(cov_sigma_name))
        f.Close()

    return fnQ, covQ, fn_sigma, cov_sigma

def read_lspars(fname, tree_name="fitParamTree", rename_cols={}):
    """ Read the lineshape parameters from the given root file. Usually this contains
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
        'SubPeakRatio2': 'er_ratio'
    }
    cols = {**_NEW_COLS, **rename_cols}

    rdf = RDF(tree_name, fname)
    npdf = rdf.AsNumpy(columns=list(str(col) for col in rdf.GetColumnNames()))
    # TODO: can we do this without pandas? or two different conversions?
    df = pd.DataFrame(npdf).rename(columns=cols).set_index(['ch', 'ds'])

    return df.to_xarray()
    
def sample_pars(fn, cov, skip_last_par=True, skip_last_cov=False, nsamples=1, sample_point=None):
    """ Sample parameters from a function and covariance matrix """
    # TODO: TODO
    # cov.Print()
    means = np.array([ fn.GetParameter(i) for i in range(fn.GetNpar()) ])
    means = means[:-1] if skip_last_par else means

    n = cov.GetNrows()
    if skip_last_cov: n -= 1
    cv = np.array([ [ cov(i, j) for j in range(n) ] for i in range(n) ])

    if(len(means) != cv.shape[0]):
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
        errs = np.sqrt(np.diag(cv))*psigma
        sampled = means + errs
    return sampled

