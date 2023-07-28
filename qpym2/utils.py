""" Contains utulity classes and routines. 

TODO: Move things into appropriate modules.

    + Move the logs into a seperate module and implement proper logging.
    + Math utilities should go into a seperate module. (include hist routines too)
"""

_DEBUG_ = 5

def log(level:int, *args, **kwargs):
    """ log output (print) if log `level` is less than the `_DEBUG_` level 
    
    level < 0 is skipped. Use this for specific debugging while keeping print statements

    Parameters
    ----------
    level : int
        The log level.
    *args and **kwargs : passed to print
    """ 
    if level < 0:
        # negative levels are skipped.
        return
    if _DEBUG_ >= level:
        print(*args, **kwargs)

def debug(*args, **kwargs):
    log(5, *args, **kwargs)

class TempConfig:
    """ A temporary configuration class that will have all the keyword arguments 
        implemented as data members.

        Example:
            >>> cfg = TempConfig(name='hmu', threshold=threshold)
            >>> cfg.name
            'hmu'  

        TODO: find an existing implementation for this. 
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        self._kwd_names = kwds.keys()

    def __repr__(self):
        return f'TempConfig({self._kwd_names})'

def get_ncores(print_state=False):
    """ Check the number of cores available on the node. Copied from NERSC doc.

    Returns:
        int: number of cores available on the node.
    """
    import multiprocessing as mp
    import os
    import psutil

    nthreads = psutil.cpu_count(logical=True)
    ncores = psutil.cpu_count(logical=False)
    nthreads_per_core = nthreads // ncores
    nthreads_available = len(os.sched_getaffinity(0))
    ncores_available = nthreads_available // nthreads_per_core

    assert nthreads == os.cpu_count()
    assert nthreads == mp.cpu_count()

    if print_state:
        print(f'{nthreads=}')
        print(f'{ncores=}')
        print(f'{nthreads_per_core=}')
        print(f'{nthreads_available=}')
        print(f'{ncores_available=}')

    return ncores_available

from scipy.stats import gaussian_kde
from scipy.optimize import minimize_scalar
import numpy as np
def get_mode(arr, bw_method=None):
    """ Get the mode of the array `arr` using a gaussian kernel density estimation. 
    
    Parameters
    ----------
    arr : array-like
        The array to find the mode of.
    bw_method : str, scalar or callable, optional
        The method used to calculate the estimator bandwidth. This can be 'scott', 'silverman',
        a scalar constant or a callable. See `scipy.stats.gaussian_kde` for more details. 
        If None, the bandwidth is estimated using the scott's rule of thumb.
            `bw` = `n`**(-1./(arr.ndim+4))
    """
    # TODO: error checking and performance!!
    arr = np.array(arr)
    mean = arr.mean()
    marr = arr#/mean - 1

    if not bw_method:
        bw_method = marr.size**(-1./(marr.ndim+4))
    kde = gaussian_kde(marr, bw_method=bw_method)
    xmax = minimize_scalar( lambda x: -kde.evaluate(x), bounds=(marr.min(), marr.max()), method='bounded').x.item()

    # TODO: check with other methods
    # print(mean, (xmax + 1)*mean)
    return xmax #(xmax + 1)*mean

def find_mode(trace, varnames=None):
    """ Find the mode of the variables in `varnames` from the `trace` object. """
    if varnames is None:
        varnames = trace.posterior.data_vars
    return [ get_mode(trace.posterior[varname].values.flatten()) for varname in varnames ]

def test_mode(arr, ax, **kwargs):
    """ Test the mode finding method. """
    
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
    from scipy.optimize import minimize_scalar
    
    mode = get_mode(arr)
    arr = np.array(arr)
    mean = arr.mean()
    marr = arr#/mean - 1

    bw_method = marr.size**(-1./(marr.ndim+4))
    kwa_ = dict(bw_method=bw_method)
    kwa_.update(kwargs)
    
    kde = gaussian_kde(marr, **kwa_)
    xmax = minimize_scalar( lambda x: -kde.evaluate(x), bounds=(marr.min(), marr.max()), method='bounded').x.item()

    ax.hist(marr, bins=50, density=True)
    ax.axvline(x=xmax, color='red')

    ls = np.linspace(marr.min(), marr.max(), 100)
    ax.plot(ls, kde.evaluate(ls))

    return ax