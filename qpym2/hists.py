import math
from enum import Enum

import numpy as np
from ROOT import TH1D, TH2D

from qpym2.utils import log, debug

HistType = Enum('HistType', ['m2sum', 'm2diff', 'm2e2', 'h2ee', 'h2esume2', 'h2esumediff', 'h2uv'])

def get_h2proj(h2, axis=0):
    """ Return the projection of a 2D histogram along the specified axis 
    
    Parameters
    ----------
    h2 : numpy.ndarray
        The 2D histogram.
    axis : int
        The axis along which to project. Default is 0 (x-axis).

    Returns
    -------
    h1 : numpy.ndarray
        The 1D histogram.
    """
    return np.sum(h2, axis)

def get_sum_hist(comps, norms):
    """ Return the sum histogram from the components and their norms

    Parameters
    ----------
    comps : pandas.DataFrame
        The components table.
    norms : list
        The list of norms.

    Returns
    -------
    sum_hist : numpy.ndarray
        The sum histogram.

    """
    return np.sum(norms * comps['pdf_hist'], axis=0)

def get_hist_settings(hm):
    """ Return the nbins and range parameters for histogram according to the model

    Args:
        hm (hist_model): hist model

    Raises:
        TypeError: if the hist_type is not supported

    Returns:
        nbins, range: numbr of bins per axis and range (min, max) per axis

    TODO: clean w/ copilot
    """
    ht = hm.hist_type.value
    debug(ht, HistType.m2sum.value)
    if ht == HistType.m2sum.value:
        nbins = math.ceil((hm.esum_max - hm.esum_min)/hm.binsize)
        range = [hm.esum_min, hm.esum_max]
    elif ht == HistType.m2diff.value:
        nbins = math.ceil((hm.ediff_max - hm.ediff_min)/hm.binsize)
        range = [hm.ediff_min, hm.ediff_max]
    elif ht == HistType.m2e2.value:
        nbins = math.ceil((hm.e2max - hm.e2min)/hm.binsize)
        range = [hm.e2min, hm.e2max]
    elif ht == HistType.h2ee.value:
        nbins = [ math.ceil((hm.e1max - hm.e1min)/hm.binsize),
                  math.ceil((hm.e2max - hm.e2min)/hm.binsize) ]
        range = [[hm.e1min, hm.e1max], [hm.e2min, hm.e2max]]
    elif ht == HistType.h2esume2.value:
        nbins = [ math.ceil((hm.esum_max - hm.esum_min)/hm.binsize),
                 math.ceil((hm.e2max - hm.e2min)/hm.binsize) ]
        range = [[hm.esum_min, hm.esum_max], [hm.e2min, hm.e2max]]
    elif ht == HistType.h2esumediff.value:
        nbins = [math.ceil((hm.esum_max - hm.esum_min)/hm.binsize),
                 math.ceil((hm.ediff_max - hm.ediff_min)/hm.binsize) ]
        range = [[hm.esum_min, hm.esum_max], [hm.ediff_min, hm.ediff_max]]
    elif ht == HistType.h2uv.value:
        nbins = [ math.ceil((hm.umax - hm.umin)/hm.binsize),
                  math.ceil((hm.vmax - hm.vmin)/hm.binsize) ]
        range = [[hm.umin, hm.umax], [hm.vmin, hm.vmax]]
    else:
        raise TypeError("invalid hist_type in hm.")

    return nbins, range

def get_empty_hist(hm, return_numpy=True):
    """ Return an empty histogram according to hm

    TODO: clean with copilot
    """
    nbins, range = get_hist_settings(hm)
    if len(nbins) == 1:
        if return_numpy:
            hist = np.histogram([], bins=nbins, range=range)
        else:
             hist= TH1D(hm.name, hm.name, nbins, range[0], range[1])

    elif len(nbins) == 2:
        if return_numpy:
            hist = np.histogram2d([], [], bins=nbins, range=range)
        else:
             hist= TH2D(hm.name, hm.name, nbins[0], range[0][0], range[0][1],
                                          nbins[1], range[1][0], range[1][1])

    else:
        raise TypeError("nbins shape not supported. 'nbins` Must be 1D or 2D")

    return hist

def smooth_nph2(h2, smooth=(1, 'k5b')):
    log(-5, f'smoothing...: {smooth}')
    nx, ny = h2.shape
    h2root = TH2D('_QCOMP_h2root', '', nx, 0, nx, ny, 0, ny)
    for i in range(nx):
        for j in range(ny):
            h2root.SetBinContent(i+1,j+1, h2[i][j])
    h2root.Smooth(smooth[0], smooth[1])

    mchist_smoothed = np.array([
        [ h2root.GetBinContent(i+1, j+1) for j in range(ny) ] for i in range(nx) 
    ])

    return mchist_smoothed

# TODO: replace nbins, range calculation with the function above <15mins (after testing)
def create_hist_m2sum(rdf, hm, rtype='numpy'):
    """ create a 1D histogram of `esum` (`m2sum`) """
    nbins = math.ceil((hm.esum_max - hm.esum_min)/hm.binsize)
    if rtype == 'numpy':
        range = [hm.esum_min, hm.esum_max]
        h1sum = np.histogram(rdf.AsNumpy(['esum']), bins=nbins, range=range)
    elif rtype == 'root':
        h1sum = rdf.Histo1D( (hm.name, hm.name, nbins, hm.esum_min, hm.esum_max), 'esum').GetValue()
    else:
        raise TypeError("rtype must be 'numpy' or 'root'")
    
    return h1sum

def create_hist_m2diff(rdf, hm,  rtype='numpy'):
    """ create a 1D histogram of `ediss` (`v*\sqrt(2)`) """
    nbins = math.ceil((hm.ediff_max - hm.ediff_min)/hm.binsize)
    if rtype == 'numpy':
        range = [hm.ediff_min, hm.ediff_max]
        h1diff = np.histogram(rdf.AsNumpy(['ediff']), bins=nbins, range=range)
    elif rtype == 'root':
        h1diff = rdf.Histo1D( (hm.name, hm.name, nbins, hm.ediff_min, hm.ediff_max), 'ediff').GetValue()
    else:
        raise TypeError("rtype must be 'numpy' or 'root'")
    
    return h1diff

def create_hist_m2e2(rdf, hm, rtype='numpy'):
    """ create a 1D histogram of `e2` """
    nbins = math.ceil((hm.e2max - hm.e2min)/hm.binsize)
    if rtype == 'numpy':
        range = [hm.e2min, hm.e2max]
        h1 = np.histogram(rdf.AsNumpy(['e2']), bins=nbins, range=range)
    elif rtype == 'root':
        h1 = rdf.Histo1D( (hm.name, hm.name, nbins, hm.e2min, hm.e2max), 'e2').GetValue()
    else:
        raise TypeError("rtype must be 'numpy' or 'root'")
    return h1

def create_hist_h2ee(rdf, hm,  rtype='numpy'):
    """ create a  2D histogram of e1 vs e2 """
    nbins2 = math.ceil((hm.e2max - hm.e2min)/hm.binsize)
    nbins1 = math.ceil((hm.e1max - hm.e1min)/hm.binsize)
    if rtype == 'numpy':
        events = rdf.AsNumpy(['e1', 'e2'])
        range = [[hm.e1min, hm.e1max], [hm.e2min, hm.e2max]]
        h2 = np.histogram2d(events['e1'], events['e2'], bins=[nbins1, nbins2], range=range)
    elif rtype == 'root':
        h2 = rdf.Histo2D( (hm.name, hm.name, nbins1, hm.e1min, hm.e1max, nbins2, hm.e2min, hm.e2max), 'e1', 'e2').GetValue()
    else:
        raise TypeError("rtype must be 'numpy' or 'root'")

    return h2

def create_hist_h2esumediff(rdf, hm,  rtype='numpy'):
    """ create 2D histogram of ediff vs esum (scaled u,v) """
    nbins2 = math.ceil((hm.ediff_max - hm.ediff_min)/hm.binsize)
    nbins1 = math.ceil((hm.esum_max - hm.esum_min)/hm.binsize)
    if rtype == 'numpy':
        events = rdf.AsNumpy(['esum', 'ediff'])
        range = [[hm.esum_min, hm.esum_max], [hm.ediff_min, hm.ediff_max]]
        h2 = np.histogram2d(events['esum'], events['ediff'], bins=[nbins1, nbins2], range=range)
    elif rtype == 'root':
        h2 = rdf.Histo2D( (hm.name, hm.name, nbins1, hm.esum_min, hm.esum_max,
                                            nbins2, hm.ediff_min, hm.ediff_max), 'esum', 'ediff').GetValue()
    else:
        raise TypeError("rtype must be 'numpy' or 'root'")
    
    return h2

def create_hist_h2esume2(rdf, hm, rtype='numpy'):
    """ create 2D histogram of e2 vs esum """
    nbins1 = math.ceil((hm.esum_max - hm.esum_min)/hm.binsize)
    nbins2 = math.ceil((hm.e2max - hm.e2min)/hm.binsize)
    if rtype == 'numpy':
        events = rdf.AsNumpy(['esum', 'e2'])
        range = [[hm.esum_min, hm.esum_max], [hm.e2min, hm.e2max]]
        h2 = np.histogram2d(events['u'], events['v'], bins=[nbins1, nbins2], range=range)
    elif rtype == 'root':
        h2 = rdf.Histo2D( (hm.name, hm.name, nbins1, hm.esum_min, hm.esum_max,
                                         nbins2, hm.e2min, hm.e2max), 'esum', 'e2').GetValue()
    else:
        raise TypeError("rtype must be 'numpy' or 'root'")

    return h2

def create_hist_h2uv(rdf, hm,  rtype='numpy'):
    """ create 2D histogram of u vs v """
    nbins1 = math.ceil((hm.umax - hm.umin)/hm.binsize)
    nbins2 = math.ceil((hm.vmax - hm.vmin)/hm.binsize)

    if rtype == 'numpy':
        events = rdf.AsNumpy(['u', 'v'])
        range = [[hm.umin, hm.umax], [hm.vmin, hm.vmax]]
        h2 = np.histogram2d(events['u'], events['v'], bins=[nbins1, nbins2], range=range)
    elif rtype == 'root':
        h2 = rdf.Histo2D( (hm.name, hm.name, nbins1, hm.umin, hm.umax,
                                             nbins2, hm.vmin, hm.vmax), 'u', 'v').GetValue()
    else:
        raise TypeError("rtype must be 'numpy' or 'root'")
    
    return h2

def get_hist_from_rdf(rdf, hm, rtype='numpy'):
    """ Return a histogram according to hm. Assumes 'u' and 'v' branches are present in rdf.
    
    Args:
        rdf (RDataFrame): RDataFrame object
        hm (hist_model): hist model
        rtype (str): return type. 'numpy' or 'root'. Default is 'numpy'.

    Returns:
        histogram: (histo, xedges, yedges) if return_numpy is True. otherwise return a THist variant (ROOT)

    TODO: there must be a better version than if-else with enum
    """
    if hm.hist_type == HistType.m2sum:
        return create_hist_m2sum(rdf, hm, rtype)
    elif hm.hist_type == HistType.m2diff:
        return create_hist_m2diff(rdf, hm, rtype)
    elif hm.hist_type == HistType.m2e2:
        return create_hist_m2e2(rdf, hm, rtype)
    elif hm.hist_type == HistType.h2ee:
        return create_hist_h2ee(rdf, hm, rtype)
    elif hm.hist_type == HistType.h2esume2:
        return create_hist_h2esume2(rdf, hm, rtype)
    elif hm.hist_type == HistType.h2esumediff:
        return create_hist_h2esumediff(rdf, hm, rtype)
    elif hm.hist_type == HistType.h2uv:
        return create_hist_h2uv(rdf, hm, rtype)
    else:
        raise TypeError("invalid hist_type in hm.")
