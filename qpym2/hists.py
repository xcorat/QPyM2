import math
from enum import Enum

import numpy as np
from ROOT import TH1D, TH2D

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

    if hm.hist_type == 'm2sum':
        nbins = math.ceil((hm.esum_max - hm.esum_min)/hm.binsize)
        range = [hm.esum_min, hm.esum_max]
    elif hm.hist_type == 'm2diff':
        nbins = math.ceil((hm.ediff_max - hm.ediff_min)/hm.binsize)
        range = [hm.ediff_min, hm.ediff_max]
    elif hm.hist_type == 'm2e2':
        nbins = math.ceil((hm.e2max - hm.e2min)/hm.binsize)
        range = [hm.e2min, hm.e2max]
    elif hm.hist_type == 'h2ee':
        nbins = [ math.ceil((hm.e1max - hm.e1min)/hm.binsize),
                  math.ceil((hm.e2max - hm.e2min)/hm.binsize) ]
        range = [[hm.e1min, hm.e1max], [hm.e2min, hm.e2max]]
    elif hm.hist_type == 'h2esume2':
        nbins = [ math.ceil((hm.esum_max - hm.esum_min)/hm.binsize),
                 math.ceil((hm.e2max - hm.e2min)/hm.binsize) ]
        range = [[hm.esum_min, hm.esum_max], [hm.e2min, hm.e2max]]
    elif hm.hist_type == 'h2esumediff':
        nbins = [math.ceil((hm.esum_max - hm.esum_min)/hm.binsize),
                 math.ceil((hm.ediff_max - hm.ediff_min)/hm.binsize) ]
        range = [[hm.esum_min, hm.esum_max], [hm.ediff_min, hm.ediff_max]]
    elif hm.hist_type == 'h2uv':
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

HistType = Enum('HistType', 'm2sum m2diff m2e2 h2ee h2esume2 h2esumediff h2uv')

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
    rdf = rdf.Define('e1', '(u+v)/sqrt(2)')\
            .Define('e2', '(u-v)/sqrt(2)')\
            .Define('esum', 'u*sqrt(2)')\
            .Define('ediff', 'v*sqrt(2)')
    
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
