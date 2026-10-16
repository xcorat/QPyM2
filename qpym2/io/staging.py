""" Implement input/output functions related to staging files """

import numpy as np

from qpym2.utils.logging import log, debug
from qpym2.utils import hists

def read_hist(fpath, hm, treename='uvtree', defs=[], filters=[], rtype='numpy'):
    """ Read a ROOT file and return the histogram. Implements the branch
    definitions and filters.

    Args:
        fpath (str): path to the ROOT file.
        hm (hist_model): hist model
        treename (str): name of the tree. Defaults to 'uvtree'.
        defs (list of (alias, def)): list of branch definitions.
        filters (list of string): list of filters. Default is None.
        rtype (str): return type. 'numpy' or 'root'. Default is 'numpy'.

    Returns:
        histogram: (histo, xedges, yedges) if return_numpy is True. otherwise return a THist variant (ROOT)
    """
    from ROOT import RDataFrame as RDF
    import os.path

    # Check if file exists. ROOT error is ugly, so skip it.
    if not os.path.isfile(fpath):
        debug('[INFO] file not found: ', fpath)
        retnp = True if rtype == 'numpy' else False 
        return hists.get_empty_hist(hm, retnp)
    
    debug('reading :', fpath)
    # TODO: combine or remove if we change the syst staging code
    try:
        rdf = RDF(treename, fpath)

        if rdf.Count().GetValue() == 0:
            log(2, 'empty rdf: ', fpath)
            # TODO: fo we need this type change?
            retnp = True if rtype == 'numpy' else False 
            return hists.get_empty_hist(hm, retnp)

        for def_ in defs:
            debug('defining: ', def_[0], def_[1])
            rdf = rdf.Define(def_[0], def_[1])

        # add secondary values
        # TODO: we might want to check if these are already defined
        #       and redefine as necessary (and Warn)
        rdf = rdf.Define('e1', '(u+v)/sqrt(2)')\
            .Define('e2', '(u-v)/sqrt(2)')\
            .Define('esum', 'u*sqrt(2)')\
            .Define('ediff', 'v*sqrt(2)')
        
        for f in filters:
            if f:
                debug('filtering: ', f)
                rdf = rdf.Filter(f)

        return hists.get_hist_from_rdf(rdf, hm, rtype=rtype)
    except Exception as e:
        debug('[DEBUG] error reading: ', fpath)
        debug('[DEBUG] error: ', e)
        # debug('DEBUG:', e)
        return hists.get_empty_hist(hm)

def read_evlist(fpath, treename, defs=[], filters=[], rtype='numpy'):
    """ Read a ROOT file and return the event list. Implements the branch
    definitions and filters.

    Args:
        fpath (str): path to the ROOT file.
        treename (str): name of the tree.
        defs (list of (alias, def)): list of branch definitions.
        filters (list of string): list of filters. Default is None.
        rtype (str): return type. Default is 'numpy', others not implementd.

    Returns:
        event list: (u, v) numpy array
    """
    # TODO: NOTE: are we using this? why? 
    from ROOT import RDataFrame as RDF

    debug('reading :', fpath)
    rdf = RDF(treename, fpath)

    if rdf.Count().GetValue() == 0:
        log(2, 'empty rdf: ', fpath)
        return np.array([])

    for alias, def_ in defs:
        debug('defining: ', alias, def_)
        rdf = rdf.Define(alias, def_)
    
    rdf = rdf.Define('e1', '(u+v)/sqrt(2)')\
        .Define('e2', '(u-v)/sqrt(2)')\
        .Define('esum', 'u*sqrt(2)')\
        .Define('ediff', 'v*sqrt(2)')
    for f in filters:
        if f:
            debug('filtering: ', f)
            rdf = rdf.Filter(f)

    if rtype == 'numpy':
        # TODO: this is not efficient, we should be able to get the data directly from the rdf
        npdata = rdf.AsNumpy(columns=['u', 'v'])
        data = [ [npdata['u'][i], npdata['v'][i]] for i in range(len(npdata['u']))]
        return np.array(data, dtype=np.float64)
    else:
        raise NotImplementedError
