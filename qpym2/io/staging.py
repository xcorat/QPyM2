""" input output related to staging data. """

import numpy as np
from qpym2.utils import debug, log
from qpym2.io.rdf import create_rdf
from qpym2 import hists

def create_rdf_m2mc(fname, mcpath, out_fname='', read_config={}, write_config={}, filters=[]):
    """ Create the RDataFrame from the MC files for M2 analysis.

    Args:
        fname (str): mc file name (from the mctable, e.g. '0vbb.root')
        mcpath (str): path to the MC files.
        out_fname (str): path to the output file. If empty, use the default
            as `write_config['outpath']`/`fname`.
        read_config (dict): configuration for the read_rdf function.
            Default is {
                'treename': 'outTree',
                'evar': 'Energy',
                'esum_var': 'TotalEnergy',
                'ch1': 'Channel',
                'ch2': 'Multiplet[1]'
            }.
        write_config (dict): configuration for the create_rdf function.
            Default is {
                'outpath': 'output',
                'treename': 'uvtree',
                'out_cols': ['u', 'v', 'ch1', 'ch2']
            }.
        filters (list of string): list of filters. `Multiplicity==2` added by default.

    Returns:
        RDataFrame: the RDataFrame object.

    TODO:
        - see if we need ` && cutsTree.Included && cutsTree.FullMultiplet`
    """
    _READ_CONFIG_DEFAULT = {
        'treename': 'outTree',
        'evar': 'Energy',
        'esum_var': 'TotalEnergy',
        'ch1': 'Channel',
        'ch2': 'Multiplet[1]',
        'add_friends': [], # list of (fname, treename) of friend trees to be added
        'ds': 'Dataset',
        't': 'Time'
    }
    _FILTERS_DEFAULT = ['Multiplicity==2']
    _WRITE_CONFIG_DEFAULT = {
        'outpath': 'output',
        'treename': 'uvtree',
        'out_cols': ['u', 'v', 'ch1', 'ch2']
    }


    read_config = {**_READ_CONFIG_DEFAULT, **read_config} # merge the two dicts with read_config overriding the default
    defs = [('e1', read_config['evar']),
            ('esum', read_config['esum_var']),
            ('ch1', read_config['ch1']),
            ('ch2', read_config['ch2']),
            ('ds', read_config['ds']),
            ('t', read_config['t']),]
    # using u0, v0 as aliases for the unchanged u,v variables.
    defs.append(('u0', f'{read_config["esum_var"]}/sqrt(2)'))
    defs.append(('v0', f'(2*{read_config["evar"]} - {read_config["esum_var"]})/sqrt(2)'))
    
    input_add_friends = read_config.get('add_friends', None)

    filters = _FILTERS_DEFAULT + filters
    write_config = {**_WRITE_CONFIG_DEFAULT, **write_config} 
    if out_fname == '':
        out_fname = f'{write_config["outpath"]}/{fname}'
    
    debug('in: ', f'{mcpath}/{fname}')
    create_rdf(input_fname=f'{mcpath}/{fname}', input_treename=read_config['treename'],
               output_fname=out_fname, output_treename=write_config['treename'],
               defs=defs, filters=filters, output_cols=write_config['out_cols'],
               input_add_friends=input_add_friends)
    debug('out: ', out_fname)

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

    # print(len(defs))
    # for def_ in defs:
    #     print(def_)
    #     print(len(def_))
    #     # for d in def_:
    #     #     print(d)
    # return None
    debug('reading :', fpath)
    rdf = RDF(treename, fpath)

    if rdf.Count().GetValue() == 0:
        log(2, 'empty rdf: ', fpath)
        return hists.get_empty_hist(hm)

    for def_ in defs:
        debug('defining: ', def_[0])
        rdf = rdf.Define(def_[0], def_[1])

    # add secondary values
    rdf = rdf.Define('e1', '(u+v)/sqrt(2)')\
        .Define('e2', '(u-v)/sqrt(2)')\
        .Define('esum', 'u*sqrt(2)')\
        .Define('ediff', 'v*sqrt(2)')
    
    for f in filters:
        debug('filtering: ', f)
        rdf = rdf.Filter(f)

    return hists.get_hist_from_rdf(rdf, hm, rtype=rtype)

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
        debug('filtering: ', f)
        rdf = rdf.Filter(f)

    if rtype == 'numpy':
        # TODO: this is not efficient, we should be able to get the data directly from the rdf
        npdata = rdf.AsNumpy(columns=['u', 'v'])
        data = [ [npdata['u'][i], npdata['v'][i]] for i in range(len(npdata['u']))]
        return np.array(data, dtype=np.float64)
    else:
        raise NotImplementedError
    