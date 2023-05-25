""" input output related to staging data. """

from qpym2.utils import debug # as print

def read_rdf(fname, treename, defs=[], filters=[]):
    """ Read a ROOT file and return the RDataFrame. Implements the branch
    definitions and filters.

    Args:
        fname (str): path to the ROOT file.
        treename (str): name of the tree.
        defs (list of (alias, def)): list of branch definitions.
        filters (list of string): list of filters. Default is None.
    
    Returns:
        RDataFrame: the RDataFrame object.
    """
    from ROOT import RDataFrame as RDF

    rdf = RDF(treename, fname)
    for alias, def_ in defs:
        rdf = rdf.Define(alias, def_)
    for f in filters:
        rdf = rdf.Filter(f)
    return rdf

def create_rdf(input_fname, input_treename, output_fname, output_treename,
               defs=[], filters=[] , transforms=[], output_cols=[]):
    """ Create a ROOT file with the RDataFrame. Implements the branch
    definitions and filters.

    Args:
        input_fname (str): path to the ROOT file.
        input_treename (str): name of the tree.
        output_fname (str): path to the output ROOT file.
        output_treename (str): name of the output tree.
        defs (list of (alias, def)): list of branch definitions.
        filters (list of string): list of filters. 
        output_cols (list of string): list of columns to be saved in the output
            tree. 
    """
    rdf = read_rdf(input_fname, input_treename, defs, filters)
    rdf.Snapshot(output_treename, output_fname, output_cols)

def create_rdf_m2mc(fname, mcpath, read_config={}, write_config={}, filters=[]):
    """ Create the RDataFrame from the MC files for M2 analysis.

    Args:
        fname (str): mc file name (from the mctable, e.g. '0vbb.root')
        mcpath (str): path to the MC files.
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
        'ch2': 'Multiplet[1]'   
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
            ('ch2', read_config['ch2'])]
    # using u0, v0 as aliases for the unchanged u,v variables.
    defs.append(('u0', f'{read_config["esum_var"]}/sqrt(2)'))
    defs.append(('v0', f'(2*{read_config["evar"]} - {read_config["esum_var"]})/sqrt(2)'))
    
    filters = _FILTERS_DEFAULT + filters
    write_config = {**_WRITE_CONFIG_DEFAULT, **write_config} # merge the two dicts with write_config overriding the default
    
    debug('in: ', f'{mcpath}/{fname}')
    create_rdf(input_fname=f'{mcpath}/{fname}', input_treename=read_config['treename'],
               output_fname=f'{write_config["outpath"]}/{fname}', output_treename=write_config['treename'],
               defs=defs, filters=filters, output_cols=write_config['out_cols'])
    debug('out: ', f'{write_config["outpath"]}/{fname}')