""" input output related to staging data. """

from qpym2.utils import debug
from qpym2.io.rdf import read_rdf, create_rdf

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
    
    input_add_friends = read_config.get('add_friends', None)

    filters = _FILTERS_DEFAULT + filters
    write_config = {**_WRITE_CONFIG_DEFAULT, **write_config} # merge the two dicts with write_config overriding the default
    if out_fname == '':
        out_fname = f'{write_config["outpath"]}/{fname}'
    
    debug('in: ', f'{mcpath}/{fname}')
    create_rdf(input_fname=f'{mcpath}/{fname}', input_treename=read_config['treename'],
               output_fname=out_fname, output_treename=write_config['treename'],
               defs=defs, filters=filters, output_cols=write_config['out_cols'],
               input_add_friends=input_add_friends)
    debug('out: ', out_fname)
