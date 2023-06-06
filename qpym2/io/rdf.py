from qpym2.utils import debug # as print

def read_rdf(fname, treename, defs=[], filters=[], add_friends=[]):
    """ Read a ROOT file and return the RDataFrame. Implements the branch
    definitions and filters.

    Args:
        fname (str): path to the ROOT file.
        treename (str): name of the tree.
        defs (list of (alias, def)): list of branch definitions.
        filters (list of string): list of filters. Default is None.
        add_friends (list of (fname, treename)): list of friend trees to be added. 
            If fname is None, use the same file.
    
    Returns:
        RDataFrame: the RDataFrame object.
    """
    from ROOT import RDataFrame as RDF
    from ROOT import TChain

    if len(add_friends) > 0:
        tc = TChain(treename)
        tc.Add(fname)

        tcfs = []
        for friend in add_friends:
            tcf = TChain(friend[1])
            tcf.Add(friend[0] or fname)
            
            tcfs.append(tcf) # Do we need to keep this list?
            tc.AddFriend(tcf)

        rdf = RDF(tc)
    else:
        rdf = RDF(treename, fname)

    for alias, def_ in defs:
        rdf = rdf.Define(alias, def_)
    for f in filters:
        rdf = rdf.Filter(f)

    return rdf

def create_rdf(input_fname, input_treename, output_fname, output_treename,
               defs=[], filters=[] , transforms=[], output_cols=[], input_add_friends=[]):
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
        input_add_friends (list of (fname, treename)): list of friend trees to be added.
            while reading the input. If fname is None, use the same file.
    """
    rdf = read_rdf(input_fname, input_treename, defs, filters, add_friends=input_add_friends)
    rdf.Snapshot(output_treename, output_fname, output_cols)