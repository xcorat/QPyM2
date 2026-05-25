from ROOT import RDataFrame as RDF
from ROOT import TChain

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

    if add_friends and len(add_friends) > 0:
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

def write_rdf(rdf, output_treename, output_fname, output_cols):
    """
    Write the RDataFrame to a ROOT file.

    Args:
        rdf (RDataFrame): the RDataFrame object.
        output_treename (str): name of the output tree.
        output_fname (str): path to the output ROOT file.
        output_cols (list of str): list of columns to be saved in the output tree.

    Returns:
        RDataFrame: the RDataFrame object after writing.
    """
    out_rdf = rdf.Snapshot(output_treename, output_fname, output_cols)
    return out_rdf