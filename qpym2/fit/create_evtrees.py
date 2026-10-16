"""RDataFrame-handle cache — plural-file variant.

Main API
--------
create_evtrees(staging_cfg, comp_fnames, treename='uvtree')
    -> Dict[str, RDF | None]
    Open staging ROOT files and return RDF handles keyed by component name.

create_evlist
    Re-exported from ``create_m2uv`` for convenience.
"""

import os
from typing import Dict, List, Optional, Union

import ROOT
from ROOT import RDataFrame as RDF

from qpym2.utils.logging import debug, log


TREENAME = 'uvtree'

# ---------------------------------------------------------------------------
# Internal: open one or many files as a single RDataFrame
# ---------------------------------------------------------------------------

def _open_rdf(treename: str, fpaths: List[str]) -> Optional[RDF]:
    """Open one or more ROOT files as a single RDataFrame.

    Parameters
    ----------
    treename:
        TTree name inside the files.
    fpaths:
        Absolute paths to existing ROOT files.  Multiple files are merged
        using ROOT's multi-file constructor.

    Returns
    -------
    RDataFrame, or *None* when the resulting tree contains zero entries.
    """
    if len(fpaths) == 1:
        rdf = RDF(treename, fpaths[0])
    else:
        fvec = ROOT.std.vector['string']()
        for fp in fpaths:
            fvec.push_back(fp)
        rdf = RDF(treename, fvec)

    if rdf.Count().GetValue() == 0:
        return None
    return rdf


# ---------------------------------------------------------------------------
# Public: create RDF cache (all components)
# ---------------------------------------------------------------------------

def create_evtrees(
    staging_cfg: dict,
    comp_fnames: Dict[str, Union[str, List[str]]],
    treename: str = TREENAME,
) -> Dict[str, Optional[RDF]]:
    """Open staging ROOT files as RDataFrames and return a dict keyed by component.

    Parameters
    ----------
    staging_cfg:
        Staging configuration dict; must contain ``'outpath'`` (str) — the
        directory holding the staging ROOT files.
    comp_fnames:
        ``{comp_name: fname_or_list_of_fnames}``.  Values may be a single
        filename string or a list of filenames; multiple files are merged
        into one RDataFrame via ROOT's multi-file constructor.
    treename:
        TTree name inside the ROOT files (default ``'uvtree'``).

    Returns
    -------
    Dict mapping component name → RDataFrame (or *None* when all corresponding
    files are missing / tree is empty).
    """
    staging_dir = staging_cfg['outpath']
    evtrees: Dict[str, Optional[RDF]] = {}

    for comp_name, fnames in comp_fnames.items():
        if isinstance(fnames, str):
            fnames = [fnames]

        fpaths = [os.path.join(staging_dir, fn) for fn in fnames]
        existing = [fp for fp in fpaths if os.path.isfile(fp)]

        if not existing:
            debug(f"[create_evtrees] no files found for {comp_name}: {fpaths}")
            evtrees[comp_name] = None
            continue

        try:
            rdf = _open_rdf(treename, existing)
            if rdf is None:
                debug(f"[create_evtrees] empty tree for {comp_name}")
            else:
                debug(f"[create_evtrees] opened {len(existing)} file(s) for {comp_name}")
            evtrees[comp_name] = rdf
        except Exception as exc:
            log(2, f"[create_evtrees] error opening {comp_name}: {exc}")
            evtrees[comp_name] = None

    return evtrees
