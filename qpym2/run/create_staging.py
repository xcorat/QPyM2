import logging
from logging import debug

from qpym2.io.rdf import read_rdf, write_rdf


def _create_rdf_m2mc(
    fname, mcpath,
    treename, evar, esum_var, ch1, ch2, ds, t, add_friends,
    base_cuts, outpath, output_treename, output_columns,
    out_fname="", extra_filters=None,
):
    """Internal: create and write a staged RDataFrame for one MC file.

    All parameters are explicit primitives/lists — no cfg object.
    Use :func:`create_rdf_m2ndbd` or :func:`create_rdf_m2bkg` instead.

    TODO:
        - see if we need ` && cutsTree.Included && cutsTree.FullMultiplet`
    """
    defs = [
        ('e1',   evar),
        ('esum', esum_var),
        ('ch1',  ch1),
        ('ch2',  ch2),
        ('ds',   ds),
        ('t',    t),
    ]
    defs.append(('u0', f'{esum_var}/sqrt(2)'))
    defs.append(('v0', f'(2*{evar} - {esum_var})/sqrt(2)'))
    debug(defs)

    filters = [base_cuts] + (extra_filters or [])

    if out_fname == '':
        out_fname = f'{outpath}/{fname}'

    debug('in: ', f'{mcpath}/{fname}')
    debug('cuts: ', filters)

    try:
        rdf = read_rdf(
            f'{mcpath}/{fname}',
            treename,
            defs, filters,
            add_friends=add_friends,
        )
        out_rdf = write_rdf(rdf, output_treename, out_fname, output_columns)
        n_final = out_rdf.Count().GetValue()
        print(f"[staging] final events: {n_final}  -> {out_fname}")
    except Exception as e:
        debug("Error while creating the RDataFrame: {}".format(e))
    debug('out: ', out_fname)


def _cfg_kwargs(cfg):
    """Extract _create_rdf_m2mc keyword args from a Hydra cfg object."""
    return dict(
        treename=cfg.qroot.treename,
        evar=cfg.qroot.evar,
        esum_var=cfg.qroot.esum_var,
        ch1=cfg.qroot.ch1,
        ch2=cfg.qroot.ch2,
        ds=cfg.qroot.ds,
        t=cfg.qroot.t,
        add_friends=list(cfg.qroot.add_friends),
        base_cuts=cfg.cuts.base_cuts,
        outpath=cfg.staging.outpath,
        output_treename=cfg.staging.output_treename,
        output_columns=list(cfg.staging.output_columns),
    )


def create_rdf_m2ndbd(fname, cfg, out_fname="", extra_filters=None):
    """Stage one 0vbb MC file.

    Args:
        fname (str): MC file name (e.g. ``'0vbb_1M_001.root'``).
        cfg (DictConfig): Hydra config. MC path read from
            ``cfg.infrastructure.directories.ndbd_mc_dir``.
        out_fname (str): Output path override. Defaults to
            ``cfg.staging.outpath/fname``.
        extra_filters (list of str, optional): Additional filters on top of
            ``cfg.cuts.base_cuts``.
    """
    _create_rdf_m2mc(
        fname,
        mcpath=cfg.infrastructure.directories.ndbd_mc_dir,
        **_cfg_kwargs(cfg),
        out_fname=out_fname,
        extra_filters=extra_filters,
    )


def create_rdf_m2bkg(fname, cfg, out_fname="", extra_filters=None):
    """Stage one background MC file.

    Args:
        fname (str): MC file name (e.g. ``'232Th_TeO2_001.root'``).
        cfg (DictConfig): Hydra config. MC path read from
            ``cfg.infrastructure.directories.bkg_mc_dir``.
        out_fname (str): Output path override. Defaults to
            ``cfg.staging.outpath/fname``.
        extra_filters (list of str, optional): Additional filters on top of
            ``cfg.cuts.base_cuts``.
    """
    _create_rdf_m2mc(
        fname,
        mcpath=cfg.infrastructure.directories.bkg_mc_dir,
        **_cfg_kwargs(cfg),
        out_fname=out_fname,
        extra_filters=extra_filters,
    )


# --------------------------------------------------------------------------- #
# Batch entry point
# --------------------------------------------------------------------------- #

def create_rdf_m2data(cfg, out_fname="", extra_filters=None):
    """Stage the CUORE reduced data file.

    Args:
        cfg (DictConfig): Hydra config. Input path read from
            ``cfg.infrastructure.directories.reduced_data_dir`` /
            ``cfg.infrastructure.data_files.reduced_data``.
            Output name taken from
            ``cfg.infrastructure.naming_conventions.data_name``.
        out_fname (str): Output path override. Defaults to
            ``cfg.staging.outpath/<data_name>.root``.
        extra_filters (list of str, optional): Additional filters on top of
            ``cfg.cuts.base_cuts``.
    """
    fname = cfg.infrastructure.data_files.reduced_data
    data_name = cfg.infrastructure.naming_conventions.data_name
    if out_fname == '':
        out_fname = f'{cfg.staging.outpath}/{data_name}.root'
    _create_rdf_m2mc(
        fname,
        mcpath=cfg.infrastructure.directories.reduced_data_dir,
        **_cfg_kwargs(cfg),
        out_fname=out_fname,
        extra_filters=extra_filters,
    )


# --------------------------------------------------------------------------- #
# Batch entry point
# --------------------------------------------------------------------------- #

def _bkg_fnames(cfg):
    """Background MC filenames, read from the JAGS HDF5 component table."""
    from qpym2.io.bkg_model import read_jagsh5
    table, _ = read_jagsh5(cfg.infrastructure.data_files.jags_output)
    return list(table.fname)


def _ndbd_fnames(cfg):
    """Signal-MC filenames matching ``ndbd_pattern`` in ``ndbd_mc_dir``."""
    import os
    from glob import glob
    mc_dir = cfg.infrastructure.directories.ndbd_mc_dir
    pattern = cfg.infrastructure.file_patterns.ndbd_pattern
    return sorted(os.path.basename(p) for p in glob(f'{mc_dir}/{pattern}'))


def _run_bkg(fname, cfg):
    create_rdf_m2bkg(fname, cfg)


def _run_ndbd(fname, cfg):
    create_rdf_m2ndbd(fname, cfg)


def _run_data(cfg):
    create_rdf_m2data(cfg)


def build_staging_jobs(cfg, include_bkg=True, include_ndbd=True, include_data=True):
    """Return ``[(kind, fname), ...]`` covering bkg, ndbd, and data files."""
    jobs = []
    if include_bkg:
        jobs += [('bkg', f) for f in _bkg_fnames(cfg)]
    if include_ndbd:
        jobs += [('ndbd', f) for f in _ndbd_fnames(cfg)]
    if include_data:
        fname = cfg.infrastructure.data_files.reduced_data
        jobs += [('data', fname)]
    return jobs


def main():
    """Stage all configured files (bkg MC, ndbd MC, data) in parallel.

    Usage::

        python -m qpym2.run.create_staging [--environment nersc]
                                           [--analysis m2]
                                           [--ncores N]
                                           [--ndbd-only | --bkg-only | --data-only]
    """
    import argparse
    from qpym2.cfg.hydra_adapter import get_config
    from qpym2.utils.pool import JobPool
    from qpym2.utils.logging import setup_logging

    parser = argparse.ArgumentParser(
        description='Stage all configured files (bkg MC, ndbd MC, data).',
    )
    parser.add_argument('--environment', '-e', default='nersc')
    parser.add_argument('--analysis', '-a', default='m2')
    parser.add_argument('--ncores', '-n', type=int, default=None,
                        help='Pool size; default = os.cpu_count()')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--ndbd-only', action='store_true')
    group.add_argument('--bkg-only', action='store_true')
    group.add_argument('--data-only', action='store_true')
    args = parser.parse_args()

    cfg = get_config(environment=args.environment, analysis=args.analysis)
    setup_logging(cfg)

    pool = JobPool()
    if not args.ndbd_only and not args.data_only:
        pool.extend(_run_bkg, [(f, cfg) for f in _bkg_fnames(cfg)])
    if not args.bkg_only and not args.data_only:
        pool.extend(_run_ndbd, [(f, cfg) for f in _ndbd_fnames(cfg)])
    if not args.bkg_only and not args.ndbd_only:
        pool.add(_run_data, cfg)

    print(f'[staging] dispatching {len(pool)} jobs on {args.ncores or "all"} cores')
    pool.run(ncores=args.ncores)
    print('[staging] done.')


if __name__ == '__main__':
    main()
