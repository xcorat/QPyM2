"""
Create staging ROOT files from raw data for M2 analysis.

Staging files are intermediate, standardized ROOT files with u,v transformed
coordinates and quality-filtered events.

Usage:
    python run_create_staging.py [--environment ENV] [--analysis ANALYSIS]
    python run_create_staging.py --ndbd-only  # Create only signal MC
"""

import argparse
import csv
import numpy as np
import pandas as pd
from math import sqrt
import ROOT

from multiprocessing import Pool

from qpym2.cfg.hydra_adapter import get_config
from qpym2.cfg.config_builders import (
    build_staging_config, build_mc_config, build_ndbd_config,
    build_data_config, get_jags_path, get_data_names
)
from qpym2.utils.hw import get_ncores
from qpym2.io.rdf import create_rdf, read_rdf
from qpym2.io.bkg_model import read_jagsh5
from qpym2.utils.mathy import mergeto_random
from qpym2.utils.logging import log, debug

STAGING_RUN_CMD = 'python qpym2/run/run_create_staging.py '


def _build_job_cuts(read_config, base_cuts):
    cuts = [base_cuts]
    if 'extra_cuts' in read_config:
        cuts.append(read_config['extra_cuts'])
    return cuts


def _create_staging_job(job):
    create_rdf_m2mc(
        fname=job['fname'],
        mcpath=job['mcpath'],
        out_fname=job['out_fname'],
        read_config=job['read_config'],
        write_config=job['write_config'],
        filters=job['filters'],
    )
    return job['kind'], job['fname']


def build_staging_jobs(jags_fpath, mc_cfg, ndbd_cfg, data_cfg, staging_cfg,
                       include_mc=True, include_ndbd=True, include_data=True):
    """Build a flat list of staging jobs with explicit input/output params."""
    jobs = []
    cuts = staging_cfg['cuts']

    if include_mc:
        table, _ = read_jagsh5(jags_fpath)
        cuts_mc = _build_job_cuts(mc_cfg, cuts)
        for fname in table.fname:
            jobs.append({
                'kind': 'mc',
                'fname': fname,
                'mcpath': mc_cfg['inpath'],
                'out_fname': '',
                'read_config': mc_cfg,
                'write_config': staging_cfg,
                'filters': cuts_mc,
            })

    if include_ndbd:
        cuts_ndbd = _build_job_cuts(ndbd_cfg, cuts)
        jobs.append({
            'kind': 'ndbd',
            'fname': ndbd_cfg['fname'],
            'mcpath': ndbd_cfg['inpath'],
            'out_fname': f'{staging_cfg["outpath"]}/{ndbd_cfg["out_fname"]}',
            'read_config': ndbd_cfg,
            'write_config': staging_cfg,
            'filters': cuts_ndbd,
        })

    if include_data:
        cuts_data = _build_job_cuts(data_cfg, cuts)
        jobs.append({
            'kind': 'data',
            'fname': data_cfg['fname'],
            'mcpath': data_cfg['inpath'],
            'out_fname': f'{staging_cfg["outpath"]}/{data_cfg["out_fname"]}',
            'read_config': data_cfg,
            'write_config': staging_cfg,
            'filters': cuts_data,
        })

    return jobs


def run_jobs_local(jobs, multiproc=True, ncores=None):
    """Run staging jobs locally (parallel for non-SLURM execution)."""
    if not jobs:
        debug('No staging jobs to run.')
        return []

    if multiproc and len(jobs) > 1:
        if ncores is None or ncores < 1:
            ncores = get_ncores()
        debug(f'starting local pool with ncores={ncores}, njobs={len(jobs)}')
        with Pool(processes=ncores) as pool:
            results = pool.map(_create_staging_job, jobs)
        debug('end pool...', len(results))
        return results

    debug('creating staging files without parallelization...')
    results = []
    for job in jobs:
        results.append(_create_staging_job(job))
    return results


def write_jobslist_tsv(jobs, jobslist_path, environment, analysis):
    """Write staging jobs as a TSV list for SLURM array submission."""
    fieldnames = ['task_id', 'job_type', 'fname', 'command']
    rows = []

    for idx, job in enumerate(jobs, start=1):
        command = (
            STAGING_RUN_CMD +
            f'--environment {environment} '
            f'--analysis {analysis} '
            f'--job-type {job["kind"]} '
            f'--job-fname {job["fname"]}'
        )
        rows.append({
            'task_id': idx,
            'job_type': job['kind'],
            'fname': job['fname'],
            'command': command,
        })

    with open(jobslist_path, 'w', newline='') as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter='\t',
            lineterminator='\n',
        )
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


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
        - **REFACTOR**: separate creating the rdf from writing atleast.
    """

    # DEFAULTS: TODO: Move to __init__?
    # Would it mess up when importing?
    _READ_CONFIG_DEFAULT = {
        'treename': 'outTree',
        'evar': 'Energy',
        'esum_var': 'TotalEnergy',
        'ch1': 'Channel',
        'ch2': 'Multiplet[1]',
        # TODO: `add_friends` is not used, so maybe we can remove it..
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
    debug(defs)
    input_add_friends = read_config.get('add_friends', None)

    filters = _FILTERS_DEFAULT + filters
    write_config = {**_WRITE_CONFIG_DEFAULT, **write_config} 
    if out_fname == '':
        out_fname = f'{write_config["outpath"]}/{fname}'
    
    debug('in: ', f'{mcpath}/{fname}')
    debug('cuts: ', filters)
    create_rdf(input_fname=f'{mcpath}/{fname}', input_treename=read_config['treename'],
               output_fname=out_fname, output_treename=write_config['treename'],
               defs=defs, filters=filters, output_cols=write_config['out_cols'],
               input_add_friends=input_add_friends)
    debug('out: ', out_fname)


def create_blinded_rdf(data_cfg, blinding_cfg, staging_cfg, write=False):
    """ Create blinded data.

    TODO: complete this, add docs
    """
    _READ_CONFIG_DEFAULT = {
        'treename': 'outTree',
        'evar': 'Energy',
        'esum_var': 'TotalEnergy',
        'ch1': 'Channel',
        'ch2': 'Multiplet[1]',
        # TODO: `add_friends` is not used, so maybe we can remove it..
        'add_friends': [], # list of (fname, treename) of friend trees to be added 
        'ds': 'Dataset',
        't': 'Time'
    }

    _BLINDING_CONFIG_DEFAULT = {
        'n_inject': 0,
    }

    # TODO: ummm... meh?
    blinding_cfg = {**_BLINDING_CONFIG_DEFAULT, **blinding_cfg}
    blinding_cfg['read_cfg'] = {**_READ_CONFIG_DEFAULT, **blinding_cfg['read_cfg']}
    ndbd_cfg = blinding_cfg['read_cfg']

    cuts = staging_cfg['cuts']
    cuts_data = [cuts]
    data_cfg = {**_READ_CONFIG_DEFAULT, **data_cfg} # merge the two dicts with read_config overriding the default

    if 'extra_cuts' in data_cfg:
        cuts_data = [cuts, data_cfg['extra_cuts']]
    _FILTERS_DEFAULT = ['Multiplicity==2']
    cuts_data = _FILTERS_DEFAULT + cuts_data

    defs = [('e1', data_cfg['evar']),
            ('esum', data_cfg['esum_var']),
            ('ch1', data_cfg['ch1']),
            ('ch2', data_cfg['ch2']),
            ('ds', data_cfg['ds']),
            ('t', data_cfg['t']),]
    # using u0, v0 as aliases for the unchanged u,v variables.
    defs.append(('u0', f'{data_cfg["esum_var"]}/sqrt(2)'))
    defs.append(('v0', f'(2*{data_cfg["evar"]} - {data_cfg["esum_var"]})/sqrt(2)'))
    debug(defs)

    add_friends_data = None
    input_fname_data = f"{data_cfg['inpath']}/{data_cfg['fname']}"
    debug('in: ', input_fname_data)
    debug('cuts: ', cuts_data)
    rdf_data = read_rdf(input_fname_data, data_cfg['treename'], defs, 
                        cuts_data, add_friends=add_friends_data)

    add_friends_0v = None
    cuts_ndbd = [cuts, ndbd_cfg['extra_cuts']]
    input_fname_0v = f'{ndbd_cfg["inpath"]}/{ndbd_cfg["fname"]}'
    debug('in: ', f'{input_fname_0v}')
    debug('cuts: ', cuts_ndbd)
    rdf_0v = read_rdf(input_fname_0v, ndbd_cfg['treename'],
                      defs, cuts_ndbd, add_friends=add_friends_0v)

    def _get_inj_frac_(n_inject, nselected_0v):
        n_inject_min, n_inject_max = int(n_inject-sqrt(n_inject)), \
                                     int(n_inject+sqrt(n_inject))
        n = np.random.uniform(n_inject_min, n_inject_max)
        return n/nselected_0v

    selected_fraction = _get_inj_frac_(n_inject=blinding_cfg['n_inject'],
                                       nselected_0v=rdf_0v.Count().GetValue())

    """ declare a ROOT function to select a fraction of events randomly """
    if not hasattr(ROOT, 'pass_random'):
        ROOT.gInterpreter.Declare( """
        bool pass_random(double fraction){
            TRandom3 rnd;   
            rnd.SetSeed(0);
            return rnd.Uniform() < fraction;
        }
        """ )

    # create a list of random events to insert
    rdf_uv0v_rnd = rdf_0v.Filter(f'pass_random({selected_fraction})')

    df_uv0v_rnd = pd.DataFrame(rdf_uv0v_rnd.AsNumpy(columns=staging_cfg['out_cols']))
    df_uvdata = pd.DataFrame(rdf_data.AsNumpy(columns=staging_cfg['out_cols']))
    df_blinded = mergeto_random(df_uvdata, df_uv0v_rnd)

    log(0, 'data+0v:', len(df_uvdata), len(df_uv0v_rnd))

    np_blinded = {key: df_blinded[key].values for key in df_blinded.columns}
    rdf_blinded = ROOT.RDF.FromNumpy(np_blinded)

    if write:

        out_cols = staging_cfg['out_cols']
        outdir = staging_cfg['outpath']
        out_fpath = f'{outdir}/{blinding_cfg["out_fname"]}'
        if out_cols:
            out_rdf = rdf_blinded.Snapshot('uvtree', out_fpath, out_cols)
        else:
            out_rdf = rdf_blinded.Snapshot('uvtree', out_fpath)
        n_final = out_rdf.Count().GetValue()
        print(f"[staging] final events: {n_final}  -> {out_fpath}")

    return rdf_blinded


def create_ndbd(ndbd_cfg, staging_cfg):
    """ Create 0νββ signal MC staging file. 
    
    Args:
        ndbd_cfg (dict): NDBD input configuration from config_builders
        staging_cfg (dict): Staging output configuration from config_builders
    """
    read_config_ndbd = ndbd_cfg

    outpath = staging_cfg['outpath']
    cuts = staging_cfg['cuts']
    
    out_fname_ndbd = f'{outpath}/{ndbd_cfg["out_fname"]}'

    fname_ndbd = read_config_ndbd['fname']

    cuts_ndbd = [cuts]
    if 'extra_cuts' in read_config_ndbd:
        cuts_ndbd = [cuts, read_config_ndbd['extra_cuts']]


    debug('creating staging ndbd file...')
    ret_0v = create_rdf_m2mc(fname_ndbd, 
                             mcpath=read_config_ndbd['inpath'],
                             read_config=read_config_ndbd,
                             write_config=staging_cfg,
                             filters=cuts_ndbd,
                             out_fname=out_fname_ndbd)


def create_single_job_from_kind(jags_fpath, mc_cfg, ndbd_cfg, data_cfg,
                                staging_cfg, job_type, job_fname):
    jobs = build_staging_jobs(
        jags_fpath, mc_cfg, ndbd_cfg, data_cfg, staging_cfg,
        include_mc=(job_type == 'mc'),
        include_ndbd=(job_type == 'ndbd'),
        include_data=(job_type == 'data'),
    )

    selected = [job for job in jobs if job['fname'] == job_fname]
    if len(selected) != 1:
        raise ValueError(
            f'Expected exactly one job for type={job_type}, fname={job_fname}; '
            f'found {len(selected)}'
        )
    _create_staging_job(selected[0])


def create_all(jags_fpath, mc_cfg, ndbd_cfg, data_cfg,
                staging_cfg, multiproc=True, ncores=None):
    """ Create MC, 0νββ and data staging files in parallel. 
    
    Args:
        jags_fpath (str): Path to JAGS HDF5 file with MC component table
        mc_cfg (dict): Background MC configuration from config_builders
        ndbd_cfg (dict): NDBD configuration from config_builders
        data_cfg (dict): Data configuration from config_builders
        staging_cfg (dict): Staging output configuration from config_builders
        multiproc (bool): if True, use multiprocessing.
        ncores (int): number of cores to use. If None, use all available cores.

    TODO: also need separate functions to run these independently
    TODO: this actually does not run in parallel, we need to replace `map` method with 
        an async method and read the output later.
    """
    jobs = build_staging_jobs(jags_fpath, mc_cfg, ndbd_cfg, data_cfg, staging_cfg)
    return run_jobs_local(jobs, multiproc=multiproc, ncores=ncores)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Create staging ROOT files from raw data for M2 analysis.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default tech_note environment
  python run_create_staging.py
  
  # Run with specific environment
  python run_create_staging.py --environment production
  
  # Create only signal MC (NDBD)
  python run_create_staging.py --ndbd-only
  
    # Use specific analysis
    python run_create_staging.py --environment tech_note --analysis tech_note
        """
    )
    
    parser.add_argument(
        '--environment', '-e',
        default='tech_note',
        help='Environment name (default: tech_note)'
    )
    
    parser.add_argument(
        '--analysis', '-a',
        default=None,
        help='Analysis name (default: same as environment)'
    )

    parser.add_argument(
        '--ndbd-only',
        action='store_true',
        help='Create only NDBD (signal MC) staging file'
    )
    
    parser.add_argument(
        '--no-multiproc',
        action='store_true',
        help='Disable multiprocessing (run sequentially)'
    )
    
    parser.add_argument(
        '--ncores', '-n',
        type=int,
        default=None,
        help='Number of cores for multiprocessing (default: all available)'
    )

    parser.add_argument(
        '--create-jobslist',
        action='store_true',
        help='Create staging jobs TSV and exit (for SLURM submission)'
    )

    parser.add_argument(
        '--jobslist',
        default='scripts/staging_params.tsv',
        help='Path to staging jobs TSV (default: scripts/staging_params.tsv)'
    )

    parser.add_argument(
        '--job-type',
        choices=['mc', 'ndbd', 'data'],
        default=None,
        help='Run a single staging job type (used by SLURM array commands)'
    )

    parser.add_argument(
        '--job-fname',
        default=None,
        help='Input ROOT filename to run for --job-type (used by SLURM array commands)'
    )
    
    return parser.parse_args()


def main():
    """ Running this as a notebook took ~22mins while running it with in parallel with 
        `multiproc=True` in shell took ~27mins  
    """
    args = parse_args()
    
    # Load Hydra configuration
    analysis = args.analysis if args.analysis else args.environment
    cfg = get_config(
        environment=args.environment,
        analysis=analysis
    )
    
    # Build configuration dicts using builders
    staging_cfg = build_staging_config(cfg)
    jags_fpath = get_jags_path(cfg)
    ndbd_cfg = build_ndbd_config(cfg)
    data_cfg = build_data_config(cfg, blinded=False)
    mc_cfg = build_mc_config(cfg)

    if args.job_type:
        if args.job_fname is None:
            raise ValueError('--job-fname is required when --job-type is provided')
        create_single_job_from_kind(
            jags_fpath=jags_fpath,
            mc_cfg=mc_cfg,
            ndbd_cfg=ndbd_cfg,
            data_cfg=data_cfg,
            staging_cfg=staging_cfg,
            job_type=args.job_type,
            job_fname=args.job_fname,
        )
        return

    if args.create_jobslist:
        jobs = build_staging_jobs(
            jags_fpath=jags_fpath,
            mc_cfg=mc_cfg,
            ndbd_cfg=ndbd_cfg,
            data_cfg=data_cfg,
            staging_cfg=staging_cfg,
            include_mc=not args.ndbd_only,
            include_ndbd=True,
            include_data=not args.ndbd_only,
        )
        n_jobs = write_jobslist_tsv(
            jobs=jobs,
            jobslist_path=args.jobslist,
            environment=args.environment,
            analysis=analysis,
        )
        log(0, f'Wrote staging jobs list: {args.jobslist} ({n_jobs} jobs)')
        return
    
    if args.ndbd_only:
        # Create only NDBD staging file
        create_ndbd(ndbd_cfg, staging_cfg)
    else:
        jobs = build_staging_jobs(jags_fpath, mc_cfg, ndbd_cfg, data_cfg, staging_cfg)
        run_jobs_local(jobs, multiproc=not args.no_multiproc, ncores=args.ncores)
    
    # Create blinded data (optional, commented out by default)
    # data_name, data_blinded_name = get_data_names(cfg)
    # blinding_cfg = {
    #     'read_cfg': ndbd_cfg,
    #     'n_inject': 175,
    #     'out_fname': f'{data_blinded_name}.root'
    # }
    # create_blinded_rdf(data_cfg, blinding_cfg, staging_cfg, write=True)


def main_ndbd_only():
    """ Create only the ndbd staging file (legacy entry point). """
    # Load Hydra config with default environment
    cfg = get_config(
        environment='tech_note',
        analysis='tech_note'
    )
    
    ndbd_cfg = build_ndbd_config(cfg)
    staging_cfg = build_staging_config(cfg)
    
    create_ndbd(ndbd_cfg, staging_cfg)


if __name__ == '__main__':
    main()
