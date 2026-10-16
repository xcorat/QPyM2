"""run_syst_lss.py — Parallel LSS systematic fits (sscl / dq).

Two independent subcommands:

  ``python -m qpym2.run.run_syst_lss sscl``
      Sigma-scaling: fixed +1σ parameter set, one fit per dataset.
      Per-event stochasticity from ``rng.normal()`` in ``sigma_scaling_np``.

  ``python -m qpym2.run.run_syst_lss dq --n-dq 100``
      Q-residual: N random parameter draws from each dataset's covariance
      matrix, each producing a deterministic polynomial transform.

Results are saved to separate files:
  ``lss_sscl_results.nc`` and ``lss_dq_results.nc``.
"""

import argparse
import gc
import multiprocessing
import os
import resource
import sys
import time

import numpy as np
import pandas as pd
import xarray as xr

from qpym2.utils.logging import log, start_timer, log_elapsed


# ---------------------------------------------------------------------------
# Pool worker (module-level; data passed via spawn initializer)
# ---------------------------------------------------------------------------

_WORKER_SHARED = {}


def _worker_pool_init(shared):
    """Initializer for spawn pool workers: receive shared data by pickle.

    Using spawn (not fork) avoids libffi MAP_SHARED corruption: libffi 3.4+
    backs closures with MAP_SHARED memfd regions that are shared between
    parent and child after fork.  When either process frees a closure the
    shared allocator state is corrupted, causing dlfree to abort in the
    other process.  spawn creates an independent Python interpreter so each
    worker imports ROOT cleanly without inheriting corrupted closures.
    """
    global _WORKER_SHARED
    _WORKER_SHARED = shared


def _lss_fit_single(job, ws=None):
    """Run one LSS systematic fit inside a Pool worker.

    Thin wrapper that unpacks shared state, pre-indexes LSS parameters,
    and delegates to :mod:`qpym2.fit.syst_lss` functions.
    """
    if ws is None:
        ws = _WORKER_SHARED

    # Per-worker PyTensor / ArviZ setup
    _worker_dir = os.path.join(
        os.environ.get('TMPDIR', '/tmp'),
        f'.pytensor_worker_{os.getpid()}',
    )
    os.makedirs(_worker_dir, exist_ok=True)
    os.environ['PYTENSOR_FLAGS'] = f'compiledir={_worker_dir}'
    import pytensor  # noqa: F401
    import pathlib
    (pathlib.Path.home() / '.cache' / 'arviz').mkdir(parents=True, exist_ok=True)

    # Unpack — dq jobs carry an extra q_pars_sample element
    if len(job) == 4:
        sname, scfg, ds, q_pars_sample = job
    else:
        sname, scfg, ds = job
        q_pars_sample = None

    import time as _time
    from qpym2.fit.syst_lss import (
        run_sscl_fit, run_dq_fit,
        extract_posterior_stats, apply_bias_correction,
    )

    label = f"{sname} ds={ds}"
    _t0 = _time.perf_counter()
    print(f"  [{label}] fitting ...", flush=True)

    # Pre-compute cuts from fit_cfg
    cuts = ws['fit_cfg']['cuts']
    cuts_mc = [f for f in [cuts.get('common'), cuts.get('mc_only')] if f]
    cuts_ndbd = [f for f in [cuts.get('common'), cuts.get('ndbd_only')] if f]

    try:
        # Pre-index lss_params for this dataset
        ids = ds - ws['ds_start']
        stype = scfg['type']

        if stype == 'sscl':
            spt = scfg.get('spt', 1)
            dp1_key = 'delta_p1_1r' if spt >= 0 else 'delta_p1_1l'
            trace = run_sscl_fit(
                ws['evdata'], ws['mcfit_table'], ws['data'],
                ws['hm'], ws['comp_cfg'], ws['nominal_priors'],
                target_ds=ds,
                ch_sigma0=ws['lss_params']['sigma0_2d'][ids],
                delta_p1=ws['lss_params'][dp1_key][ids],
                chains=ws['chains'], draws=ws['draws'],
                cuts_mc=cuts_mc, cuts_ndbd=cuts_ndbd,
            )
        elif stype == 'dq':
            if q_pars_sample is not None:
                qps = q_pars_sample[ids]
            else:
                ptn = scfg['ptn']
                qps = ws['lss_params'][f'q_pars_sample_{ptn}'][ids]
            trace = run_dq_fit(
                ws['evdata'], ws['mcfit_table'], ws['data'],
                ws['hm'], ws['comp_cfg'], ws['nominal_priors'],
                target_ds=ds,
                q_pars_mean=ws['lss_params']['q_pars_arr'][ids],
                q_pars_sample=qps,
                chains=ws['chains'], draws=ws['draws'],
                cuts_mc=cuts_mc, cuts_ndbd=cuts_ndbd,
            )
        else:
            raise ValueError(f"Unknown syst type: {stype}")

        # Save trace if requested
        if ws.get('save_trace_dir'):
            trace_path = os.path.join(
                ws['save_trace_dir'], f'trace_{sname}_ds{ds}.nc',
            )
            trace.to_netcdf(trace_path)
            print(f"  [{label}] saved trace to {trace_path}", flush=True)

        stats = extract_posterior_stats(trace)
        component_modes = stats['component_modes']
        p16, p90 = stats['p16'], stats['p90']

        mode_bc = p16_bc = p90_bc = np.nan
        if ws.get('bias_popt') is not None:
            bc = apply_bias_correction(
                stats['ndbd_samples'], ws['bias_popt'], ws['bias_pcov'],
            )
            mode_bc, p16_bc, p90_bc = bc['mode_bc'], bc['p16_bc'], bc['p90_bc']

        _elapsed = _time.perf_counter() - _t0
        print(f"  [{label}] ndbd_mode={component_modes.get('ndbd', np.nan):.4f}"
              f"  ({_elapsed:.1f}s)", flush=True)
        return (sname, ds, component_modes, p16, p90, mode_bc, p16_bc, p90_bc)

    except Exception as e:
        _elapsed = _time.perf_counter() - _t0
        print(f"  [{label}] FAILED ({_elapsed:.1f}s): {e}", flush=True)
        return (sname, ds, {}, np.nan, np.nan, np.nan, np.nan, np.nan)


# ---------------------------------------------------------------------------
# Shared setup
# ---------------------------------------------------------------------------

def _setup_common(args):
    """Load config, data, LSS params, and materialize RDFs to numpy.

    Returns a dict with everything both sscl and dq need.
    """
    from qpym2.cfg.hydra_adapter import get_config
    from qpym2.cfg.config_builders import (
        get_h2uv_model, get_jags_path, get_fit_dir,
        build_input_config, get_comp_cfg,
    )
    from qpym2.cfg.model_utils import build_model_dir, get_shift_suffix
    from qpym2.cfg.syst_config import get_systematic_dir, get_lss_data_path
    from qpym2.fit.params import model_cfg
    from qpym2.fit.create_evtrees import create_evtrees
    from qpym2.fit.create_m2uv import create_evlist
    from qpym2.io.qls_scaling import load_lss_params
    from qpym2.io.bkg_model import read_jagsh5
    from qpym2.fit.bias import read_bias_csv

    analysis = args.analysis or args.environment
    cfg = get_config(
        environment=args.environment,
        analysis=analysis,
        hyperparams=args.hyperparams,
    )

    jags_fpath = get_jags_path(cfg)
    h2uv_model = get_h2uv_model(cfg)
    comp_cfg_base = get_comp_cfg(cfg)
    input_cfg = build_input_config(cfg)
    syst_dir = get_systematic_dir(cfg)

    mcfg = model_cfg(
        h2uv_model,
        co_cut=args.co_cut,
        ndbd_cutoff=args.ndbd_cutoff,
        co_shift=args.co_shift,
    )
    fit_cfg = {'h2uv_model': mcfg['h2uv_model'], 'cuts': mcfg['cuts']}
    hm = fit_cfg['h2uv_model']

    # LSS data path
    lss_data_path = args.lss_data_path or get_lss_data_path(cfg)
    if lss_data_path is None:
        print("[setup] ERROR: --lss-data-path not provided and "
              "cfg.infrastructure.data_files.lss_data_path not set.",
              file=sys.stderr)
        sys.exit(1)

    datasets = list(range(args.ds_start, args.ds_end))
    n_ds = len(datasets)

    # LSS parameters
    print("[setup] loading LSS parameters ...", flush=True)
    lss_params = load_lss_params(cfg, n_ds=n_ds)

    # Background model (for mcfit_table metadata only — priors from nominal fit)
    print("[setup] reading background model ...", flush=True)
    mcfit_table, _ = read_jagsh5(jags_fpath)

    # Staging RDFs
    print("[setup] opening staging ROOT files ...", flush=True)
    comp_fnames = dict(zip(mcfit_table.index, mcfit_table.fname))
    comp_fnames['ndbd'] = 'ndbd.root'
    evtrees = create_evtrees(input_cfg, comp_fnames)

    # Data
    print("[setup] reading data ...", flush=True)
    data = create_evlist(mcfg, input_cfg)

    # Read priors from nominal fit (HDF5 only — skip trace to avoid
    # NetCDF C-library state that corrupts on fork)
    print("[setup] reading nominal fit priors ...", flush=True)
    fit_dir = get_fit_dir(cfg)
    model_dir = build_model_dir(
        co_cut=args.co_cut, ndbd_cutoff=args.ndbd_cutoff,
    )
    shift_suffix = get_shift_suffix(co_shift=args.co_shift)
    nominal_fit_h5 = os.path.join(fit_dir, model_dir, f'fit{shift_suffix}.h5')
    nominal_comps = pd.read_hdf(nominal_fit_h5, key='components')
    nominal_priors = list(nominal_comps['prior'])
    print(f"[setup] nominal priors ({len(nominal_priors)} components): "
          f"{[p[0] if isinstance(p, tuple) else type(p).__name__ for p in nominal_priors]}",
          flush=True)

    comp_cfg = list(comp_cfg_base)

    # Materialize RDFs to numpy
    print("[setup] materializing RDFs to numpy ...", flush=True)
    _NP_COLS = ['u0', 'v0', 'ch1', 'ch2', 'ds']
    evdata = {}
    for comp, rdf in evtrees.items():
        if rdf is None:
            evdata[comp] = None
        else:
            evdata[comp] = {c: np.array(a) for c, a
                            in rdf.AsNumpy(columns=_NP_COLS).items()}
    n_mat = sum(1 for v in evdata.values() if v is not None)
    print(f"[setup] materialized {n_mat} components", flush=True)
    del evtrees

    # Force GC to free cppyy/ctypes libffi closures created during
    # rdf.AsNumpy() before fork().  Closures surviving into child processes
    # trigger a double-free in libffi's allocator (dlfree abort) when the
    # child's GC collects them.
    gc.collect()
    gc.collect()

    # Bias correction
    bias_popt = bias_pcov = None
    if args.bias_csv:
        model_pars = {
            'co_cut': args.co_cut,
            'ndbd_cutoff': args.ndbd_cutoff,
            'co_shift': args.co_shift,
        }
        bias_popt, bias_pcov = read_bias_csv(
            args.bias_csv, model_pars, xmin=args.xmin_bias,
        )
        print(f"[setup] bias correction: a={bias_popt[0]:.4f}, "
              f"b={bias_popt[1]:.4f}", flush=True)

    return {
        'cfg': cfg,
        'datasets': datasets,
        'n_ds': n_ds,
        'lss_params': lss_params,
        'mcfit_table': mcfit_table,
        'data': data,
        'hm': hm,
        'fit_cfg': fit_cfg,
        'comp_cfg': comp_cfg,
        'nominal_priors': nominal_priors,
        'evdata': evdata,
        'syst_dir': syst_dir,
        'bias_popt': bias_popt,
        'bias_pcov': bias_pcov,
    }


def _set_worker_shared(ctx, args):
    """Populate _WORKER_SHARED from setup context."""
    global _WORKER_SHARED
    _WORKER_SHARED = {
        'evdata':          ctx['evdata'],
        'mcfit_table':     ctx['mcfit_table'],
        'fit_cfg':         ctx['fit_cfg'],
        'lss_params':      ctx['lss_params'],
        'data':            ctx['data'],
        'hm':              ctx['hm'],
        'comp_cfg':        ctx['comp_cfg'],
        'nominal_priors':  ctx['nominal_priors'],
        'chains':          args.chains,
        'draws':           args.draws,
        'ds_start':        ctx['datasets'][0],
        'bias_popt':       ctx['bias_popt'],
        'bias_pcov':       ctx['bias_pcov'],
        'save_trace_dir':  ctx.get('save_trace_dir'),
    }


def _build_xr_datavars(ndbd_mode, ndbd_p16, ndbd_p90,
                       ndbd_mode_bc, ndbd_p16_bc, ndbd_p90_bc,
                       comp_modes, dim_names, has_bias):
    """Build xarray data_vars dict from result arrays."""
    dvars = {
        'ndbd_mode': (dim_names, ndbd_mode),
        'ndbd_p16':  (dim_names, ndbd_p16),
        'ndbd_p90':  (dim_names, ndbd_p90),
    }
    if has_bias:
        dvars['ndbd_mode_bc'] = (dim_names, ndbd_mode_bc)
        dvars['ndbd_p16_bc']  = (dim_names, ndbd_p16_bc)
        dvars['ndbd_p90_bc']  = (dim_names, ndbd_p90_bc)
    for var_name, mode_arr in comp_modes.items():
        if var_name != 'ndbd':
            dvars[f'{var_name}_mode'] = (dim_names, mode_arr)
    return dvars


def _peak_mem_mb():
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return kb / 1024.0


# ---------------------------------------------------------------------------
# SSCL runner
# ---------------------------------------------------------------------------

def run_sscl(args):
    """Sigma-scaling systematic: fixed +1σ, one fit per dataset."""
    _t0 = start_timer()
    ctx = _setup_common(args)

    if getattr(args, 'save_traces', False):
        out_dir = args.output_dir or ctx['syst_dir']
        cs_tag = str(args.co_shift).replace('.', 'p')
        ctx['save_trace_dir'] = os.path.join(out_dir, f'traces_sscl_cs{cs_tag}')
        os.makedirs(ctx['save_trace_dir'], exist_ok=True)

    _set_worker_shared(ctx, args)

    datasets = ctx['datasets']
    n_ds = ctx['n_ds']

    sscl_scfg = {'type': 'sscl', 'spt': 1}
    jobs = [('sscale_1r', sscl_scfg, ds) for ds in datasets]
    n_total = len(jobs)

    njobs = args.njobs if args.njobs > 0 else max(1, (os.cpu_count() or 1) // 2)
    njobs = min(njobs, n_total)
    print(f"[sscl] {n_total} fits with {njobs} workers", flush=True)

    ndbd_mode    = np.full((n_ds, 1), np.nan)
    ndbd_p16     = np.full((n_ds, 1), np.nan)
    ndbd_p90     = np.full((n_ds, 1), np.nan)
    ndbd_mode_bc = np.full((n_ds, 1), np.nan)
    ndbd_p16_bc  = np.full((n_ds, 1), np.nan)
    ndbd_p90_bc  = np.full((n_ds, 1), np.nan)
    comp_modes_all = {}

    ds_to_i = {d: i for i, d in enumerate(datasets)}

    mp_ctx = multiprocessing.get_context('spawn')
    with mp_ctx.Pool(processes=njobs, maxtasksperchild=1,
                     initializer=_worker_pool_init,
                     initargs=(_WORKER_SHARED,)) as pool:
        for _n, (sname, ds, comp_modes, p16, p90,
                 m_bc, p16_bc, p90_bc) in enumerate(
            pool.imap_unordered(_lss_fit_single, jobs), start=1
        ):
            i = ds_to_i[ds]
            for vn, mv in comp_modes.items():
                if vn not in comp_modes_all:
                    comp_modes_all[vn] = np.full((n_ds, 1), np.nan)
                comp_modes_all[vn][i, 0] = mv

            ndbd_mode[i, 0]    = comp_modes.get('ndbd', np.nan)
            ndbd_p16[i, 0]     = p16
            ndbd_p90[i, 0]     = p90
            ndbd_mode_bc[i, 0] = m_bc
            ndbd_p16_bc[i, 0]  = p16_bc
            ndbd_p90_bc[i, 0]  = p90_bc

            if _n % max(1, n_total // 5) == 0 or _n == n_total:
                print(f"[sscl] {_n}/{n_total}", flush=True)

    # Save
    out_dir = args.output_dir or ctx['syst_dir']
    os.makedirs(out_dir, exist_ok=True)

    dvars = _build_xr_datavars(
        ndbd_mode, ndbd_p16, ndbd_p90,
        ndbd_mode_bc, ndbd_p16_bc, ndbd_p90_bc,
        comp_modes_all, ['ds', 'systs'], has_bias=bool(args.bias_csv),
    )
    res = xr.Dataset(
        dvars,
        coords={'ds': datasets, 'systs': ['sscale_1r']},
        attrs={
            'co_cut': args.co_cut, 'ndbd_cutoff': args.ndbd_cutoff,
            'co_shift': args.co_shift, 'environment': args.environment,
        },
    )
    cs_tag = str(args.co_shift).replace('.', 'p')
    out_path = os.path.join(out_dir, f'lss_sscl_results_cs{cs_tag}.nc')
    res.to_netcdf(out_path)
    print(f"\n[sscl] saved to {out_path}", flush=True)
    print(res, flush=True)
    log_elapsed("sscl total", _t0)
    print(f"[sscl] peak memory: {_peak_mem_mb():.0f} MiB", flush=True)


# ---------------------------------------------------------------------------
# DQ runner
# ---------------------------------------------------------------------------

def run_dq(args):
    """Q-residual systematic: N random parameter draws per dataset.

    With ``--pm1`` uses the pre-computed ±1σ parameter sets
    (``dqres_1r`` / ``dqres_1l``) instead of random draws.
    """
    _t0 = start_timer()
    ctx = _setup_common(args)

    if getattr(args, 'save_traces', False):
        out_dir = args.output_dir or ctx['syst_dir']
        cs_tag = str(args.co_shift).replace('.', 'p')
        ctx['save_trace_dir'] = os.path.join(out_dir, f'traces_dq_cs{cs_tag}')
        os.makedirs(ctx['save_trace_dir'], exist_ok=True)

    _set_worker_shared(ctx, args)

    datasets = ctx['datasets']
    n_ds = ctx['n_ds']

    if getattr(args, 'pm1', False):
        # Use pre-loaded fixed ±1σ Q-residual parameter sets
        lss_params = ctx['lss_params']
        dq_random_samples = [
            lss_params['q_pars_sample_1r'],
            lss_params['q_pars_sample_1l'],
        ]
        n_dq = 2
        _dq_snames = ['dqres_1r', 'dqres_1l']
        print("[dq] using fixed ±1σ parameter sets (pm1 mode)", flush=True)
    else:
        from qpym2.io.qls_scaling import sample_dq_random
        n_dq = args.n_dq
        print(f"[dq] sampling {n_dq} random parameter draws ...", flush=True)
        dq_random_samples = sample_dq_random(
            ctx['cfg'], n_samples=n_dq, seed=args.dq_seed,
        )
        _dq_snames = None

    # Build jobs: n_dq samples × n_ds datasets
    dq_scfg = {'type': 'dq'}
    jobs = []
    sample_names = []
    for si in range(n_dq):
        sname = _dq_snames[si] if _dq_snames else f'dq_r{si:04d}'
        sample_names.append(sname)
        qps = dq_random_samples[si]
        for ds in datasets:
            jobs.append((sname, dq_scfg, ds, qps))

    n_total = len(jobs)
    njobs = args.njobs if args.njobs > 0 else max(1, (os.cpu_count() or 1) // 2)
    njobs = min(njobs, n_total)
    print(f"[dq] {n_total} fits ({n_dq} samples × {n_ds} datasets) "
          f"with {njobs} workers", flush=True)

    ndbd_mode    = np.full((n_ds, n_dq), np.nan)
    ndbd_p16     = np.full((n_ds, n_dq), np.nan)
    ndbd_p90     = np.full((n_ds, n_dq), np.nan)
    ndbd_mode_bc = np.full((n_ds, n_dq), np.nan)
    ndbd_p16_bc  = np.full((n_ds, n_dq), np.nan)
    ndbd_p90_bc  = np.full((n_ds, n_dq), np.nan)
    comp_modes_all = {}

    ds_to_i = {d: i for i, d in enumerate(datasets)}
    sname_to_j = {s: j for j, s in enumerate(sample_names)}
    _log_interval = max(1, n_total // 10)
    _t_pool = time.perf_counter()

    mp_ctx = multiprocessing.get_context('spawn')
    with mp_ctx.Pool(processes=njobs, maxtasksperchild=1,
                     initializer=_worker_pool_init,
                     initargs=(_WORKER_SHARED,)) as pool:
        for _n, (sname, ds, comp_modes, p16, p90,
                 m_bc, p16_bc, p90_bc) in enumerate(
            pool.imap_unordered(_lss_fit_single, jobs), start=1
        ):
            i, j = ds_to_i[ds], sname_to_j[sname]

            for vn, mv in comp_modes.items():
                if vn not in comp_modes_all:
                    comp_modes_all[vn] = np.full((n_ds, n_dq), np.nan)
                comp_modes_all[vn][i, j] = mv

            ndbd_mode[i, j]    = comp_modes.get('ndbd', np.nan)
            ndbd_p16[i, j]     = p16
            ndbd_p90[i, j]     = p90
            ndbd_mode_bc[i, j] = m_bc
            ndbd_p16_bc[i, j]  = p16_bc
            ndbd_p90_bc[i, j]  = p90_bc

            if _n % _log_interval == 0 or _n == n_total:
                _elapsed = time.perf_counter() - _t_pool
                print(f"[dq] {_n}/{n_total}  elapsed={_elapsed:.1f}s"
                      f"  peak_mem={_peak_mem_mb():.0f} MiB", flush=True)

    # Save
    out_dir = args.output_dir or ctx['syst_dir']
    os.makedirs(out_dir, exist_ok=True)

    dvars = _build_xr_datavars(
        ndbd_mode, ndbd_p16, ndbd_p90,
        ndbd_mode_bc, ndbd_p16_bc, ndbd_p90_bc,
        comp_modes_all, ['ds', 'samples'], has_bias=bool(args.bias_csv),
    )
    res = xr.Dataset(
        dvars,
        coords={'ds': datasets, 'samples': sample_names},
        attrs={
            'co_cut': args.co_cut, 'ndbd_cutoff': args.ndbd_cutoff,
            'co_shift': args.co_shift, 'environment': args.environment,
            'n_dq': n_dq, 'dq_seed': args.dq_seed if args.dq_seed is not None else -1,
        },
    )
    cs_tag = str(args.co_shift).replace('.', 'p')
    fname = 'lss_dq_pm1_results' if getattr(args, 'pm1', False) else 'lss_dq_results'
    out_path = os.path.join(out_dir, f'{fname}_cs{cs_tag}.nc')
    res.to_netcdf(out_path)
    print(f"\n[dq] saved to {out_path}", flush=True)
    print(res, flush=True)
    log_elapsed("dq total", _t0)
    print(f"[dq] peak memory: {_peak_mem_mb():.0f} MiB", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_common_args(parser):
    """Add arguments shared by both subcommands."""
    parser.add_argument('--environment', '-e', default='m2v5-1')
    parser.add_argument('--analysis', '-a', default='tech_note')
    parser.add_argument('--hyperparams', default=None)
    parser.add_argument('--co-cut', default='gamma10_v150')
    parser.add_argument('--ndbd-cutoff', default='u1740')
    parser.add_argument('--co-shift', type=float, default=1.0)
    parser.add_argument('--lss-data-path', default=None,
                        help='Path to per-dataset LSS ROOT files.')
    parser.add_argument('--ds-start', type=int, default=3601)
    parser.add_argument('--ds-end', type=int, default=3616)
    parser.add_argument('--draws', type=int, default=2500)
    parser.add_argument('--chains', type=int, default=2)
    parser.add_argument('--output-dir', '-o', default=None,
                        help='Output directory (default: systematic_dir from config).')
    parser.add_argument('--mkchain-samples', type=int, default=100)
    parser.add_argument('--njobs', type=int, default=0,
                        help='Parallel workers (0 = os.cpu_count() // 2).')
    parser.add_argument('--bias-csv', default=None,
                        help='Path to batch CSV with bias fit results.')
    parser.add_argument('--xmin-bias', type=float, default=15,
                        help='Minimum ndbd_inj for bias fit (default: 15).')
    parser.add_argument('--save-traces', action='store_true', default=False,
                        help='Save per-fit MCMC traces to <output_dir>/traces/')


def _parse_args():
    parser = argparse.ArgumentParser(
        description='Parallel LSS systematic fits (sscl / dq).',
    )
    sub = parser.add_subparsers(dest='subcommand', required=True)

    p_sscl = sub.add_parser('sscl',
        help='Sigma-scaling: fixed +1σ, one fit per dataset.')
    _add_common_args(p_sscl)

    p_dq = sub.add_parser('dq',
        help='Q-residual: N random parameter draws per dataset.')
    _add_common_args(p_dq)
    p_dq.add_argument('--n-dq', type=int, default=100,
                      help='Number of random dq parameter draws (default: 100).')
    p_dq.add_argument('--dq-seed', type=int, default=None,
                      help='RNG seed for reproducibility.')
    p_dq.add_argument('--pm1', action='store_true', default=False,
                      help='Use fixed ±1σ parameter sets (dqres_1r/dqres_1l) '
                           'instead of random draws. Ignores --n-dq and --dq-seed.')

    return parser.parse_args()


def main():
    args = _parse_args()
    if args.subcommand == 'sscl':
        run_sscl(args)
    elif args.subcommand == 'dq':
        run_dq(args)


if __name__ == '__main__':
    main()
