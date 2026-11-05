"""run_syst_mc.py — Parallel MC (BM-stat) systematic fits.

Samples N normalisation vectors from the JAGS Markov chain perturbed by
the background correlation matrix, then re-fits the full Bayesian model
for each sample via a ``multiprocessing.Pool`` using fork-based COW sharing.

The nominal histograms are built once and shared across all workers; only
the ``norm`` column of the mctable is varied per sample.

Usage
-----
    python -m qpym2.run.run_syst_mc

    python -m qpym2.run.run_syst_mc \\
        --environment tech_note \\
        --n 500 --draws 2500 --chains 2 --njobs 50
"""

import argparse
import multiprocessing
import os
import resource
import sys
import time

import numpy as np
import pandas as pd
import xarray as xr

from qpym2.utils.logging import debug, log, start_timer, log_elapsed


# ---------------------------------------------------------------------------
# Pool worker (module-level for fork-based COW sharing)
# ---------------------------------------------------------------------------

_WORKER_SHARED = {}

def mc_fit_single_wrapper(
        mcnorm_dict,
        mctable,
        data,
        signal_df,
        hm,
        comp_cfg,
        bias_popt=None,
        bias_pcov=None,
        chains=2,
        draws=5000,
        sample_idx=0,
        save_trace_dir=None):
    """Wrapper to call _mc_fit_single with a single argument (for Pool.map)."""

    job = (sample_idx, mcnorm_dict)
    ws = {
        'mctable':   mctable,
        'data':      data,
        'signal_df': signal_df,
        'hm':        hm,
        'comp_cfg':  comp_cfg,
        'chains':    chains,
        'draws':     draws,
        'bias_popt': bias_popt,
        'bias_pcov': bias_pcov,
        'save_trace_dir': save_trace_dir,
    }
    return _mc_fit_single(job, ws=ws)

def _mc_fit_single(job, ws=None):
    """Run one BM-stat systematic fit inside a Pool worker.

    Reads shared data from the module-level ``_WORKER_SHARED`` dict
    (inherited via ``fork``) to avoid pickling large arrays.
    """
    if ws is None:
        ws = _WORKER_SHARED
    import os
    import time
    
    # Set PYTENSOR_FLAGS before importing to avoid "already initialized" config lock.
    _worker_dir = os.path.join(
        os.environ.get('TMPDIR', '/tmp'),
        f'.pytensor_worker_{os.getpid()}',
    )
    os.makedirs(_worker_dir, exist_ok=True)
    os.environ['PYTENSOR_FLAGS'] = f'compiledir={_worker_dir}'

    import pytensor

    # Ensure ArviZ cache directory exists (prevents FileNotFoundError in forked workers).
    import pathlib
    _arviz_cache = pathlib.Path.home() / '.cache' / 'arviz'
    _arviz_cache.mkdir(parents=True, exist_ok=True)

    sample_idx, mcnorm_dict = job

    from qpym2.fit import fit
    from qpym2.utils.mathy import get_mode

    label = f"sample={sample_idx} pid={os.getpid()}"
    _t_worker_start = time.perf_counter()
    
    # Verify compiledir was set via environment variable
    if hasattr(pytensor.config, 'compiledir'):
        _actual_compiledir = str(pytensor.config.compiledir)
        if not _actual_compiledir.endswith(f'_worker_{os.getpid()}'):
            print(f"  [{label}] WARNING: compiledir mismatch: {_actual_compiledir}", flush=True)
    
    print(f"  [{label}] fitting ...", flush=True)
    _t_fit_start = time.perf_counter()
    
    # Time the shared data copy to check for contention
    _t_copy_start = time.perf_counter()
    mctable = ws['mctable'].copy()
    _t_copy_elapsed = time.perf_counter() - _t_copy_start
    # Align by index so source order never matters
    mctable['norm'] = pd.Series(mcnorm_dict).reindex(mctable.index)

    m2model, _ = fit.setup_fit(
        mctable, ws['data'], ws['signal_df'], ws['hm'],
        comp_cfg=ws['comp_cfg'],
    )
    trace = fit.run_fit(
        m2model, single_core=True,
        chains=ws['chains'], draws=ws['draws'],
    )

    if ws.get('save_trace_dir'):
        trace_path = os.path.join(ws['save_trace_dir'], f"trace_mc{sample_idx:04d}.nc")
        trace.to_netcdf(trace_path)
        print(f"  [{label}] saved trace to {trace_path}", flush=True)

    # Extract modes for all components in the posterior
    component_modes = {}
    for var_name in trace.posterior.data_vars:
        samples = trace.posterior[var_name].values.flatten()
        component_modes[var_name] = float(get_mode(samples))
    
    ndbd_samples = trace.posterior['ndbd'].values.flatten()
    p16 = float(np.percentile(ndbd_samples, 16))
    p90 = float(np.percentile(ndbd_samples, 90))

    # Bias-corrected ndbd stats (if bias params available)
    mode_bc = p16_bc = p90_bc = np.nan
    if ws.get('bias_popt') is not None:
        from qpym2.fit.bias import apply_bias_correction
        corrected = apply_bias_correction(
            ndbd_samples, ws['bias_popt'],
        )
        mode_bc = float(get_mode(corrected))
        p16_bc  = float(np.percentile(corrected, 16))
        p90_bc  = float(np.percentile(corrected, 90))

    _t_fit_elapsed = time.perf_counter() - _t_fit_start
    _t_total_elapsed = time.perf_counter() - _t_worker_start
    _startup_time = _t_fit_start - _t_worker_start
    print(f"  [{label}] ndbd_mode={component_modes.get('ndbd', np.nan):.4f} (startup={_startup_time:.1f}s, fit={_t_fit_elapsed:.1f}s, copy={_t_copy_elapsed:.3f}s, total={_t_total_elapsed:.1f}s)", flush=True)
    return (sample_idx, component_modes, p16, p90, mode_bc, p16_bc, p90_bc)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Parallel MC (BM-stat) systematic fits (Pool-dispatched).',
    )
    parser.add_argument('--environment', '-e', default=None)
    parser.add_argument('--analysis', '-a', default=None)
    parser.add_argument('--model-key', '-m', default='default')
    parser.add_argument('--co-cut', default='gamma10_v150')
    parser.add_argument('--ndbd-cutoff', default='u1740')
    parser.add_argument('--co-shift', type=float, default=1.0)
    parser.add_argument('--n', type=int, default=500,
                        help='Number of BM-stat samples to fit (default: 500).')
    parser.add_argument('--draws', type=int, default=2500)
    parser.add_argument('--chains', type=int, default=2)
    parser.add_argument('--output', '-o', default=None)
    parser.add_argument('--mkchain-samples', type=int, default=100)
    parser.add_argument('--njobs', type=int, default=0,
                        help='Parallel workers (0 = os.cpu_count() // 4).')
    parser.add_argument('--bias-csv', default=None,
                        help='Path to batch CSV with bias fit results. '
                             'If provided, bias-corrected ndbd stats are '
                             'computed inside each worker.')
    parser.add_argument('--xmin-bias', type=float, default=15,
                        help='Minimum ndbd_inj for bias fit (default: 15).')
    parser.add_argument('--save-traces', action='store_true', default=False,
                        help='Save per-sample MCMC traces to '
                             '<syst_dir>/<model_name>/trace_mc{idx:04d}.nc.')
    parser.add_argument('--stat-only', action='store_true', default=False,
                        help='Sample raw JAGS chain without covariance '
                             'perturbation (stat-only, no BM systematic).')
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _peak_mem_mb():
    """Return peak RSS memory usage of this process in MiB (Linux: ru_maxrss is kB)."""
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return kb / 1024.0


def main():
    args = _parse_args()
    _t0 = start_timer()

    # -- 0. Configuration --
    from qpym2.cfg.hydra_adapter import get_config
    from qpym2.cfg.config_builders import (
        get_h2uv_model, get_jags_path, get_fit_dir,
        build_input_config, get_comp_cfg, get_corrmat_path,
    )
    from qpym2.cfg.model_utils import build_model_dir, get_shift_suffix
    from qpym2.cfg.syst_config import get_systematic_dir
    from qpym2.fit.params import model_cfg
    from qpym2.io.bkg_model import read_mkchain, get_mcsyst_samples
    from qpym2.fit.create_m2uv import create_evlist, get_mctable
    from qpym2.fit.bias import read_bias_csv, apply_bias_correction

    cfg = get_config(
        environment=args.environment,
        analysis=args.analysis,
    )

    jags_fpath = get_jags_path(cfg)
    corrmat_fpath = None if args.stat_only else get_corrmat_path(cfg)
    h2uv_model = get_h2uv_model(cfg, model_key=args.model_key)
    comp_cfg_base = get_comp_cfg(cfg, model_key=args.model_key)
    input_cfg = build_input_config(cfg, model_key=args.model_key)
    syst_dir = get_systematic_dir(cfg)

    mcfg = model_cfg(
        h2uv_model,
        co_cut=args.co_cut,
        ndbd_cutoff=args.ndbd_cutoff,
        co_shift=args.co_shift,
    )
    fit_cfg = {'h2uv_model': mcfg['h2uv_model'], 'cuts': mcfg['cuts']}

    # -- 1. Load JAGS chain once; derive BM-stat samples and priors from it --
    print(f"[run_syst_mc] reading JAGS chain ...", flush=True)
    chain = read_mkchain(jags_fpath)

    # BM-stat systematic: MVN-perturbed normalization vectors (used run_syst_mc only)
    print(f"[run_syst_mc] sampling {args.n} BM-stat points ...", flush=True)
    mc_samples = get_mcsyst_samples(chain, corrmat_fpath=corrmat_fpath, n=args.n)
    # mc_samples: DataFrame(n, n_sources), columns = JAGS source names

    # -- 3. Read priors from nominal fit (matches run_syst_lss) --
    fit_dir = get_fit_dir(cfg)
    model_dir = build_model_dir(
        co_cut=args.co_cut, ndbd_cutoff=args.ndbd_cutoff,
    )
    shift_suffix = get_shift_suffix(co_shift=args.co_shift)
    nominal_fit_h5 = os.path.join(fit_dir, model_dir, f'fit{shift_suffix}.h5')
    print(f"[run_syst_mc] reading nominal fit priors from {nominal_fit_h5} ...", flush=True)
    nominal_comps = pd.read_hdf(nominal_fit_h5, key='components')
    nominal_priors = list(nominal_comps['prior'])
    print(f"[run_syst_mc] nominal priors ({len(nominal_priors)} components): "
          f"{[p[0] if isinstance(p, tuple) else type(p).__name__ for p in nominal_priors]}",
          flush=True)

    # -- 4. Build nominal mctable + ndbd signal_df (histograms read once) --
    print("[run_syst_mc] building nominal mctable ...", flush=True)
    staging_comps = get_mctable(fit_cfg, input_cfg, jags_fpath,
                                fit_pars=mcfg['fit_pars'])
    mctable = staging_comps.loc[staging_comps.index != 'ndbd']
    signal_df = staging_comps.loc[staging_comps.index == 'ndbd']

    # -- 5. Inject nominal priors into comp_cfg --
    comp_cfg = list(comp_cfg_base)
    for i, comp in enumerate(comp_cfg):
        prior_type, prior_params = nominal_priors[i]
        comp['mean'] = prior_params['mu']
        comp['std'] = prior_params['sigma']

    # -- 6. Read observed data --
    print("[run_syst_mc] reading data ...", flush=True)
    data = create_evlist(mcfg, input_cfg)

    # -- 7. Build job list --
    jobs = [
        (i, mc_samples.iloc[i].to_dict())
        for i in range(len(mc_samples))
    ]
    njobs = args.njobs if args.njobs > 0 else max(1, (os.cpu_count() or 1) // 4)
    njobs = min(njobs, len(jobs))
    print(f"[run_syst_mc] os.cpu_count()={os.cpu_count()}, njobs={njobs}",
          flush=True)
    print(f"[run_syst_mc] launching {len(jobs)} fits with {njobs} workers ...",
          flush=True)
    _t_pool_start = start_timer()

    ndbd_mode = np.full(len(jobs), np.nan)
    ndbd_p16  = np.full(len(jobs), np.nan)
    ndbd_p90  = np.full(len(jobs), np.nan)
    ndbd_mode_bc = np.full(len(jobs), np.nan)
    ndbd_p16_bc  = np.full(len(jobs), np.nan)
    ndbd_p90_bc  = np.full(len(jobs), np.nan)
    component_modes_all = {}  # dict of {var_name: np.array(n_samples)}

    # -- 7b. Read bias correction parameters (if requested) --
    bias_popt = None
    bias_pcov = None
    if args.bias_csv:
        model_pars = {
            'co_cut': args.co_cut,
            'ndbd_cutoff': args.ndbd_cutoff,
            'co_shift': args.co_shift,
        }
        bias_popt, bias_pcov = read_bias_csv(
            args.bias_csv, model_pars, xmin=args.xmin_bias,
        )
        print(f"[run_syst_mc] bias correction: a={bias_popt[0]:.4f}, "
              f"b={bias_popt[1]:.4f}", flush=True)

    # -- 7c. Build trace output directory (if --save-traces requested) --
    save_trace_dir = None
    if args.save_traces:
        co_shift_tag = str(args.co_shift).replace('.', 'p')
        save_trace_dir = os.path.join(syst_dir, f'traces_mc_cs{co_shift_tag}')
        os.makedirs(save_trace_dir, exist_ok=True)
        print(f"[run_syst_mc] saving traces to {save_trace_dir}", flush=True)

    global _WORKER_SHARED
    _WORKER_SHARED = {
        'mctable':   mctable,
        'signal_df': signal_df,
        'data':      data,
        'hm':        fit_cfg['h2uv_model'],
        'comp_cfg':  comp_cfg,
        'chains':    args.chains,
        'draws':     args.draws,
        'bias_popt': bias_popt,
        'bias_pcov': bias_pcov,
        'save_trace_dir': save_trace_dir,
    }

    _n_total = len(jobs)
    _log_interval = max(1, _n_total // 10)  # log every ~10 %
    ctx = multiprocessing.get_context('fork')
    with ctx.Pool(processes=njobs) as pool:
        for _n_done, (sample_idx, comp_modes, p16, p90,
                      m_bc, p16_bc, p90_bc) in enumerate(
            pool.imap_unordered(_mc_fit_single, jobs), start=1
        ):
            # Initialize arrays for new components on first encounter
            for var_name, mode_val in comp_modes.items():
                if var_name not in component_modes_all:
                    component_modes_all[var_name] = np.full(len(jobs), np.nan)
                component_modes_all[var_name][sample_idx] = mode_val
            
            ndbd_mode[sample_idx] = comp_modes.get('ndbd', np.nan)
            ndbd_p16[sample_idx]  = p16
            ndbd_p90[sample_idx]  = p90
            ndbd_mode_bc[sample_idx] = m_bc
            ndbd_p16_bc[sample_idx]  = p16_bc
            ndbd_p90_bc[sample_idx]  = p90_bc
            if _n_done % _log_interval == 0 or _n_done == _n_total:
                _elapsed = time.perf_counter() - _t_pool_start
                _mem_mb = _peak_mem_mb()
                print(
                    f"[run_syst_mc] progress: {_n_done}/{_n_total}"
                    f"  elapsed={_elapsed:.1f}s"
                    f"  peak_mem={_mem_mb:.0f} MiB",
                    flush=True,
                )

    # -- 8. Assemble and save results --
    # Build data_vars dict with all component modes
    data_vars = {
        'ndbd_mode': (['sample'], ndbd_mode),
        'ndbd_p16':  (['sample'], ndbd_p16),
        'ndbd_p90':  (['sample'], ndbd_p90),
    }
    if args.bias_csv:
        data_vars['ndbd_mode_bc'] = (['sample'], ndbd_mode_bc)
        data_vars['ndbd_p16_bc']  = (['sample'], ndbd_p16_bc)
        data_vars['ndbd_p90_bc']  = (['sample'], ndbd_p90_bc)
    # Add modes for all other components
    for var_name, mode_arr in component_modes_all.items():
        if var_name != 'ndbd':  # ndbd_mode already added above
            data_vars[f'{var_name}_mode'] = (['sample'], mode_arr)
    
    sres = xr.Dataset(
        data_vars,
        coords={'sample': np.arange(len(jobs))},
        attrs={
            'co_cut':      args.co_cut,
            'ndbd_cutoff': args.ndbd_cutoff,
            'co_shift':    args.co_shift,
            'environment': args.environment,
            'n_samples':   args.n,
        },
    )

    co_shift_tag = str(args.co_shift).replace('.', 'p')
    out_path = args.output or os.path.join(syst_dir, f'mc_results_cs{co_shift_tag}.nc')
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    sres.to_netcdf(out_path)
    print(f"[run_syst_mc] saved results to {out_path}", flush=True)
    print(sres, flush=True)
    log_elapsed("run_syst_mc total", _t0)
    print(f"[run_syst_mc] peak memory: {_peak_mem_mb():.0f} MiB", flush=True)


if __name__ == '__main__':
    main()
