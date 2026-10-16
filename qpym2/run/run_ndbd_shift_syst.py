"""run_ndbd_shift_syst.py — 0νββ signal energy-shift systematic fits.

Shifts the ndbd signal component by each value in ``ALL_CO_SHIFTS``
(0.0–1.5 keV) while keeping background at the nominal (unshifted) model,
refits, and saves traces + aggregated xarray results.

Usage:
    python -m qpym2.run.run_ndbd_shift_syst
    python -m qpym2.run.run_ndbd_shift_syst --draws 500 --chains 2
    python -m qpym2.run.run_ndbd_shift_syst --co-cut gamma10_v150 --ndbd-cutoff u1740
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
    """Initializer for spawn pool workers: receive shared data by pickle."""
    global _WORKER_SHARED
    _WORKER_SHARED = shared


def _ndbd_shift_fit_single(ndbd_shift, ws=None):
    """Run one ndbd-shift systematic fit inside a Pool worker.

    Parameters
    ----------
    ndbd_shift : float
        Energy shift in keV to apply to the ndbd signal component.
    ws : dict or None
        Shared workspace (uses module-level ``_WORKER_SHARED`` if None).

    Returns
    -------
    tuple
        ``(ndbd_shift, component_modes, p16, p90, mode_bc, p16_bc, p90_bc)``
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

    import time as _time
    from qpym2.fit.syst import make_ndbd_shift_signal_df
    from qpym2.fit.syst_lss import (
        run_lss_fit,
        extract_posterior_stats, apply_bias_correction,
    )

    ns_tag = str(ndbd_shift).replace('.', 'p')
    label = f"ndbd_shift={ndbd_shift}"
    _t0 = _time.perf_counter()
    print(f"  [{label}] fitting ...", flush=True)

    try:
        # Build shifted signal_df
        signal_df = make_ndbd_shift_signal_df(
            ws['evdata_ndbd'], ndbd_shift, ws['hm'],
            filters=ws['cuts_ndbd'],
        )

        # Run fit using nominal mctable + shifted signal
        trace = run_lss_fit(
            ws['mctable'], signal_df, ws['data'],
            ws['hm'], ws['comp_cfg'], ws['nominal_priors'],
            chains=ws['chains'], draws=ws['draws'],
        )

        # Save trace if requested
        if ws.get('save_trace_dir'):
            trace_path = os.path.join(
                ws['save_trace_dir'], f'trace_ndbd_shift_{ns_tag}.nc',
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
        return (ndbd_shift, component_modes, p16, p90, mode_bc, p16_bc, p90_bc)

    except Exception as e:
        _elapsed = _time.perf_counter() - _t0
        print(f"  [{label}] FAILED ({_elapsed:.1f}s): {e}", flush=True)
        return (ndbd_shift, {}, np.nan, np.nan, np.nan, np.nan, np.nan)


# ---------------------------------------------------------------------------
# Shared setup
# ---------------------------------------------------------------------------

def _setup_common(args):
    """Load config, nominal mctable, data, and ndbd event arrays.

    Returns a dict with everything the ndbd-shift worker needs.
    """
    from qpym2.cfg.hydra_adapter import get_config
    from qpym2.cfg.config_builders import (
        get_h2uv_model, get_jags_path, get_fit_dir,
        build_input_config, get_comp_cfg,
    )
    from qpym2.cfg.model_utils import build_model_dir, get_shift_suffix
    from qpym2.cfg.syst_config import get_systematic_dir
    from qpym2.fit.params import model_cfg
    from qpym2.fit.create_evtrees import create_evtrees
    from qpym2.fit.create_m2uv import create_evlist
    from qpym2.fit.syst import build_nominal_mctable
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

    # Build model config (no co_shift on background for this systematic)
    mcfg = model_cfg(
        h2uv_model,
        co_cut=args.co_cut,
        ndbd_cutoff=args.ndbd_cutoff,
        co_shift=0.0,  # nominal background, no shift
    )
    fit_cfg = {'h2uv_model': mcfg['h2uv_model'], 'cuts': mcfg['cuts']}
    hm = fit_cfg['h2uv_model']

    # Background model
    print("[setup] reading background model ...", flush=True)
    mcfit_table, _ = read_jagsh5(jags_fpath)

    # Build nominal mctable (background histograms from staging)
    print("[setup] building nominal mctable ...", flush=True)
    mctable, _ = build_nominal_mctable(mcfit_table, input_cfg, fit_cfg)

    # Data
    print("[setup] reading data ...", flush=True)
    data = create_evlist(mcfg, input_cfg)

    # Read nominal fit priors
    print("[setup] reading nominal fit priors ...", flush=True)
    fit_dir = get_fit_dir(cfg)
    model_dir = build_model_dir(
        co_cut=args.co_cut, ndbd_cutoff=args.ndbd_cutoff,
    )
    shift_suffix = get_shift_suffix(co_shift=0.0)
    nominal_fit_h5 = os.path.join(fit_dir, model_dir, f'fit{shift_suffix}.h5')
    nominal_comps = pd.read_hdf(nominal_fit_h5, key='components')
    nominal_priors = list(nominal_comps['prior'])
    print(f"[setup] nominal priors ({len(nominal_priors)} components): "
          f"{[p[0] if isinstance(p, tuple) else type(p).__name__ for p in nominal_priors]}",
          flush=True)

    comp_cfg = list(comp_cfg_base)

    # Materialize ndbd events to numpy (only ndbd needed)
    print("[setup] materializing ndbd events ...", flush=True)
    _NP_COLS = ['u0', 'v0', 'ch1', 'ch2', 'ds']
    comp_fnames = {'ndbd': 'ndbd.root'}
    evtrees = create_evtrees(input_cfg, comp_fnames)
    ndbd_rdf = evtrees.get('ndbd')
    if ndbd_rdf is not None:
        evdata_ndbd = {c: np.array(a) for c, a
                       in ndbd_rdf.AsNumpy(columns=_NP_COLS).items()}
    else:
        evdata_ndbd = None
    del evtrees

    gc.collect()
    gc.collect()

    # Pre-compute ndbd cuts
    cuts = fit_cfg['cuts']
    cuts_ndbd = [f for f in [cuts.get('common'), cuts.get('ndbd_only')] if f]

    # Bias correction
    bias_popt = bias_pcov = None
    if args.bias_csv:
        model_pars = {
            'co_cut': args.co_cut,
            'ndbd_cutoff': args.ndbd_cutoff,
            'co_shift': 0.0,
        }
        bias_popt, bias_pcov = read_bias_csv(
            args.bias_csv, model_pars, xmin=args.xmin_bias,
        )
        print(f"[setup] bias correction: a={bias_popt[0]:.4f}, "
              f"b={bias_popt[1]:.4f}", flush=True)

    return {
        'cfg': cfg,
        'mctable': mctable,
        'data': data,
        'hm': hm,
        'fit_cfg': fit_cfg,
        'comp_cfg': comp_cfg,
        'nominal_priors': nominal_priors,
        'evdata_ndbd': evdata_ndbd,
        'cuts_ndbd': cuts_ndbd,
        'syst_dir': syst_dir,
        'bias_popt': bias_popt,
        'bias_pcov': bias_pcov,
    }


def _peak_mem_mb():
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return kb / 1024.0


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_ndbd_shift(args):
    """Run ndbd energy-shift systematic: one fit per shift value."""
    from qpym2.cfg.params_grid import ALL_CO_SHIFTS

    _t0 = start_timer()
    ctx = _setup_common(args)

    shift_values = ALL_CO_SHIFTS
    n_shifts = len(shift_values)

    # Trace saving
    out_dir = args.output_dir or ctx['syst_dir']
    os.makedirs(out_dir, exist_ok=True)
    save_trace_dir = None
    if getattr(args, 'save_traces', True):
        save_trace_dir = os.path.join(out_dir, 'traces_ndbd_shift')
        os.makedirs(save_trace_dir, exist_ok=True)

    # Build shared workspace for workers
    global _WORKER_SHARED
    _WORKER_SHARED = {
        'mctable':         ctx['mctable'],
        'data':            ctx['data'],
        'hm':              ctx['hm'],
        'comp_cfg':        ctx['comp_cfg'],
        'nominal_priors':  ctx['nominal_priors'],
        'evdata_ndbd':     ctx['evdata_ndbd'],
        'cuts_ndbd':       ctx['cuts_ndbd'],
        'chains':          args.chains,
        'draws':           args.draws,
        'bias_popt':       ctx['bias_popt'],
        'bias_pcov':       ctx['bias_pcov'],
        'save_trace_dir':  save_trace_dir,
    }

    njobs = args.njobs if args.njobs > 0 else max(1, (os.cpu_count() or 1) // 2)
    njobs = min(njobs, n_shifts)
    print(f"[ndbd_shift] {n_shifts} fits with {njobs} workers", flush=True)

    # Result arrays
    ndbd_mode    = np.full(n_shifts, np.nan)
    ndbd_p16     = np.full(n_shifts, np.nan)
    ndbd_p90     = np.full(n_shifts, np.nan)
    ndbd_mode_bc = np.full(n_shifts, np.nan)
    ndbd_p16_bc  = np.full(n_shifts, np.nan)
    ndbd_p90_bc  = np.full(n_shifts, np.nan)
    comp_modes_all = {}

    shift_to_i = {s: i for i, s in enumerate(shift_values)}

    mp_ctx = multiprocessing.get_context('spawn')
    with mp_ctx.Pool(processes=njobs, maxtasksperchild=1,
                     initializer=_worker_pool_init,
                     initargs=(_WORKER_SHARED,)) as pool:
        for result in pool.imap_unordered(
            _ndbd_shift_fit_single, shift_values
        ):
            (shift_val, comp_modes, p16, p90,
             m_bc, p16_bc, p90_bc) = result
            i = shift_to_i[shift_val]

            for vn, mv in comp_modes.items():
                if vn not in comp_modes_all:
                    comp_modes_all[vn] = np.full(n_shifts, np.nan)
                comp_modes_all[vn][i] = mv

            ndbd_mode[i]    = comp_modes.get('ndbd', np.nan)
            ndbd_p16[i]     = p16
            ndbd_p90[i]     = p90
            ndbd_mode_bc[i] = m_bc
            ndbd_p16_bc[i]  = p16_bc
            ndbd_p90_bc[i]  = p90_bc

    # Build xarray Dataset
    dvars = {
        'ndbd_mode': ('ndbd_shift', ndbd_mode),
        'ndbd_p16':  ('ndbd_shift', ndbd_p16),
        'ndbd_p90':  ('ndbd_shift', ndbd_p90),
    }
    if args.bias_csv:
        dvars['ndbd_mode_bc'] = ('ndbd_shift', ndbd_mode_bc)
        dvars['ndbd_p16_bc']  = ('ndbd_shift', ndbd_p16_bc)
        dvars['ndbd_p90_bc']  = ('ndbd_shift', ndbd_p90_bc)
    for var_name, mode_arr in comp_modes_all.items():
        if var_name != 'ndbd':
            dvars[f'{var_name}_mode'] = ('ndbd_shift', mode_arr)

    res = xr.Dataset(
        dvars,
        coords={'ndbd_shift': shift_values},
        attrs={
            'co_cut': args.co_cut,
            'ndbd_cutoff': args.ndbd_cutoff,
            'environment': args.environment,
        },
    )

    out_path = os.path.join(out_dir, 'ndbd_shift_results.nc')
    res.to_netcdf(out_path)
    print(f"\n[ndbd_shift] saved to {out_path}", flush=True)
    print(res, flush=True)
    log_elapsed("ndbd_shift total", _t0)
    print(f"[ndbd_shift] peak memory: {_peak_mem_mb():.0f} MiB", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description='0νββ signal energy-shift systematic fits.',
    )
    parser.add_argument('--environment', '-e', default='m2v5-1')
    parser.add_argument('--analysis', '-a', default='tech_note')
    parser.add_argument('--hyperparams', default=None)
    parser.add_argument('--co-cut', default='gamma10_v150')
    parser.add_argument('--ndbd-cutoff', default='u1740')
    parser.add_argument('--draws', type=int, default=4000)
    parser.add_argument('--chains', type=int, default=2)
    parser.add_argument('--output-dir', '-o', default=None,
                        help='Output directory (default: systematic_dir from config).')
    parser.add_argument('--njobs', type=int, default=0,
                        help='Parallel workers (0 = os.cpu_count() // 2).')
    parser.add_argument('--bias-csv', default=None,
                        help='Path to batch CSV with bias fit results.')
    parser.add_argument('--xmin-bias', type=float, default=15,
                        help='Minimum ndbd_inj for bias fit (default: 15).')
    parser.add_argument('--save-traces', action='store_true', default=True,
                        help='Save per-shift MCMC traces (default: True).')
    parser.add_argument('--no-save-traces', action='store_false', dest='save_traces',
                        help='Do not save per-shift MCMC traces.')
    return parser.parse_args()


def main():
    args = _parse_args()
    run_ndbd_shift(args)


if __name__ == '__main__':
    main()
