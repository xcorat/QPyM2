"""
Run toy fits for likelihood validation and bias studies.

This script runs toy, sensitivity, and bias fits using the Hydra
configuration system. It replaces the old run_old/run_toy_fits.py
which used hardcoded imports and the old model_cfg dict pattern.

Usage:
    python run_toy_fits.py [--environment ENV] [--analysis ANALYSIS]
    python run_toy_fits.py --type sens --n 200
    python run_toy_fits.py --type bias --bias-min 5 --bias-max 200
"""

import argparse
import numpy as np
import pandas as pd
import xarray as xr

from qpym2.cfg.hydra_adapter import get_config
from qpym2.cfg.config_builders import (
    get_fit_dir, get_h2uv_model, get_data_names
)
from qpym2.cfg.params_grid import ALL_CO_CUTS, ALL_NDBD_CUTS, ALL_CO_SHIFTS
from qpym2.utils.mathy import get_mode, find_mode
from qpym2.utils.logging import debug
from qpym2.io import fit as fit_io
from qpym2.io.toys import write_toys
from qpym2.io.csv import append_csv, flatten_columns
from qpym2.fit.posterior_checks import toy_fits, sens_fits, bias_fits
from qpym2.fit.results import _get_stats
from qpym2.fit.params import model_cfg


# --- Core fit functions -------------------------------------------------------

def _get_model_outdir(mcfg, fit_dir):
    """Return the model output directory for the current file layout."""
    return f"{fit_dir}/{mcfg.get('dir_name', mcfg['name'])}"


def _get_fit_fname(mcfg, include_blinding=False):
    """Build fit filename stem using current shift-suffix layout."""
    fname = 'fit'
    if mcfg['skip_roi']:
        fname += '_noroi'
    if include_blinding and mcfg['blinded']:
        fname += '_blinded'
    fname += mcfg.get('shift_suffix', '')
    return fname


def _get_toys_fname(mcfg):
    """Build toy output filename stem using current shift-suffix layout."""
    fname = 'toys'
    if mcfg['skip_roi']:
        fname += '_noroi'
    if mcfg['blinded']:
        fname += '_blinded'
    fname += mcfg.get('shift_suffix', '')
    return fname

def run_toy_fit(mcfg, fit_dir, n=1, write=False, csv_path=None):
    """Run toy fits for likelihood validation.

    Generates toy data from the posterior mode and re-fits to validate
    the likelihood.

    Args:
        mcfg (dict): Model configuration dict (from ``model_cfg`` or
            ``get_model_cfg``).
        fit_dir (str): Base directory containing fit results.
        n (int): Number of toy fits to run.
        write (bool): If True, persist results to HDF5.
        csv_path (str | None): If given, append flattened stats to this
            CSV file (concurrency-safe via ``fcntl`` locking).

    Returns:
        pd.DataFrame: Fit results with multi-level stat columns per
            variable, plus ``logp_mode``, ``ndbd_inj``, ``model``, and
            ``type='toys'``.
    """
    model_outdir = _get_model_outdir(mcfg, fit_dir)
    comps, trace, _ = fit_io.read_fit(model_outdir, fname=_get_fit_fname(mcfg))
    vars = list(comps.index)

    toys = toy_fits(comps, trace, mcfg, nsamples=n)

    fit_vals = np.array([[_get_stats(t.posterior[var].values.flatten())
                          for var in vars]
                         for t in toys])

    fv_multi = []
    for i, fv in enumerate(fit_vals):
        sd = {}
        for var, stat in zip(vars, fv):
            for s, v in stat.items():
                sd[(var, s)] = v
        fv_multi.append(pd.Series(sd, name=i))

    fit_res = pd.concat(fv_multi, axis=1).T
    fit_res['logp_mode'] = [get_mode(t.log_likelihood.pik.values.flatten()) for t in toys]
    fit_res['ndbd_inj'] = None
    fit_res['model'] = mcfg['name']
    fit_res['type'] = 'toys'

    if write:
        write_toys(fit_res, model_outdir, type='fit', fname=_get_toys_fname(mcfg))

    if csv_path is not None:
        append_csv(flatten_columns(fit_res), csv_path)

    return fit_res


def run_sens_fit(mcfg, fit_dir, n=1, signal_comp='ndbd', write=False, csv_path=None):
    """Run sensitivity fits for likelihood validation.

    Like toy fits but with the signal component set to zero when
    generating toy data.

    Args:
        mcfg (dict): Model configuration dict.
        fit_dir (str): Base directory containing fit results.
        n (int): Number of toy fits to run.
        signal_comp (str): Name of the signal component.
        write (bool): If True, persist results to HDF5.
        csv_path (str | None): If given, append flattened stats to this
            CSV file (concurrency-safe via ``fcntl`` locking).

    Returns:
        pd.DataFrame: Fit results with multi-level stat columns,
            ``type='sens'``.
    """
    model_outdir = _get_model_outdir(mcfg, fit_dir)
    comps, trace, _ = fit_io.read_fit(model_outdir, fname=_get_fit_fname(mcfg))
    vars = list(comps.index)

    toys = sens_fits(comps, trace, mcfg, nsamples=n, signal_comp=signal_comp)
    print('sens', mcfg['name'])

    fit_vals = np.array([[_get_stats(t.posterior[var].values.flatten())
                          for var in vars]
                         for t in toys])

    fv_multi = []
    for i, fv in enumerate(fit_vals):
        sd = {}
        for var, stat in zip(vars, fv):
            for s, v in stat.items():
                sd[(var, s)] = v
        fv_multi.append(pd.Series(sd, name=i))

    fit_res = pd.concat(fv_multi, axis=1).T
    fit_res['logp_mode'] = [get_mode(t.log_likelihood.pik.values.flatten()) for t in toys]
    fit_res['ndbd_inj'] = 0
    fit_res['model'] = mcfg['name']
    fit_res['type'] = 'sens'

    if write:
        write_toys(fit_res, model_outdir, type='sens', fname=_get_toys_fname(mcfg))

    if csv_path is not None:
        append_csv(flatten_columns(fit_res), csv_path)

    return fit_res


def run_bias_fit(mcfg, fit_dir, n=1, signal_comp='ndbd',
                 bias_min=0, bias_max=None, write=False, csv_path=None):
    """Run bias fits for signal injection studies.

    Generates toy data with a random uniform signal injection and re-fits
    to check for potential biases.

    Args:
        mcfg (dict): Model configuration dict.
        fit_dir (str): Base directory containing fit results.
        n (int): Number of toy fits to run.
        signal_comp (str): Name of the signal component.
        bias_min (float): Minimum injected signal counts.
        bias_max (float | None): Maximum injected signal counts.
            If None, defaults to 2× the posterior mean.
        write (bool): If True, persist results to HDF5.
        csv_path (str | None): If given, append flattened stats to this
            CSV file (concurrency-safe via ``fcntl`` locking).

    Returns:
        pd.DataFrame: Fit results with multi-level columns for stats,
            plus ``logp_mode``, ``ndbd_inj``, ``model``, and ``type='bias'``.
    """
    model_outdir = _get_model_outdir(mcfg, fit_dir)
    comps, trace, _ = fit_io.read_fit(
        model_outdir,
        fname=_get_fit_fname(mcfg, include_blinding=True),
    )
    vars = list(comps.index)

    if bias_max is None:
        bias_max = 2 * trace.posterior[signal_comp].mean().item()

    toys, sig_inj = bias_fits(comps, trace, mcfg, nsamples=n,
                              signal_comp=signal_comp,
                              bias_min=bias_min, bias_max=bias_max)
    print('bias', mcfg['name'])

    fit_vals = np.array([[_get_stats(t.posterior[var].values.flatten())
                          for var in vars]
                         for t in toys])

    # Extract multi-level stats into a dataframe
    fv_multi = []
    for i, fv in enumerate(fit_vals):
        sd = {}
        dd = [{(var, s): v for s, v in stat.items()} for var, stat in zip(vars, fv)]
        for d in dd:
            sd.update(d)
        fv_multi.append(pd.Series(sd, name=i))

    fit_res = pd.concat(fv_multi, axis=1).T

    fit_res['logp_mode'] = [get_mode(t.log_likelihood.pik.values.flatten()) for t in toys]
    fit_res['ndbd_inj'] = sig_inj
    fit_res['model'] = mcfg['name']
    fit_res['type'] = 'bias'

    if write:
        write_toys(fit_res, model_outdir, type='bias', fname=_get_toys_fname(mcfg))

    if csv_path is not None:
        append_csv(flatten_columns(fit_res), csv_path)

    return fit_res


# --- Entry points using Hydra config -----------------------------------------

def _build_model_grid(cfg, h2uv_model, co_cuts, ndbd_cutoffs, co_shifts,
                      skip_roi=False, blinded=False):
    """Build an xarray Dataset of model_cfg dicts over a parameter grid.

    Uses ``model_cfg`` with cut dicts extracted from the Hydra config.

    Args:
        cfg: Hydra DictConfig.
        h2uv_model: H2UV model object / dict.
        co_cuts: List of co calibration cut keys (or None).
        ndbd_cutoffs: List of ndbd cutoff keys (or None).
        co_shifts: List of co shift values.
        skip_roi: ROI to skip (False or tuple).
        blinded: Blinding flag.

    Returns:
        xr.Dataset with a ``'cfg'`` variable.
    """
    co_cuts_dict = dict(cfg.cuts.co_cuts) if 'cuts' in cfg and 'co_cuts' in cfg.cuts else None
    ndbd_cuts_dict = dict(cfg.cuts.ndbd_cutoff) if 'cuts' in cfg and 'ndbd_cutoff' in cfg.cuts else None

    fit_models = xr.Dataset(
        coords={
            'ndbd_cutoff': ndbd_cutoffs,
            'co_cuts': co_cuts,
            'co_shift': co_shifts,
        },
    )

    smooth = getattr(h2uv_model, 'smooth', True)

    fit_models['cfg'] = xr.apply_ufunc(
        lambda coc, ndc, cos: model_cfg(
            h2uv_model,
            co_cut=coc,
            ndbd_cutoff=ndc,
            co_shift=cos,
            skip_roi=skip_roi,
            blinded=blinded,
            smooth=smooth,
            co_cuts_dict=co_cuts_dict,
            ndbd_cuts_dict=ndbd_cuts_dict),
        fit_models.co_cuts, fit_models.ndbd_cutoff, fit_models.co_shift,
        vectorize=True)

    return fit_models


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run toy / sensitivity / bias fits for M2 analysis.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            python run_toy_fits.py
            python run_toy_fits.py --environment tech_note --type sens --n 200
            python run_toy_fits.py --type bias --bias-min 5 --bias-max 200 --n 500
        """
    )

    parser.add_argument('--environment', '-e', default='tech_note',
                        help='Environment name (default: tech_note)')
    parser.add_argument('--analysis', '-a', default=None,
                        help='Analysis profile name (default: same as environment)')
    parser.add_argument('--hyperparams', default=None,
                        help='Hyperparameter config name')
    parser.add_argument('--type', '-t', default='all',
                        choices=['toys', 'sens', 'bias', 'all'],
                        help='Type of toy fit to run (default: all)')
    parser.add_argument('--n', type=int, default=200,
                        help='Number of toy fits per model (default: 200)')
    parser.add_argument('--bias-min', type=float, default=5,
                        help='Min injected signal for bias fits (default: 5)')
    parser.add_argument('--bias-max', type=float, default=200,
                        help='Max injected signal for bias fits (default: 200)')
    parser.add_argument('--write', action='store_true', default=True,
                        help='Write results to HDF5 (default: True)')

    return parser.parse_args()


def main():
    """Run toy fits over a model parameter grid using Hydra config."""
    args = parse_args()

    # Load Hydra configuration
    analysis = args.analysis if args.analysis else args.environment
    cfg = get_config(
        environment=args.environment,
        analysis=analysis,
        hyperparams=args.hyperparams,
    )

    # Derive paths and model from config
    fit_dir = get_fit_dir(cfg)
    h2uv_model = get_h2uv_model(cfg)

    from qpym2.utils.logging import DEBUG
    DEBUG = 0

    # Parameter grid — matches old main_all_shifts default
    co_cuts      = ALL_CO_CUTS[3:4]    # ['gamma10_v150']
    ndbd_cutoffs = ALL_NDBD_CUTS[2:3]  # ['u1740']
    co_shifts    = ALL_CO_SHIFTS

    fit_models = _build_model_grid(cfg, h2uv_model, co_cuts, ndbd_cutoffs,
                                   co_shifts)

    n = args.n
    for mcfg in fit_models.cfg.values.flatten():
        debug('running:', mcfg['name'])

        if args.type in ('toys', 'all'):
            res = run_toy_fit(mcfg, fit_dir, n=n, write=args.write)
            print('toy dataframe:', res.shape)

        if args.type in ('sens', 'all'):
            res_s = run_sens_fit(mcfg, fit_dir, n=n, write=args.write)
            print('sens dataframe:', res_s.shape)

        if args.type in ('bias', 'all'):
            res_b = run_bias_fit(mcfg, fit_dir, n=n,
                                 bias_min=args.bias_min,
                                 bias_max=args.bias_max,
                                 write=args.write)
            print('bias dataframe:', res_b.shape)


if __name__ == '__main__':
    main()
