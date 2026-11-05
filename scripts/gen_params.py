#!/usr/bin/env python3
"""Unified parameter-grid generator for all QPyM2 run types.

Subcommands:  staging, input, fit, toys, sens, bias
Each writes a TSV to scripts/params/<run_type>_params.tsv by default.

Pass -i / --interactive to be prompted for grid values; otherwise the
defaults from qpym2.cfg.params_grid are used.  Every option can also
be set directly via flags (flags override interactive selections).

Usage:
    python scripts/gen_params.py staging
    python scripts/gen_params.py fit -i
    python scripts/gen_params.py toys --co-cut gamma10_v150 v150 --n-fits 200
    python scripts/gen_params.py sens --co-cut gamma10_v150 --n-fits 100
    python scripts/gen_params.py input --co-shift 0.8 1.0 1.2
"""

import argparse
import csv
import itertools
import os
import random
import subprocess
import sys
import textwrap

# Ensure the project root is importable (scripts/ lives one level below)
_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from qpym2.cfg.params_grid import (
    ALL_CO_CUTS, ALL_NDBD_CUTS, ALL_CO_SHIFTS,
    CO_CUTS, NDBD_CUTS, CO_SHIFTS,
)

# --------------------------------------------------------------------------- #
#  Constants                                                                   #
# --------------------------------------------------------------------------- #
PARAMS_DIR    = os.path.join('scripts', 'params')
MASTER_SEED   = 47

TOY_SCRIPT    = 'qpym2/run/run_toy_batch.py'
FIT_SCRIPT    = 'qpym2/run/run_single_fit.py'
INPUT_SCRIPT  = 'qpym2/run/run_create_input.py'
STAGING_SCRIPT = 'qpym2/run/run_create_staging.py'
SYST_LSS_SCRIPT = 'qpym2/run/run_syst_lss.py'
SYST_MC_SCRIPT  = 'qpym2/run/run_syst_mc.py'

FIELDNAMES_GRID = [
    'task_id', 'seed', 'co_cut', 'ndbd_cutoff',
    'co_shift', 'fit_type', 'rep', 'command',
]
FIELDNAMES_INPUT = [
    'task_id', 'co_cut', 'ndbd_cutoff', 'co_shifts', 'command',
]
FIELDNAMES_SYST_LSS = [
    'task_id', 'co_shift', 'command',
]
FIELDNAMES_SYST_MC = [
    'task_id', 'co_shift', 'command',
]


# --------------------------------------------------------------------------- #
#  Interactive helpers                                                         #
# --------------------------------------------------------------------------- #

def _format_value(v):
    """Display-friendly label for a grid value (None → 'None')."""
    return str(v) if v is not None else 'None'


def _prompt_multiselect(label, choices, defaults):
    """Print choices, let user pick by number.  Enter keeps *defaults*.

    Returns a list of selected values from *choices*.
    """
    # Map display labels
    labels = [_format_value(c) for c in choices]
    default_indices = [i for i, c in enumerate(choices) if c in defaults]

    print(f'\n  {label}')
    for i, lab in enumerate(labels):
        marker = '*' if i in default_indices else ' '
        print(f'    [{marker}] {i + 1}) {lab}')

    default_str = ','.join(str(i + 1) for i in default_indices)
    sel = input(f'  Select (comma-separated, default={default_str}): ').strip()

    if sel == '':
        return [choices[i] for i in default_indices]

    try:
        indices = [int(x.strip()) - 1 for x in sel.split(',')]
        return [choices[i] for i in indices if 0 <= i < len(choices)]
    except (ValueError, IndexError):
        print('  Invalid input — keeping defaults.')
        return [choices[i] for i in default_indices]


def _prompt_single(label, choices, default):
    """Single-value numbered prompt.  Enter keeps *default*."""
    labels = [_format_value(c) for c in choices]
    default_idx = choices.index(default) if default in choices else 0

    print(f'\n  {label}')
    for i, lab in enumerate(labels):
        marker = '*' if i == default_idx else ' '
        print(f'    [{marker}] {i + 1}) {lab}')

    sel = input(f'  Select (default={default_idx + 1}): ').strip()

    if sel == '':
        return default

    try:
        idx = int(sel) - 1
        return choices[idx] if 0 <= idx < len(choices) else default
    except (ValueError, IndexError):
        print('  Invalid input — keeping default.')
        return default


def _prompt_int(label, default):
    """Prompt for an integer value."""
    sel = input(f'\n  {label} (default={default}): ').strip()
    if sel == '':
        return default
    try:
        return int(sel)
    except ValueError:
        print('  Invalid input — keeping default.')
        return default


def _interactive_grid(args, default_co_cuts=None, default_ndbd_cuts=None,
                      default_co_shifts=None):
    """Run interactive prompts for the three grid axes.

    Only called when -i is set.  Flags already provided on the CLI
    take priority (those axes are skipped in the prompt).
    """
    dcc = default_co_cuts if default_co_cuts is not None else CO_CUTS
    dnc = default_ndbd_cuts if default_ndbd_cuts is not None else NDBD_CUTS
    dcs = default_co_shifts if default_co_shifts is not None else CO_SHIFTS

    if args.co_cut is None:
        args.co_cut = _prompt_multiselect('co_cuts', ALL_CO_CUTS, dcc)
    if args.ndbd_cutoff is None:
        args.ndbd_cutoff = _prompt_multiselect('ndbd_cutoffs', ALL_NDBD_CUTS, dnc)
    if args.co_shift is None:
        args.co_shift = _prompt_multiselect(
            'co_shifts', ALL_CO_SHIFTS, dcs)

    return args


# --------------------------------------------------------------------------- #
#  Grid builders                                                               #
# --------------------------------------------------------------------------- #

def _resolve_grid(args, default_co_cuts=None, default_ndbd_cuts=None,
                  default_co_shifts=None):
    """Resolve the three grid axes from flags/interactive/defaults."""
    dcc = default_co_cuts if default_co_cuts is not None else CO_CUTS
    dnc = default_ndbd_cuts if default_ndbd_cuts is not None else NDBD_CUTS
    dcs = default_co_shifts if default_co_shifts is not None else CO_SHIFTS

    co_cuts      = args.co_cut   if args.co_cut   is not None else dcc
    ndbd_cutoffs = args.ndbd_cutoff if args.ndbd_cutoff is not None else dnc
    co_shifts    = args.co_shift if args.co_shift  is not None else dcs

    return co_cuts, ndbd_cutoffs, co_shifts


def build_toy_grid(co_cuts, ndbd_cutoffs, co_shifts, fit_types,
                   n_reps, master_seed, environment, toys_config,
                   analysis=None, model_key='default'):
    rows, task_id = [], 1
    rng = random.Random(master_seed)

    for co_cut, ndbd_cut, co_shift, fit_type, rep in itertools.product(
        co_cuts, ndbd_cutoffs, co_shifts, fit_types, range(n_reps),
    ):
        seed = rng.randint(0, 2**31 - 1)
        cmd = (
            f'python {TOY_SCRIPT}'
            f' --toys {toys_config}'
            f' --environment {environment}'
        )
        if analysis is not None:
            cmd += f' --analysis {analysis}'
        if model_key:
            cmd += f' --model-key {model_key}'
        cmd += (
            f' --type {fit_type}'
            f' --n 1'
            f' --co-cut {co_cut}'
            f' --ndbd-cutoff {ndbd_cut}'
            f' --co-shift {co_shift}'
            f' --seed {seed}'
        )
        rows.append({
            'task_id': task_id, 'seed': seed,
            'co_cut': co_cut, 'ndbd_cutoff': ndbd_cut,
            'co_shift': co_shift, 'fit_type': fit_type,
            'rep': rep, 'command': cmd,
        })
        task_id += 1

    return rows


def build_fixed_toy_grid(co_cuts, ndbd_cutoffs, co_shifts, fit_type,
                         n_fits, master_seed, environment, toys_config,
                         analysis=None, model_key='default'):
    """Build one row per model point for a fixed toy-like fit type."""
    rows, task_id = [], 1
    rng = random.Random(master_seed)

    for co_cut, ndbd_cut, co_shift in itertools.product(
        co_cuts, ndbd_cutoffs, co_shifts,
    ):
        seed = rng.randint(0, 2**31 - 1)
        cmd = (
            f'python {TOY_SCRIPT}'
            f' --toys {toys_config}'
            f' --environment {environment}'
        )
        if analysis is not None:
            cmd += f' --analysis {analysis}'
        if model_key:
            cmd += f' --model-key {model_key}'
        cmd += (
            f' --type {fit_type}'
            f' --n {n_fits}'
            f' --co-cut {co_cut}'
            f' --ndbd-cutoff {ndbd_cut}'
            f' --co-shift {co_shift}'
            f' --seed {seed}'
        )
        rows.append({
            'task_id': task_id, 'seed': seed,
            'co_cut': co_cut, 'ndbd_cutoff': ndbd_cut,
            'co_shift': co_shift, 'fit_type': fit_type,
            'rep': 0, 'command': cmd,
        })
        task_id += 1

    return rows


def build_fit_grid(co_cuts, ndbd_cutoffs, co_shifts,
                   environment, chains, draws, analysis=None, prior_width='5sig',
                   model_key='default'):
    rows, task_id = [], 1

    for co_cut, ndbd_cut, co_shift in itertools.product(
        co_cuts, ndbd_cutoffs, co_shifts,
    ):
        cmd = (
            f'python {FIT_SCRIPT}'
            f' --environment {environment}'
        )
        if analysis is not None:
            cmd += f' --analysis {analysis}'
        if model_key:
            cmd += f' --model-key {model_key}'
        cmd += (
            f' --co-cut {co_cut}'
            f' --ndbd-cutoff {ndbd_cut}'
            f' --co-shift {co_shift}'
            f' --chains {chains}'
            f' --draws {draws}'
            f' --prior-width {prior_width}'
        )
        rows.append({
            'task_id': task_id, 'seed': 0,
            'co_cut': co_cut, 'ndbd_cutoff': ndbd_cut,
            'co_shift': co_shift, 'fit_type': 'fit',
            'rep': 0, 'command': cmd,
        })
        task_id += 1

    return rows


def build_input_grid(co_cuts, ndbd_cutoffs, co_shifts, environment, analysis,
                     model_key='default'):
    """One task per (co_cut, ndbd_cutoff). All co_shifts handled in-process."""
    rows, task_id = [], 1
    shifts_str = ' '.join(str(s) for s in co_shifts)

    for co_cut, ndbd_cut in itertools.product(
        co_cuts, ndbd_cutoffs,
    ):
        cmd = (
            f'python {INPUT_SCRIPT}'
            f' --environment {environment}'
        )
        if analysis is not None:
            cmd += f' --analysis {analysis}'
        if model_key:
            cmd += f' --model-key {model_key}'
        cmd += (
            f' --co-cut {co_cut}'
            f' --ndbd-cutoff {ndbd_cut}'
            f' --co-shifts {shifts_str}'
        )
        rows.append({
            'task_id': task_id,
            'co_cut': co_cut, 'ndbd_cutoff': ndbd_cut,
            'co_shifts': shifts_str, 'command': cmd,
        })
        task_id += 1

    return rows


# --------------------------------------------------------------------------- #
#  Writers                                                                     #
# --------------------------------------------------------------------------- #

def _write_tsv(rows, path, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='\n') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def _write_tsv_append(rows, path, fieldnames):
    """Append rows to an existing TSV, re-numbering task_id sequentially."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        # Read existing rows to find max task_id
        with open(path, 'r', newline='\n') as f:
            reader = csv.DictReader(f, delimiter='\t')
            existing = list(reader)
        max_id = max((int(r['task_id']) for r in existing), default=0)
        for i, row in enumerate(rows, start=max_id + 1):
            row['task_id'] = i
        with open(path, 'a', newline='\n') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writerows(rows)
    else:
        _write_tsv(rows, path, fieldnames)


def _print_summary(subcommand, rows, out):
    print(f'Subcommand : {subcommand}')
    print(f'Generated  : {len(rows)} tasks')
    print(f'Output     : {out}')
    if rows:
        print(f'Array range: 1-{len(rows)}')


# --------------------------------------------------------------------------- #
#  Subcommand: staging                                                         #
# --------------------------------------------------------------------------- #

def add_staging_parser(subparsers):
    p = subparsers.add_parser(
        'staging',
        help='Generate staging params (ROOT → staging ROOT)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--environment', '-e', default='tech_note')
    p.add_argument('--analysis', '-a', default=None)
    p.add_argument('--ndbd-only', action='store_true')
    p.set_defaults(func=run_staging)


def run_staging(args):
    out = args.out or os.path.join(PARAMS_DIR, 'staging_params.tsv')
    os.makedirs(os.path.dirname(out), exist_ok=True)

    command = [
        sys.executable, STAGING_SCRIPT,
        '--environment', args.environment,
        '--create-jobslist',
        '--jobslist', out,
    ]
    if args.analysis is not None:
        command.extend(['--analysis', args.analysis])
    if args.ndbd_only:
        command.append('--ndbd-only')

    subprocess.check_call(command)

    # Count rows for summary
    with open(out) as f:
        n_rows = sum(1 for _ in f) - 1
    print('Subcommand : staging')
    print(f'Generated  : {n_rows} tasks')
    print(f'Output     : {out}')
    if n_rows > 0:
        print(f'Array range: 1-{n_rows}')


# --------------------------------------------------------------------------- #
#  Subcommand: input                                                           #
# --------------------------------------------------------------------------- #

def add_input_parser(subparsers):
    p = subparsers.add_parser(
        'input',
        help='Generate input-creation params (staging ROOT → HDF5)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--environment', '-e', default='tech_note')
    p.add_argument('--analysis', '-a', default=None)
    p.add_argument('--model-key', '-m', default='default')
    _add_grid_args(p)
    p.set_defaults(func=run_input)


def run_input(args):
    if args.interactive:
        _interactive_grid(args,
                          default_co_cuts=ALL_CO_CUTS,
                          default_ndbd_cuts=ALL_NDBD_CUTS,
                          default_co_shifts=ALL_CO_SHIFTS)

    co_cuts, ndbd_cutoffs, co_shifts = _resolve_grid(
        args,
        default_co_cuts=ALL_CO_CUTS,
        default_ndbd_cuts=ALL_NDBD_CUTS,
        default_co_shifts=ALL_CO_SHIFTS,
    )

    out = args.out or os.path.join(PARAMS_DIR, 'input_params.tsv')
    rows = build_input_grid(co_cuts, ndbd_cutoffs, co_shifts,
                            args.environment, args.analysis,
                            model_key=args.model_key)
    _write_tsv(rows, out, FIELDNAMES_INPUT)
    _print_summary('input', rows, out)


# --------------------------------------------------------------------------- #
#  Subcommand: fit                                                             #
# --------------------------------------------------------------------------- #

def add_fit_parser(subparsers):
    p = subparsers.add_parser(
        'fit',
        help='Generate main Bayesian-fit params',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--environment', '-e', default='tech_note')
    p.add_argument('--analysis', '-a', default=None)
    p.add_argument('--model-key', '-m', default='default')
    p.add_argument('--chains', type=int, default=8)
    p.add_argument('--draws', type=int, default=20000)
    p.add_argument('--prior-width', choices=['5sig', '10sig', 'syst', '2syst'], default='5sig')
    _add_grid_args(p)
    p.set_defaults(func=run_fit)


def run_fit(args):
    if args.interactive:
        _interactive_grid(args)
        args.chains = _prompt_int('chains', args.chains)
        args.draws  = _prompt_int('draws',  args.draws)

    co_cuts, ndbd_cutoffs, co_shifts = _resolve_grid(args)

    out = args.out or os.path.join(PARAMS_DIR, 'fit_params.tsv')
    rows = build_fit_grid(co_cuts, ndbd_cutoffs, co_shifts,
                          args.environment, args.chains, args.draws,
                          analysis=args.analysis, prior_width=args.prior_width,
                          model_key=args.model_key)
    _write_tsv(rows, out, FIELDNAMES_GRID)
    _print_summary('fit', rows, out)


# --------------------------------------------------------------------------- #
#  Subcommands: toys / sens / bias                                             #
# --------------------------------------------------------------------------- #

def add_toylike_parser(subparsers, run_type):
    p = subparsers.add_parser(
        run_type,
        help=f'Generate {run_type} batch-fit params',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--environment', '-e', default='tech_note')
    p.add_argument('--analysis', '-a', default=None)
    p.add_argument('--model-key', '-m', default='default')
    p.add_argument('--n-fits', type=int, default=200,
                   help='Number of fits to run per model-point task (default: 200)')
    p.add_argument('--toys-config', default='batch',
                   help='Toy config name from config/toys/ (default: batch)')
    p.add_argument('--seed', type=int, default=MASTER_SEED,
                   help=f'Master seed (default: {MASTER_SEED})')
    _add_grid_args(p)
    p.set_defaults(func=run_toylike, run_type=run_type)


def run_toylike(args):
    if args.interactive:
        _interactive_grid(args)
        args.n_fits = _prompt_int('n_fits', args.n_fits)

    co_cuts, ndbd_cutoffs, co_shifts = _resolve_grid(args)

    out = args.out or os.path.join(PARAMS_DIR, f'{args.run_type}_params.tsv')
    rows = build_fixed_toy_grid(co_cuts, ndbd_cutoffs, co_shifts,
                                args.run_type, args.n_fits, args.seed,
                                args.environment, args.toys_config,
                                analysis=args.analysis,
                                model_key=args.model_key)
    _write_tsv(rows, out, FIELDNAMES_GRID)
    _print_summary(args.run_type, rows, out)


# --------------------------------------------------------------------------- #
#  Subcommand: syst_lss                                                        #
# --------------------------------------------------------------------------- #

# Fixed model for the LSS systematic study
_SYST_LSS_CO_CUT     = 'gamma10_v150'
_SYST_LSS_NDBD_CUT   = 'u1740'
_SYST_LSS_DRAWS      = 4000
_SYST_LSS_CHAINS     = 2
_SYST_LSS_MKSAMPLES  = 100
_SYST_LSS_OUTPUT_DIR = 'data/syst_results'

# Fixed model for the MC (BM-stat) systematic study
_SYST_MC_CO_CUT      = 'gamma10_v150'
_SYST_MC_NDBD_CUT    = 'u1740'
_SYST_MC_DRAWS       = 4000
_SYST_MC_CHAINS      = 2
_SYST_MC_MKSAMPLES   = 100
_SYST_MC_N           = 200
_SYST_MC_OUTPUT_DIR  = 'data/syst_results'


def add_syst_mc_parser(subparsers):
    p = subparsers.add_parser(
        'syst_mc',
        help='Generate MC (BM-stat) systematic fit params (one task per co_shift)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--environment', '-e', default='m2v5-1')
    p.add_argument('--analysis', '-a', default='tech_note')
    p.add_argument('--model-key', '-m', default='default')
    p.add_argument('--co-shift', nargs='+', type=float, default=None,
                   help=f'co_shift values (default: all {ALL_CO_SHIFTS})')
    p.add_argument('--n', type=int, default=_SYST_MC_N,
                   help=f'BM-stat samples per task (default: {_SYST_MC_N})')
    p.add_argument('--draws', type=int, default=_SYST_MC_DRAWS)
    p.add_argument('--chains', type=int, default=_SYST_MC_CHAINS)
    p.add_argument('--mkchain-samples', type=int, default=_SYST_MC_MKSAMPLES)
    p.add_argument('--output-dir', default=_SYST_MC_OUTPUT_DIR,
                   help='Directory for per-task output NetCDF files')
    p.add_argument('--njobs', type=int, default=0,
                   help='Parallel workers for the runner (0 = auto)')
    p.add_argument('--bias-csv', default=None,
                   help='Path to batch CSV with bias fit results')
    p.add_argument('--xmin-bias', type=float, default=15,
                   help='Minimum ndbd_inj for bias fit (default: 15)')
    p.add_argument('--save-traces', action='store_true', default=False,
                   help='Save per-sample MCMC traces to <syst_dir>/<model_dir>/trace_mc{idx:04d}.nc')
    p.add_argument('--stat-only', action='store_true', default=False,
                   help='Sample raw JAGS chain without covariance perturbation (stat-only)')
    p.set_defaults(func=run_syst_mc)


def run_syst_mc(args):
    co_shifts = args.co_shift if args.co_shift is not None else ALL_CO_SHIFTS

    rows = []
    for task_id, co_shift in enumerate(co_shifts, start=1):
        co_shift_tag = str(co_shift).replace('.', 'p')
        out_path = os.path.join(args.output_dir, f'mc_results_cs{co_shift_tag}.nc')

        cmd = (
            f'python -m qpym2.run.run_syst_mc'
            f' --environment {args.environment}'
        )
        if args.analysis is not None:
            cmd += f' --analysis {args.analysis}'
        if getattr(args, 'model_key', None):
            cmd += f' --model-key {args.model_key}'
        cmd += (
            f' --co-cut {_SYST_MC_CO_CUT}'
            f' --ndbd-cutoff {_SYST_MC_NDBD_CUT}'
            f' --co-shift {co_shift}'
            f' --n {args.n}'
            f' --draws {args.draws}'
            f' --chains {args.chains}'
            f' --mkchain-samples {args.mkchain_samples}'
            f' --output {out_path}'
        )
        if args.njobs > 0:
            cmd += f' --njobs {args.njobs}'
        if args.bias_csv is not None:
            cmd += f' --bias-csv {args.bias_csv}'
            cmd += f' --xmin-bias {args.xmin_bias}'
        if args.save_traces:
            cmd += ' --save-traces'
        if args.stat_only:
            cmd += ' --stat-only'
        rows.append({
            'task_id': task_id,
            'co_shift': co_shift,
            'command': cmd,
        })

    out = args.out or os.path.join(PARAMS_DIR, 'syst_mc_params.tsv')
    _write_tsv(rows, out, FIELDNAMES_SYST_MC)
    _print_summary('syst_mc', rows, out)


def add_syst_lss_parser(subparsers):
    p = subparsers.add_parser(
        'syst_lss',
        help='Generate LSS systematic fit params (one task per co_shift)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--syst-type', required=True, choices=['sscl', 'dq'],
                   help='LSS systematic type: sscl (sigma-scaling) or dq (Q-residual)')
    p.add_argument('--environment', '-e', default='m2v5-1')
    p.add_argument('--analysis', '-a', default='tech_note')
    p.add_argument('--model-key', '-m', default='default')
    p.add_argument('--lss-data-path', required=True,
                   help='Path to per-dataset LSS ROOT files')
    p.add_argument('--co-shift', nargs='+', type=float, default=None,
                   help=f'co_shift values (default: all {ALL_CO_SHIFTS})')
    p.add_argument('--draws', type=int, default=_SYST_LSS_DRAWS)
    p.add_argument('--chains', type=int, default=_SYST_LSS_CHAINS)
    p.add_argument('--mkchain-samples', type=int, default=_SYST_LSS_MKSAMPLES)
    p.add_argument('--njobs', type=int, default=0,
                   help='Parallel workers for the runner (0 = auto)')
    p.add_argument('--output-dir', default=_SYST_LSS_OUTPUT_DIR,
                   help='Directory for per-task output NetCDF files')
    p.add_argument('--bias-csv', default=None,
                   help='Path to batch CSV with bias fit results')
    p.add_argument('--xmin-bias', type=float, default=15,
                   help='Minimum ndbd_inj for bias fit (default: 15)')
    # dq-specific args
    p.add_argument('--n-dq', type=int, default=100,
                   help='Number of random dq parameter draws (default: 100)')
    p.add_argument('--dq-seed', type=int, default=None,
                   help='RNG seed for dq sampling reproducibility')
    p.add_argument('--save-traces', action='store_true', default=False,
                   help='Save per-dataset MCMC traces to output-dir/traces_<syst_type>/')
    p.add_argument('--pm1', action='store_true', default=False,
                   help='dq only: use fixed ±1σ parameter sets instead of random draws')
    p.add_argument('--append', action='store_true', default=False,
                   help='Append rows to existing TSV instead of overwriting')
    p.set_defaults(func=run_syst_lss)


def run_syst_lss(args):
    co_shifts = args.co_shift if args.co_shift is not None else ALL_CO_SHIFTS
    syst_type = args.syst_type

    rows = []
    for task_id, co_shift in enumerate(co_shifts, start=1):
        cmd = (
            f'python -m qpym2.run.run_syst_lss {syst_type}'
            f' --environment {args.environment}'
        )
        if args.analysis is not None:
            cmd += f' --analysis {args.analysis}'
        if getattr(args, 'model_key', None):
            cmd += f' --model-key {args.model_key}'
        cmd += (
            f' --lss-data-path {args.lss_data_path}'
            f' --co-cut {_SYST_LSS_CO_CUT}'
            f' --ndbd-cutoff {_SYST_LSS_NDBD_CUT}'
            f' --co-shift {co_shift}'
            f' --draws {args.draws}'
            f' --chains {args.chains}'
            f' --mkchain-samples {args.mkchain_samples}'
            f' --output-dir {args.output_dir}'
        )
        if args.njobs > 0:
            cmd += f' --njobs {args.njobs}'
        if syst_type == 'dq':
            if getattr(args, 'pm1', False):
                cmd += ' --pm1'
            else:
                cmd += f' --n-dq {args.n_dq}'
                if args.dq_seed is not None:
                    cmd += f' --dq-seed {args.dq_seed}'
        if args.bias_csv is not None:
            cmd += f' --bias-csv {args.bias_csv}'
            cmd += f' --xmin-bias {args.xmin_bias}'
        if args.save_traces:
            cmd += ' --save-traces'
        rows.append({
            'task_id': task_id,
            'co_shift': co_shift,
            'command': cmd,
        })

    out = args.out or os.path.join(PARAMS_DIR, 'syst_lss_params.tsv')
    if getattr(args, 'append', False):
        _write_tsv_append(rows, out, FIELDNAMES_SYST_LSS)
    else:
        _write_tsv(rows, out, FIELDNAMES_SYST_LSS)
    _print_summary(f'syst_lss ({syst_type})', rows, out)


# --------------------------------------------------------------------------- #
#  Shared arg helpers                                                          #
# --------------------------------------------------------------------------- #

def _add_grid_args(parser):
    """Add --co-cut, --ndbd-cutoff, --co-shift multi-value flags."""
    parser.add_argument(
        '--co-cut', nargs='+', default=None,
        help=f'co_cut values (default: {[_format_value(c) for c in CO_CUTS]})')
    parser.add_argument(
        '--ndbd-cutoff', nargs='+', default=None,
        help=f'ndbd_cutoff values (default: {[_format_value(c) for c in NDBD_CUTS]})')
    parser.add_argument(
        '--co-shift', nargs='+', type=float, default=None,
        help=f'co_shift values (default: {CO_SHIFTS})')


# --------------------------------------------------------------------------- #
#  Main CLI                                                                    #
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        prog='gen_params',
        description=textwrap.dedent('''\
            Unified parameter-grid generator for QPyM2 run types.

            Subcommands:
              staging   ROOT → staging ROOT (delegates to run_create_staging.py)
              input     staging ROOT → fit input HDF5
              fit       main Bayesian fits
              toys      likelihood-validation toy fits
              sens      background-only sensitivity toy fits
              bias      injected-signal bias toy fits
              syst_lss  LSS systematic fits (one task per co_shift)
              syst_mc   MC (BM-stat) systematic fits (one task per co_shift)
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-i', '--interactive', action='store_true',
                        help='Prompt interactively for grid values')
    parser.add_argument('--out', default=None,
                        help='Override output TSV path')

    subparsers = parser.add_subparsers(dest='subcommand', required=True)
    add_staging_parser(subparsers)
    add_input_parser(subparsers)
    add_fit_parser(subparsers)
    add_toylike_parser(subparsers, 'toys')
    add_toylike_parser(subparsers, 'sens')
    add_toylike_parser(subparsers, 'bias')
    add_syst_lss_parser(subparsers)
    add_syst_mc_parser(subparsers)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
