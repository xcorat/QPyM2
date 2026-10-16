#!/usr/bin/env python3
"""Systematic-study parameter-grid generator for QPyM2.

Subcommands:  syst_lss, syst_mc, ndbd_shift
Each writes a TSV to scripts/params/<subcommand>_params.tsv by default.

Moved from syst_lss / syst_mc sections of gen_params.py to keep that
file focused on the primary pipeline grid generation.

Usage:
    python scripts/gen_syst_params.py syst_lss --syst-type sscl --lss-data-path /path/to/lss
    python scripts/gen_syst_params.py syst_mc --n 200
    python scripts/gen_syst_params.py ndbd_shift
    python scripts/gen_syst_params.py ndbd_shift --co-shift 0.5 1.0
"""

import argparse
import csv
import os
import sys
import textwrap

# Ensure the project root is importable (scripts/ lives one level below)
_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from qpym2.cfg.params_grid import ALL_CO_SHIFTS


# --------------------------------------------------------------------------- #
#  Constants                                                                   #
# --------------------------------------------------------------------------- #
PARAMS_DIR = os.path.join('scripts', 'params')

SYST_LSS_SCRIPT     = 'qpym2.run.run_syst_lss'
SYST_MC_SCRIPT      = 'qpym2.run.run_syst_mc'
NDBD_SHIFT_SCRIPT   = 'qpym2.run.run_ndbd_shift_syst'

FIELDNAMES_SYST = [
    'task_id', 'co_shift', 'command',
]

# Default systematic model settings
_SYST_CO_CUT       = 'gamma10_v150'
_SYST_NDBD_CUT     = 'u1740'
_SYST_DRAWS        = 4000
_SYST_CHAINS       = 2
_SYST_MKSAMPLES    = 100
_SYST_OUTPUT_DIR   = 'data/syst_results'
_SYST_MC_N         = 200


# --------------------------------------------------------------------------- #
#  TSV helpers                                                                 #
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
#  Subcommand: syst_lss                                                        #
# --------------------------------------------------------------------------- #

def add_syst_lss_parser(subparsers):
    p = subparsers.add_parser(
        'syst_lss',
        help='Generate LSS systematic fit params (one task per co_shift)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--syst-type', required=True, choices=['sscl', 'dq'],
                   help='LSS systematic type: sscl (sigma-scaling) or dq (Q-residual)')
    p.add_argument('--environment', '-e', default='tech_note')
    p.add_argument('--analysis', '-a', default=None)
    p.add_argument('--lss-data-path', required=True,
                   help='Path to per-dataset LSS ROOT files')
    p.add_argument('--co-shift', nargs='+', type=float, default=None,
                   help=f'co_shift values (default: all {ALL_CO_SHIFTS})')
    p.add_argument('--draws', type=int, default=_SYST_DRAWS)
    p.add_argument('--chains', type=int, default=_SYST_CHAINS)
    p.add_argument('--mkchain-samples', type=int, default=_SYST_MKSAMPLES)
    p.add_argument('--njobs', type=int, default=0,
                   help='Parallel workers for the runner (0 = auto)')
    p.add_argument('--output-dir', default=_SYST_OUTPUT_DIR,
                   help='Directory for per-task output NetCDF files')
    p.add_argument('--bias-csv', default=None,
                   help='Path to batch CSV with bias fit results')
    p.add_argument('--xmin-bias', type=float, default=15,
                   help='Minimum ndbd_inj for bias fit (default: 15)')
    p.add_argument('--n-dq', type=int, default=100,
                   help='Number of random dq parameter draws (default: 100)')
    p.add_argument('--dq-seed', type=int, default=None,
                   help='RNG seed for dq sampling reproducibility')
    p.add_argument('--save-traces', action='store_true', default=False,
                   help='Save per-dataset MCMC traces')
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
            f'python -m {SYST_LSS_SCRIPT} {syst_type}'
            f' --environment {args.environment}'
        )
        if args.analysis is not None:
            cmd += f' --analysis {args.analysis}'
        cmd += (
            f' --lss-data-path {args.lss_data_path}'
            f' --co-cut {_SYST_CO_CUT}'
            f' --ndbd-cutoff {_SYST_NDBD_CUT}'
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

    out = args.out or os.path.join(PARAMS_DIR, f'syst_lss_params.tsv')
    if getattr(args, 'append', False):
        _write_tsv_append(rows, out, FIELDNAMES_SYST)
    else:
        _write_tsv(rows, out, FIELDNAMES_SYST)
    _print_summary(f'syst_lss ({syst_type})', rows, out)


# --------------------------------------------------------------------------- #
#  Subcommand: syst_mc                                                         #
# --------------------------------------------------------------------------- #

def add_syst_mc_parser(subparsers):
    p = subparsers.add_parser(
        'syst_mc',
        help='Generate MC (BM-stat) systematic fit params (one task per co_shift)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--environment', '-e', default='tech_note')
    p.add_argument('--analysis', '-a', default=None)
    p.add_argument('--co-shift', nargs='+', type=float, default=None,
                   help=f'co_shift values (default: all {ALL_CO_SHIFTS})')
    p.add_argument('--n', type=int, default=_SYST_MC_N,
                   help=f'BM-stat samples per task (default: {_SYST_MC_N})')
    p.add_argument('--draws', type=int, default=_SYST_DRAWS)
    p.add_argument('--chains', type=int, default=_SYST_CHAINS)
    p.add_argument('--mkchain-samples', type=int, default=_SYST_MKSAMPLES)
    p.add_argument('--output-dir', default=_SYST_OUTPUT_DIR,
                   help='Directory for per-task output NetCDF files')
    p.add_argument('--njobs', type=int, default=0,
                   help='Parallel workers for the runner (0 = auto)')
    p.add_argument('--bias-csv', default=None,
                   help='Path to batch CSV with bias fit results')
    p.add_argument('--xmin-bias', type=float, default=15,
                   help='Minimum ndbd_inj for bias fit (default: 15)')
    p.add_argument('--save-traces', action='store_true', default=False,
                   help='Save per-sample MCMC traces')
    p.add_argument('--stat-only', action='store_true', default=False,
                   help='Sample raw JAGS chain without covariance perturbation')
    p.set_defaults(func=run_syst_mc)


def run_syst_mc(args):
    co_shifts = args.co_shift if args.co_shift is not None else ALL_CO_SHIFTS

    rows = []
    for task_id, co_shift in enumerate(co_shifts, start=1):
        co_shift_tag = str(co_shift).replace('.', 'p')
        out_path = os.path.join(args.output_dir, f'mc_results_cs{co_shift_tag}.nc')

        cmd = (
            f'python -m {SYST_MC_SCRIPT}'
            f' --environment {args.environment}'
        )
        if args.analysis is not None:
            cmd += f' --analysis {args.analysis}'
        cmd += (
            f' --co-cut {_SYST_CO_CUT}'
            f' --ndbd-cutoff {_SYST_NDBD_CUT}'
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
    _write_tsv(rows, out, FIELDNAMES_SYST)
    _print_summary('syst_mc', rows, out)


# --------------------------------------------------------------------------- #
#  Subcommand: ndbd_shift                                                      #
# --------------------------------------------------------------------------- #

def add_ndbd_shift_parser(subparsers):
    p = subparsers.add_parser(
        'ndbd_shift',
        help='Generate 0νββ signal energy-shift systematic params (single task)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--environment', '-e', default='tech_note')
    p.add_argument('--analysis', '-a', default=None)
    p.add_argument('--draws', type=int, default=_SYST_DRAWS)
    p.add_argument('--chains', type=int, default=_SYST_CHAINS)
    p.add_argument('--output-dir', default=_SYST_OUTPUT_DIR,
                   help='Output directory for results')
    p.add_argument('--njobs', type=int, default=0,
                   help='Parallel workers for the runner (0 = auto)')
    p.add_argument('--bias-csv', default=None,
                   help='Path to batch CSV with bias fit results')
    p.add_argument('--xmin-bias', type=float, default=15,
                   help='Minimum ndbd_inj for bias fit (default: 15)')
    p.add_argument('--save-traces', action='store_true', default=True,
                   help='Save per-shift MCMC traces (default: True)')
    p.add_argument('--no-save-traces', action='store_false', dest='save_traces',
                   help='Do not save per-shift MCMC traces')
    p.set_defaults(func=run_ndbd_shift)


def run_ndbd_shift(args):
    rows = []
    cmd = (
        f'python -m {NDBD_SHIFT_SCRIPT}'
        f' --environment {args.environment}'
    )
    if args.analysis is not None:
        cmd += f' --analysis {args.analysis}'
    cmd += (
        f' --co-cut {_SYST_CO_CUT}'
        f' --ndbd-cutoff {_SYST_NDBD_CUT}'
        f' --draws {args.draws}'
        f' --chains {args.chains}'
        f' --output-dir {args.output_dir}'
    )
    if args.njobs > 0:
        cmd += f' --njobs {args.njobs}'
    if args.bias_csv is not None:
        cmd += f' --bias-csv {args.bias_csv}'
        cmd += f' --xmin-bias {args.xmin_bias}'
    if args.save_traces:
        cmd += ' --save-traces'
    else:
        cmd += ' --no-save-traces'

    rows.append({
        'task_id': 1,
        'co_shift': 'all',
        'command': cmd,
    })

    out = args.out or os.path.join(PARAMS_DIR, 'syst_ndbd_shift_params.tsv')
    _write_tsv(rows, out, FIELDNAMES_SYST)
    _print_summary('ndbd_shift', rows, out)


# --------------------------------------------------------------------------- #
#  Main CLI                                                                    #
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        prog='gen_syst_params',
        description=textwrap.dedent('''\
            Systematic-study parameter-grid generator for QPyM2.

            Subcommands:
              syst_lss    LSS systematic fits (one task per co_shift)
              syst_mc     MC (BM-stat) systematic fits (one task per co_shift)
              ndbd_shift  0νββ signal energy-shift systematic (single task)
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--out', default=None,
                        help='Override output TSV path')

    subparsers = parser.add_subparsers(dest='subcommand', required=True)
    add_syst_lss_parser(subparsers)
    add_syst_mc_parser(subparsers)
    add_ndbd_shift_parser(subparsers)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
