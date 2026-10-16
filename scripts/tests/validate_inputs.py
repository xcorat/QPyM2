#!/usr/bin/env python3
"""Validate input HDF5 files produced by run_create_input.py.

Reads the param grid from input_params.tsv and checks that every expected
HDF5 file exists, is not corrupted, and contains the required groups/datasets.

Usage:
    python scripts/test/validate_inputs.py
    python scripts/test/validate_inputs.py --params scripts/params/input_params.tsv
    python scripts/test/validate_inputs.py --verbose
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import h5py
import pandas as pd

# Ensure the project root is importable
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from qpym2.cfg.model_utils import build_model_dir, get_shift_suffix


# Expected HDF5 structure from hdf.write_input
REQUIRED_H5_GROUPS = ['bkg_model', 'mchists', 'data']
REQUIRED_SIGNAL_NAME = 'ndbd'


def parse_args():
    parser = argparse.ArgumentParser(description='Validate input HDF5 files.')
    parser.add_argument('--params', default='scripts/params/input_params.tsv',
                        help='Path to input_params.tsv')
    parser.add_argument('--fit-dir',
                        default=str(_PROJECT_ROOT / 'data' / 'tech_note'),
                        help='Base output directory (fit_dir)')
    parser.add_argument('--verbose', '-v', action='store_true')
    return parser.parse_args()


def read_param_grid(tsv_path):
    """Read input_params.tsv and return list of (co_cut, ndbd_cutoff, [co_shifts])."""
    tasks = []
    with open(tsv_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            co_cut = row['co_cut'].strip() or None
            ndbd_cutoff = row['ndbd_cutoff'].strip() or None
            shifts = [float(s) for s in row['co_shifts'].strip().split()]
            tasks.append((co_cut, ndbd_cutoff, shifts))
    return tasks


def _norm_cut(val):
    """Normalize None-like strings to actual None."""
    if val in (None, 'None', 'null', ''):
        return None
    return val


def expected_path(fit_dir, co_cut, ndbd_cutoff, co_shift):
    """Build the expected HDF5 path for a single (co_cut, ndbd_cutoff, co_shift)."""
    cc = _norm_cut(co_cut)
    nc = _norm_cut(ndbd_cutoff)
    dir_name = build_model_dir(base_name='model', co_cut=cc, ndbd_cutoff=nc)
    suffix = get_shift_suffix(co_shift)
    fname = f'input{suffix}.h5'
    return Path(fit_dir) / dir_name / fname


def validate_h5(fpath, verbose=False):
    """Check that an HDF5 file is valid and contains expected structure.

    Returns (ok: bool, issues: list[str]).
    """
    issues = []

    if not fpath.exists():
        return False, ['MISSING']

    if fpath.stat().st_size == 0:
        return False, ['EMPTY FILE (0 bytes)']

    # Check h5py-level structure
    try:
        with h5py.File(fpath, 'r') as h5f:
            keys = list(h5f.keys())

            for g in REQUIRED_H5_GROUPS:
                if g not in keys:
                    issues.append(f'missing group: {g}')

            if 'mchists' in keys:
                mchists = h5f['mchists']
                n_mc = len(mchists.keys())
                if n_mc == 0:
                    issues.append('mchists group is empty')
                elif REQUIRED_SIGNAL_NAME not in mchists:
                    issues.append(f'mchists missing "{REQUIRED_SIGNAL_NAME}" dataset')
                if verbose:
                    print(f'    mchists: {n_mc} datasets: {list(mchists.keys())[:5]}...')

            if 'data' in keys:
                data = h5f['data']
                if data.shape[0] == 0:
                    issues.append('data dataset has 0 entries')
                if verbose:
                    print(f'    data shape: {data.shape}, sum: {data[:].sum():.1f}')

    except Exception as e:
        issues.append(f'HDF5 READ ERROR: {e}')

    # Check pandas-stored tables (bkg_model, ndbd)
    try:
        bm = pd.read_hdf(fpath, key='bkg_model')
        if bm.empty:
            issues.append('bkg_model table is empty')
        if verbose:
            print(f'    bkg_model: {len(bm)} components')
    except Exception as e:
        if 'missing group: bkg_model' not in issues:
            issues.append(f'bkg_model read error: {e}')

    try:
        sig = pd.read_hdf(fpath, key=REQUIRED_SIGNAL_NAME)
        if sig.empty:
            issues.append('ndbd table is empty')
    except Exception as e:
        issues.append(f'ndbd table read error: {e}')

    return len(issues) == 0, issues


def file_meta(fpath):
    """Return mtime and size for an existing file, or None."""
    if not fpath.exists():
        return None
    stat = fpath.stat()
    return {
        'mtime': datetime.fromtimestamp(stat.st_mtime),
        'size_mb': stat.st_size / (1024 * 1024),
    }


def main():
    args = parse_args()

    if not os.path.exists(args.params):
        print(f'ERROR: params file not found: {args.params}')
        sys.exit(1)

    tasks = read_param_grid(args.params)
    print(f'Param grid: {len(tasks)} tasks from {args.params}')
    print(f'Fit dir:    {args.fit_dir}')
    print()

    total = 0
    ok_count = 0
    fail_count = 0
    missing_count = 0
    results = []

    for co_cut, ndbd_cutoff, co_shifts in tasks:
        label = build_model_dir('model',
                                _norm_cut(co_cut),
                                _norm_cut(ndbd_cutoff))
        for co_shift in co_shifts:
            total += 1
            fpath = expected_path(args.fit_dir, co_cut, ndbd_cutoff, co_shift)
            suffix = get_shift_suffix(co_shift)
            file_label = f'{label}/input{suffix}.h5'

            if args.verbose:
                print(f'  Checking {file_label}')

            valid, issues = validate_h5(fpath, verbose=args.verbose)
            meta = file_meta(fpath)

            if valid:
                ok_count += 1
                status = 'OK'
                if meta and args.verbose:
                    print(f'    size={meta["size_mb"]:.1f}MB  mtime={meta["mtime"]}')
            elif 'MISSING' in issues:
                missing_count += 1
                status = 'MISSING'
            else:
                fail_count += 1
                status = 'FAILED'

            results.append({
                'file': file_label,
                'status': status,
                'issues': issues,
                'mtime': meta['mtime'] if meta else None,
                'size_mb': meta['size_mb'] if meta else None,
            })

    # ---- Summary -----------------------------------------------------------
    print('=' * 70)
    print(f'RESULTS: {total} files expected')
    print(f'  OK:      {ok_count}')
    print(f'  MISSING: {missing_count}')
    print(f'  FAILED:  {fail_count}')
    print('=' * 70)

    # Print failures/missing
    bad = [r for r in results if r['status'] != 'OK']
    if bad:
        print('\nIssues:')
        for r in bad:
            print(f"  [{r['status']}] {r['file']}")
            for issue in r['issues']:
                print(f'         {issue}')

    # Check for suspicious modification times (files written at the same time
    # across different directories = sign of concurrent array-bug overwrites)
    mtimes = [r for r in results if r['mtime'] is not None]
    if mtimes:
        print('\nModification time summary:')
        times_sorted = sorted(mtimes, key=lambda r: r['mtime'])
        first = times_sorted[0]['mtime']
        last = times_sorted[-1]['mtime']
        print(f'  Earliest: {first}')
        print(f'  Latest:   {last}')

        # Group by minute-bucket to detect concurrent writes
        buckets = defaultdict(list)
        for r in times_sorted:
            bucket_key = r['mtime'].strftime('%Y-%m-%d %H:%M')
            buckets[bucket_key].append(r['file'])

        # >7 files in one minute means multiple directories were written
        # simultaneously (one task only produces 7 files for its co_shifts)
        multi_write = {k: v for k, v in buckets.items() if len(v) > 7}
        if multi_write:
            print(f'\n  WARNING: {len(multi_write)} time buckets with >7 files '
                  '(possible concurrent overwrites from array bug):')
            for t, files in sorted(multi_write.items()):
                print(f'    {t}: {len(files)} files')

    sys.exit(1 if fail_count > 0 else 0)


if __name__ == '__main__':
    main()
