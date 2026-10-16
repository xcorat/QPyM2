#!/usr/bin/env python3
"""Run all commands in a params TSV sequentially in the current terminal.

Usage:
    python scripts/run_params_seq.py scripts/params/fit_params.tsv
    python scripts/run_params_seq.py scripts/params/fit_params.tsv --ids 1 3 5
    python scripts/run_params_seq.py scripts/params/fit_params.tsv --dry-run
"""
import argparse
import subprocess
import sys
import time
import csv


def parse_args():
    p = argparse.ArgumentParser(description='Run params TSV commands sequentially.')
    p.add_argument('tsv', help='Path to params TSV file')
    p.add_argument('--ids', nargs='+', type=int, metavar='ID',
                   help='Only run these task_ids (default: all)')
    p.add_argument('--dry-run', action='store_true',
                   help='Print commands without running them')
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    if args.ids:
        id_set = set(args.ids)
        rows = [r for r in rows if int(r['task_id']) in id_set]
        if not rows:
            print(f"No rows found for ids {args.ids}", file=sys.stderr)
            sys.exit(1)

    print(f"Running {len(rows)} task(s) from {args.tsv}")

    failed = []
    for row in rows:
        tid = row['task_id']
        cmd = row['command'].strip()
        print(f"\n[task {tid}] {cmd}")
        if args.dry_run:
            continue

        t0 = time.time()
        result = subprocess.run(cmd, shell=True)
        elapsed = time.time() - t0

        if result.returncode != 0:
            print(f"[task {tid}] FAILED (exit {result.returncode}) after {elapsed:.0f}s")
            failed.append(tid)
        else:
            print(f"[task {tid}] done in {elapsed:.0f}s")

    if failed:
        print(f"\n{len(failed)} task(s) failed: {failed}", file=sys.stderr)
        sys.exit(1)
    elif not args.dry_run:
        print(f"\nAll tasks completed successfully.")


if __name__ == '__main__':
    main()
