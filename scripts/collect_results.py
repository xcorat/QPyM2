#!/usr/bin/env python3
"""Merge per-task CSV scatter files into a single results CSV.

Run this once after a SLURM array of toy/sens/bias fits completes.

Usage:
    python scripts/collect_results.py --type toys
    python scripts/collect_results.py --type sens --config batch
    python scripts/collect_results.py --type bias --output custom_out.csv
    python scripts/collect_results.py --type bias --base data_2700/toy_results/toys_batch_bias.csv
    python scripts/collect_results.py --type bias --base data_2700/toy_results/toys_batch_bias.csv --output data_2700/toy_results/toys_batch_bias.csv
"""

import argparse
import os
import sys

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from qpym2.cfg.toy_config import load_toy_config
from qpym2.io.csv import collect_task_csvs


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-task toy result CSVs into a single file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python scripts/collect_results.py --type toys\n"
            "  python scripts/collect_results.py --type bias --config batch\n"
            "  python scripts/collect_results.py --type sens --output my_sens.csv\n"
            "  python scripts/collect_results.py --type bias --base data_2700/toy_results/toys_batch_bias.csv\n"
        ),
    )
    parser.add_argument("--type", "-t", required=True,
                        choices=["toys", "sens", "bias"],
                        help="Fit type (resolves base CSV path from config)")
    parser.add_argument("--config", "-c", default="batch",
                        help="Toy config name from config/toys/ (default: batch)")
    parser.add_argument("--base", "-b", default=None,
                        help="Override scatter file glob base path "
                             "(default: derived from --config and --type). "
                             "Use this for non-default analysis folders, e.g. data_2700/.")
    parser.add_argument("--output", "-o", default=None,
                        help="Override merged output path "
                             "(default: same as --base or config-derived path)")
    args = parser.parse_args()

    if args.base:
        base_path = args.base
    else:
        toy_cfg = load_toy_config(args.config)
        base_path = toy_cfg['output']['csv'].format(type=args.type)

    output_path = args.output or base_path

    print(f"Collecting per-task CSVs for base path: {base_path}")
    df = collect_task_csvs(base_path, output_path=output_path)
    print(f"Merged {len(df)} rows → {output_path}")


if __name__ == "__main__":
    main()
