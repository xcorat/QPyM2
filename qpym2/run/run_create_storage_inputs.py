"""
Standalone runner to generate storage inputs from staging ROOT files.

Creates identity-addressed source-level templates and observed data under:
fit_data/roi_<min>_<max>/thr<val>/<coord>_b<bin_size>/input/

Usage:
    # Dry-run preview
    python -m qpym2.run.run_create_storage_inputs --co-shifts 0.8 1.0 1.2 --dry-run

    # Generate single point with multiple shifts
    python -m qpym2.run.run_create_storage_inputs \
        -e local -a base -m tech_note_200 \
        --co-cut gamma10_v150 --ndbd-cutoff u1740 --co-shifts 0.0 0.8 1.0 1.2
"""

import argparse
import sys
import time
from typing import List, Optional

from qpym2.cfg.hydra_adapter_new import HydraConfigManager
from qpym2.storage.inputs import create_storage_inputs_for_point
from qpym2.utils.logging import log, debug


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate identity-addressed storage inputs from staging ROOT files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Dry-run inspection
  python -m qpym2.run.run_create_storage_inputs --dry-run

  # Run for a specific model key and shifts
  python -m qpym2.run.run_create_storage_inputs \\
      --model-key tech_note_200 \\
      --co-cut gamma10_v150 --ndbd-cutoff u1740 \\
      --co-shifts 0.8 1.0 1.2
        """,
    )

    parser.add_argument(
        '--environment', '-e',
        default='local',
        help='Hydra environment config name (default: local)',
    )
    parser.add_argument(
        '--analysis', '-a',
        default='base',
        help='Analysis profile name (default: base)',
    )
    parser.add_argument(
        '--model-key', '-m',
        default='tech_note_200',
        help='Model key in config/models/ (default: tech_note_200)',
    )
    parser.add_argument(
        '--co-cut',
        default='gamma10_v150',
        help='Co-60 cut name (default: gamma10_v150)',
    )
    parser.add_argument(
        '--ndbd-cutoff',
        default='u1740',
        help='NDBD cutoff name (default: u1740)',
    )
    parser.add_argument(
        '--co-shifts',
        nargs='+',
        type=float,
        default=None,
        help='One or more Co-60 energy shifts in keV (e.g. 0.8 1.0 1.2). '
             'Invariant sources are cached and reused across shifts.',
    )
    parser.add_argument(
        '--blinded',
        action='store_true',
        help='Generate blinded observed data (data_blinded.h5)',
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing templates even if cached',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Preview destination paths and UUIDs without reading or writing files',
    )

    return parser.parse_args()


def main():
    """Execute storage input generation."""
    args = parse_args()

    print(f"Loading configuration (environment={args.environment}, analysis={args.analysis})...")
    cfg = HydraConfigManager.get_config(
        environment=args.environment,
        analysis=args.analysis,
    )

    co_shifts: List[float] = args.co_shifts if args.co_shifts is not None else [0.0]

    print("\n" + "=" * 70)
    print("QPyM2 Storage Input Generation (v0.6)")
    print("=" * 70)
    print(f"  Model key    : {args.model_key}")
    print(f"  Co-60 cut    : {args.co_cut}")
    print(f"  NDBD cutoff  : {args.ndbd_cutoff}")
    print(f"  Co shifts    : {co_shifts}")
    print(f"  Blinded      : {args.blinded}")
    print(f"  Overwrite    : {args.overwrite}")
    print(f"  Dry run      : {args.dry_run}")
    print("=" * 70 + "\n")

    hist_cache = {}
    total_written = 0
    total_cached = 0
    t_start = time.perf_counter()

    for shift in co_shifts:
        print(f"--- Processing co_shift = {shift:.2f} keV ---")
        t0 = time.perf_counter()
        manifest = create_storage_inputs_for_point(
            cfg=cfg,
            model_key=args.model_key,
            co_cut=args.co_cut,
            ndbd_cutoff=args.ndbd_cutoff,
            co_shift=shift,
            blinded=args.blinded,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            hist_cache=hist_cache,
        )
        elapsed = time.perf_counter() - t0

        input_dir = manifest['input_dir']
        templates = manifest['templates']
        n_cached = sum(1 for t in templates.values() if t['cached'])
        n_written = len(templates) - n_cached

        total_written += n_written
        total_cached += n_cached

        print(f"  Input directory : {input_dir}")
        print(f"  Templates       : {len(templates)} total ({n_cached} cached, {n_written} written) in {elapsed:.2f}s")
        data_info = manifest['data']
        d_status = "cached" if data_info.get('cached') else ("dry-run" if args.dry_run else "written")
        print(f"  Observed data   : {data_info.get('path')} [{d_status}]")

        if args.dry_run:
            print("  Sample templates (first 5):")
            for src, info in list(templates.items())[:5]:
                print(f"    - {src:12s} -> {info['uuid']} ({'cached' if info['cached'] else 'pending'})")

    total_elapsed = time.perf_counter() - t_start
    print("\n" + "=" * 70)
    print(f"Completed in {total_elapsed:.2f}s")
    print(f"Total: {total_written} templates written, {total_cached} cache hits.")
    print("=" * 70 + "\n")


if __name__ == '__main__':
    main()
