"""
Standalone CLI runner for QPyM2 studies (v0.6).

Assembles studies from study directories, validates input prerequisites,
and triggers on-demand input generation.

Usage:
    # Check input readiness for a study
    python -m qpym2.run.run_study --study-dir studies/example_co_shift

    # Save study manifest JSON
    python -m qpym2.run.run_study --study-dir studies/example_co_shift --save-manifest

    # Generate missing inputs across all points in the study
    python -m qpym2.run.run_study --study-dir studies/example_co_shift --generate-inputs
"""

import argparse
import sys

from qpym2.studies import assemble_study


def parse_args():
    parser = argparse.ArgumentParser(
        description="QPyM2 Study Assembly & Input Verification Runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-s", "--study-dir", "--study-name",
        dest="study_dir",
        required=True,
        help="Path or name of the study folder containing study.yaml",
    )
    parser.add_argument(
        "--check-inputs",
        action="store_true",
        default=False,
        help="Check and print input readiness table for all study points",
    )
    parser.add_argument(
        "--generate-inputs",
        action="store_true",
        default=False,
        help="Materialize any missing input templates and data",
    )
    parser.add_argument(
        "--overwrite-inputs",
        action="store_true",
        default=False,
        help="Force regeneration and overwrite of inputs even if cached",
    )
    parser.add_argument(
        "--save-manifest",
        action="store_true",
        default=False,
        help="Save <study_name>.manifest JSON file to the study directory",
    )
    parser.add_argument(
        "-o", "--override",
        dest="overrides",
        action="append",
        default=None,
        help="Dotlist configuration override (e.g. -o reference.co_shift=0.5)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Verbose logging output",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"\n[INFO] Assembling study from '{args.study_dir}'...")
    try:
        study = assemble_study(args.study_dir, overrides=args.overrides)
    except Exception as e:
        print(f"[ERROR] Failed to assemble study: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Successfully assembled {study}.")

    # 1. Manifest generation if requested
    if args.save_manifest:
        mpath = study.save_manifest()
        print(f"[INFO] Saved study manifest: {mpath}")

    # 2. Input generation if requested
    if args.generate_inputs:
        print(f"[INFO] Generating inputs (overwrite={args.overwrite_inputs})...")
        study.generate_inputs(overwrite=args.overwrite_inputs, verbose=args.verbose)

    # 3. Always print input status report (unless only save_manifest was explicitly asked)
    study.print_input_status()


if __name__ == "__main__":
    main()
