"""
Single-task toy fit runner for interactive or SLURM-driven runs.

Each array task runs n toy fits for one fit type and appends
the resulting rows to a shared CSV file (with file locking).

Usage:
    python -m qpym2.run.run_toy_batch \
        --type sens --n 5 --output-csv /path/to/results.csv

    If running under SLURM, SLURM_ARRAY_TASK_ID is used as a seed
    offset and written to output for traceability.
"""

import argparse
import contextlib
import os
import sys
import numpy as np
import pandas as pd

from qpym2.cfg.hydra_adapter import get_config
from qpym2.cfg.config_builders import get_fit_dir, get_h2uv_model
from qpym2.cfg.toy_config import load_toy_config
from qpym2.utils.logging import start_timer, log_elapsed
from qpym2.io.csv import task_csv_path, flatten_columns
from qpym2.run.run_toy_fits import (
    run_toy_fit, run_sens_fit, run_bias_fit, _build_model_grid,
)

ENV_TASK_ID_NAME = "SLURM_ARRAY_TASK_ID"


class Tee:
    """Write stream output to multiple file-like objects."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def main():
    parser = argparse.ArgumentParser(
        description="Single-task toy fit runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python -m qpym2.run.run_toy_batch --type sens --n 5\n"
            "    Run sensitivity fits using default config (base.yaml).\n\n"
            "  python -m qpym2.run.run_toy_batch --type bias --toys batch\n"
            "    Run bias fits using batch.yaml config.\n\n"
            "  python -m qpym2.run.run_toy_batch --type bias --bias-min 10 --bias-max 500\n"
            "    Override config defaults with CLI args.\n"
        ),
    )
    parser.add_argument("--toys", default=None,
                        help="Toy config name from config/toys/ (default: base)")
    parser.add_argument("--type", "-t", required=True,
                        choices=["toys", "sens", "bias"],
                        help="Fit type to run")
    parser.add_argument("--n", type=int, default=5,
                        help="Number of toy samples this task runs")
    parser.add_argument("--output-csv", "-o", default=None,
                        help="Path to shared output CSV (default from config)")
    parser.add_argument("--environment", "-e", default=None,
                        help="Hydra environment name (default from config: local)")
    parser.add_argument("--analysis", default=None,
                        help="Analysis profile (default from config: base)")
    parser.add_argument("--model-key", "-m", default="default",
                        help="Model key in config/models/ (default: default)")
    parser.add_argument("--bias-min", type=float, default=None,
                        help="Min injected signal (default from config)")
    parser.add_argument("--bias-max", type=float, default=None,
                        help="Max injected signal (default from config)")
    parser.add_argument("--write", "-w", action="store_true", default=False,
                        help="Also write detailed HDF5 output via io.toys")
    # Model grid params (single point by default, matching notebook)
    parser.add_argument("--co-cut", default="gamma10_v150")
    parser.add_argument("--ndbd-cutoff", default="u1740")
    parser.add_argument("--co-shift", type=float, default=1.0)
    parser.add_argument("--seed", "-s", type=int, default=0,
                        help="Seed for reproducibility, default to task ID or 0")
    parser.add_argument("--log-file", default=None,
                        help="Optional filepath to also write stdout/stderr logs (default from config)")
    parser.add_argument("--no-append", action="store_true", default=False,
                        help="Fail if output files exist (default: append mode from config)")
    args = parser.parse_args()
    
    # Load toy config for defaults
    toy_cfg = load_toy_config(args.toys)
    
    # Apply config defaults where args weren't provided
    if args.output_csv is None:
        args.output_csv = toy_cfg['output']['csv'].format(type=args.type)
    if args.bias_min is None:
        args.bias_min = toy_cfg['bias']['min']
    if args.bias_max is None:
        args.bias_max = toy_cfg['bias']['max']
    # Note: args.log_file remains None unless explicitly provided via --log-file
    # Shell scripts handle logging by default; Python logging is opt-in
    
    # Determine append mode (CLI --no-append overrides config)
    append_mode = not args.no_append if args.no_append else toy_cfg['output']['append']
    
    # Check for existing files if not in append mode
    if not append_mode:
        if os.path.exists(args.output_csv):
            raise FileExistsError(
                f"Output CSV already exists and append=False: {args.output_csv}\n"
                f"Remove the file or use append mode (default) to continue."
            )
        if args.log_file and os.path.exists(args.log_file):
            raise FileExistsError(
                f"Log file already exists and append=False: {args.log_file}\n"
                f"Remove the file or use append mode (default) to continue."
            )

    def _run_once():
        task_id = int(os.environ.get(ENV_TASK_ID_NAME, 0))

        # Seed for reproducibility (each task gets a different seed)
        seed = args.seed or task_id
        np.random.seed(47 + seed)

        # Load config and build model
        cfg = get_config(environment=args.environment, analysis=args.analysis)
        fit_dir = get_fit_dir(cfg)
        h2uv_model = get_h2uv_model(cfg, model_key=args.model_key)

        fit_models = _build_model_grid(
            cfg, h2uv_model,
            co_cuts=[args.co_cut],
            ndbd_cutoffs=[args.ndbd_cutoff],
            co_shifts=[args.co_shift],
        )
        mcfg = fit_models.cfg.values.flatten()[0]

        print(f"[task {task_id}] running {args.n} {args.type} fits "
              f"for model {mcfg['name']}  (write={args.write})")

        # Run fits
        if args.type == "toys":
            df = run_toy_fit(mcfg, fit_dir, n=args.n, write=args.write)
        elif args.type == "sens":
            df = run_sens_fit(mcfg, fit_dir, n=args.n, write=args.write)
        elif args.type == "bias":
            df = run_bias_fit(mcfg, fit_dir, n=args.n,
                              bias_min=args.bias_min,
                              bias_max=args.bias_max, write=args.write)

        # Flatten any tuple columns (bias fits) for CSV compatibility
        df = flatten_columns(df)
        df["task_id"] = task_id
        df["co_cut"] = args.co_cut
        df["ndbd_cutoff"] = args.ndbd_cutoff
        df["co_shift"] = args.co_shift
        df.drop(columns=["model_"], errors="ignore", inplace=True)

        # Write per-task CSV (scatter; merge later with scripts/collect_results.py)
        out_path = task_csv_path(args.output_csv, task_id)
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        df.to_csv(out_path, index=False)
        print(f"[task {task_id}] wrote {len(df)} rows to {out_path}")

    if args.log_file:
        log_dir = os.path.dirname(args.log_file)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        with open(args.log_file, "a", buffering=1) as logf:
            tee_out = Tee(sys.stdout, logf)
            tee_err = Tee(sys.stderr, logf)
            with contextlib.redirect_stdout(tee_out), contextlib.redirect_stderr(tee_err):
                print(f"[logging] writing stdout/stderr to {args.log_file}")
                full_start = start_timer()
                try:
                    _run_once()
                finally:
                    log_elapsed("run_toy_batch total", full_start)
    else:
        full_start = start_timer()
        try:
            _run_once()
        finally:
            log_elapsed("run_toy_batch total", full_start)


if __name__ == "__main__":
    main()
