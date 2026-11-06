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
from omegaconf import OmegaConf

from qpym2.studies import assemble_study


def parse_args():
    parser = argparse.ArgumentParser(
        description="QPyM2 Study Assembly & Input Verification Runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--study-dir",
        "--study-name",
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
        "--check-fits",
        action="store_true",
        default=False,
        help="Check and print fit completion status for all study points",
    )
    parser.add_argument(
        "--run-fits",
        action="store_true",
        default=False,
        help="Execute MCMC Bayesian fits across all study points",
    )
    parser.add_argument(
        "--force-fits",
        action="store_true",
        default=False,
        help="Overwrite existing fits even if cached",
    )
    parser.add_argument(
        "--save-trace",
        action="store_true",
        default=False,
        help="Persist full ArviZ posterior traces (.nc)",
    )
    parser.add_argument(
        "--collect-results",
        action="store_true",
        default=False,
        help="Collect completed fit metrics into results.csv and SQLite",
    )
    # Bias study arguments
    parser.add_argument(
        "--check-bias",
        action="store_true",
        default=False,
        help="Check and print bias toy completion status for all study points",
    )
    parser.add_argument(
        "--run-bias",
        action="store_true",
        default=False,
        help="Execute bias toy fits across all study points",
    )
    parser.add_argument(
        "--collect-bias",
        action="store_true",
        default=False,
        help="Merge bias task CSVs and fit calibration curves via get_bias_from_fitdf",
    )
    # Sensitivity study arguments
    parser.add_argument(
        "--check-sens",
        action="store_true",
        default=False,
        help="Check and print sensitivity toy completion status for all study points",
    )
    parser.add_argument(
        "--run-sens",
        action="store_true",
        default=False,
        help="Execute sensitivity toy fits (H0: ndbd=0) across all study points",
    )
    parser.add_argument(
        "--collect-sens",
        action="store_true",
        default=False,
        help="Merge sensitivity task CSVs, compute median sensitivity and quantiles, and save sens_summary.json",
    )
    parser.add_argument(
        "--sens-bias-calib",
        default=None,
        help="Optional path to bias_calibration.json to compute bias-corrected sensitivities",
    )
    # Systematic study arguments
    parser.add_argument(
        "--check-syst",
        action="store_true",
        default=False,
        help="Check and print systematic task completion status for all study points",
    )
    parser.add_argument(
        "--run-syst",
        action="store_true",
        default=False,
        help="Execute systematic fits across all study points",
    )
    parser.add_argument(
        "--collect-syst",
        action="store_true",
        default=False,
        help="Merge systematic task results, combine mixture posteriors, and compute delta shifts",
    )
    parser.add_argument(
        "--category",
        default=None,
        help="Systematic category override (e.g. bmmc, lss_sscl, lss_dq)",
    )
    parser.add_argument(
        "--dq-mode",
        default=None,
        choices=["pm1", "random"],
        help="Evaluation mode for lss_dq (pm1 or random)",
    )
    parser.add_argument(
        "--n-dq",
        type=int,
        default=None,
        help="Number of random parameter draws for lss_dq random mode",
    )
    parser.add_argument(
        "--lss-data-path",
        default=None,
        help="Path to LSS ROOT calibration files directory",
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=None,
        help="Number of systematic samples to evaluate (default: from config or 100)",
    )
    parser.add_argument(
        "--stat-only",
        action="store_true",
        default=False,
        help="Sample raw JAGS chain without covariance matrix perturbation (stat-only)",
    )
    parser.add_argument(
        "--n-toys",
        type=int,
        default=None,
        help="Number of toys to run per study point (default: from config or 200)",
    )
    parser.add_argument(
        "--bias-min",
        type=float,
        default=None,
        help="Minimum injected signal count (default: from config or 5.0)",
    )
    parser.add_argument(
        "--bias-max",
        type=float,
        default=None,
        help="Maximum injected signal count (default: from config or 200.0)",
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=15.0,
        help="Minimum n_inj threshold for linear calibration fit (default: 15.0)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Master RNG seed for bias toys (default: 42)",
    )
    parser.add_argument(
        "--prior-width",
        choices=["5sig", "10sig", "syst", "2syst"],
        default="5sig",
        help="Prior width mode",
    )
    parser.add_argument(
        "--chains",
        type=int,
        default=8,
        help="Number of MCMC chains",
    )
    parser.add_argument(
        "--draws",
        type=int,
        default=20000,
        help="Number of MCMC draws per chain",
    )
    parser.add_argument(
        "-o",
        "--override",
        dest="overrides",
        action="append",
        default=None,
        help="Dotlist configuration override (e.g. -o reference.co_shift=0.5)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
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

    if getattr(args, "category", None):
        if not hasattr(study.cfg, "syst") or study.cfg.syst is None:
            study.cfg.syst = OmegaConf.create({})
        study.cfg.syst.category = args.category
    if getattr(args, "dq_mode", None):
        if not hasattr(study.cfg, "syst") or study.cfg.syst is None:
            study.cfg.syst = OmegaConf.create({})
        study.cfg.syst.dq_mode = args.dq_mode
    if getattr(args, "n_dq", None):
        if not hasattr(study.cfg, "syst") or study.cfg.syst is None:
            study.cfg.syst = OmegaConf.create({})
        study.cfg.syst.n_samples = args.n_dq
    if getattr(args, "lss_data_path", None):
        if not hasattr(study.cfg, "syst") or study.cfg.syst is None:
            study.cfg.syst = OmegaConf.create({})
        study.cfg.syst.lss_data_path = args.lss_data_path

    # 1. Manifest generation if requested
    if args.save_manifest:
        mpath = study.save_manifest()
        print(f"[INFO] Saved study manifest: {mpath}")

    # 2. Input generation if requested
    if args.generate_inputs:
        print(f"[INFO] Generating inputs (overwrite={args.overwrite_inputs})...")
        study.generate_inputs(overwrite=args.overwrite_inputs, verbose=args.verbose)

    # 3. Print input status report
    is_bias_mode = args.check_bias or args.run_bias or args.collect_bias
    is_sens_mode = args.check_sens or args.run_sens or args.collect_sens
    is_syst_mode = args.check_syst or args.run_syst or args.collect_syst
    if args.check_inputs or (
        not args.check_fits
        and not args.run_fits
        and not args.collect_results
        and not is_bias_mode
        and not is_sens_mode
        and not is_syst_mode
    ):
        study.print_input_status()

    # 4. Standard fit execution if requested
    if args.run_fits:
        print(f"[INFO] Running fits across study points (force={args.force_fits})...")
        study.run_fits(
            save_trace=args.save_trace,
            force=args.force_fits,
            chains=args.chains,
            draws=args.draws,
            prior_width=args.prior_width,
            verbose=args.verbose,
        )

    # 5. Check fit completion status if requested
    if args.check_fits or args.run_fits:
        study.print_fit_status()

    # 6. Standard result collection into results.csv
    if args.collect_results or args.run_fits:
        print("[INFO] Collecting results into results.csv...")
        res_df = study.collect_results(
            prior_width=args.prior_width,
            chains=args.chains,
            draws=args.draws,
        )
        print(f"[INFO] Result summary:\n{res_df}")

    # 7. Bias toy execution if requested
    if args.run_bias:
        print(f"[INFO] Running bias toys across study points (force={args.force_fits})...")
        study.run_bias_toys(
            n_toys=args.n_toys,
            bias_min=args.bias_min,
            bias_max=args.bias_max,
            seed=args.seed,
            chains=args.chains if args.chains != 8 else 4,
            draws=args.draws if args.draws != 20000 else 2000,
            force=args.force_fits,
            verbose=args.verbose,
        )

    # 8. Check bias toy completion status
    if args.check_bias or args.run_bias:
        study.print_bias_status(n_toys=args.n_toys)

    # 9. Collect bias results and fit calibration curves
    if args.collect_bias or args.run_bias:
        print("[INFO] Collecting bias results and fitting calibration curves...")
        study.collect_bias_results(xmin=args.xmin, verbose=args.verbose)

    # 10. Sensitivity toy execution if requested
    if args.run_sens:
        print(f"[INFO] Running sensitivity toys across study points (force={args.force_fits})...")
        study.run_sens_toys(
            n_toys=args.n_toys,
            seed=args.seed,
            chains=args.chains if args.chains != 8 else 4,
            draws=args.draws if args.draws != 20000 else 2000,
            force=args.force_fits,
            verbose=args.verbose,
        )

    # 11. Check sensitivity toy completion status
    if args.check_sens or args.run_sens:
        study.print_sens_status(n_toys=args.n_toys)

    # 12. Collect sensitivity results and compute summary metrics
    if args.collect_sens or args.run_sens:
        print("[INFO] Collecting sensitivity results and computing summary metrics...")
        study.collect_sens_results(bias_calibration=args.sens_bias_calib, verbose=args.verbose)

    # 13. Systematic fit execution if requested
    if args.run_syst:
        print(f"[INFO] Running systematic fits across study points (force={args.force_fits})...")
        study.run_syst(
            n_samples=args.n_samples,
            seed=args.seed,
            chains=args.chains if args.chains != 8 else 2,
            draws=args.draws if args.draws != 20000 else 2500,
            stat_only=args.stat_only,
            force=args.force_fits,
            save_traces=args.save_trace,
            verbose=args.verbose,
        )

    # 14. Check systematic completion status
    if args.check_syst or args.run_syst:
        study.print_syst_status(n_samples=args.n_samples)

    # 15. Collect systematic results, compute delta shifts, and update SQLite
    if args.collect_syst or args.run_syst:
        print("[INFO] Collecting systematic results and computing delta shifts...")
        study.collect_syst_results(bias_calibration=args.sens_bias_calib, verbose=args.verbose)


if __name__ == "__main__":
    main()
