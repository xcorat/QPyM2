"""
Atomic PyMC Bayesian fit runner for QPyM2 storage architecture (v0.6).

Runs a single Bayesian fit for one parameter point against the storage layout:
- Reads source-level unsmoothed HDF5 templates from fit_data/<coords>/input/
- Computes canonical deterministic <output_UUID> from parameter set and MCMC options
- Checks cache existence to skip redundant MCMC sampling
- Persists summary fit_<output_UUID>.h5 and optional trace fit_<output_UUID>.nc
- Optionally indexes the point into a study SQLite database

Usage:
    # Dry run preview of output UUID and cache status
    python -m qpym2.run.run_single_storage_fit -e local -a base -m tech_note_200 --co-shift 0.8 --dry-run

    # Execute fit with 8 chains, 20000 draws
    python -m qpym2.run.run_single_storage_fit -e local -a base -m tech_note_200 --co-shift 0.8

    # Retain full posterior trace (.nc)
    python -m qpym2.run.run_single_storage_fit -e local -a base -m tech_note_200 --co-shift 0.8 --save-trace
"""

import argparse
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

from qpym2.cfg.hydra_adapter import HydraConfigManager
from qpym2.fit import fit as fit_module
from qpym2.storage.backend import read_fit, write_fit
from qpym2.storage.fits import setup_storage_fit
from qpym2.storage.study import record_study_point


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a single Bayesian MCMC fit in QPyM2 storage architecture.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-e", "--environment", default="local",
        help="Hydra environment config name (e.g. local, nersc)",
    )
    parser.add_argument(
        "-a", "--analysis", default="base",
        help="Analysis profile name (e.g. base)",
    )
    parser.add_argument(
        "-m", "--model", default="default",
        help="Model key from models tier (e.g. default, tech_note_200)",
    )
    parser.add_argument(
        "--config-dir", default=None,
        help="Optional path to custom config directory",
    )

    # Parameter point
    parser.add_argument(
        "--co-cut", default=None,
        help="Co-60 cut key or expression (e.g. gamma10_v150)",
    )
    parser.add_argument(
        "--ndbd-cutoff", default=None,
        help="NDBD cutoff key or expression (e.g. u1740)",
    )
    parser.add_argument(
        "--co-shift", type=float, default=0.0,
        help="Co-60 energy scale shift in keV",
    )
    parser.add_argument(
        "--ndbd-shift", type=float, default=0.0,
        help="NDBD energy scale shift in keV",
    )
    parser.add_argument(
        "--blinded", action="store_true", default=False,
        help="Fit blinded observed data",
    )
    parser.add_argument(
        "--skip-roi", nargs=2, type=float, metavar=("MIN", "MAX"), default=None,
        help="Exclude energy range from fit (e.g. --skip-roi 2515 2540)",
    )

    # Priors and MCMC
    parser.add_argument(
        "--prior-width", choices=["5sig", "10sig", "syst", "2syst"], default="5sig",
        help="Prior width mode",
    )
    parser.add_argument(
        "--mkchain-samples", type=int, default=100,
        help="Number of JAGS chain samples for prior generation",
    )
    parser.add_argument(
        "--chains", type=int, default=8,
        help="Number of MCMC chains",
    )
    parser.add_argument(
        "--draws", type=int, default=20000,
        help="Number of MCMC draws per chain",
    )
    parser.add_argument(
        "--tune", type=int, default=1000,
        help="Number of tuning steps",
    )
    parser.add_argument(
        "--single-core", action="store_true", default=False,
        help="Run MCMC with cores=1 (useful in packed Slurm jobs)",
    )

    # Storage and Execution
    parser.add_argument(
        "--save-trace", action="store_true", default=False,
        help="Persist ArviZ InferenceData posterior trace (.nc)",
    )
    parser.add_argument(
        "--study-dir", default=None,
        help="Optional study directory to register completed point in SQLite DB",
    )
    parser.add_argument(
        "--force", action="store_true", default=False,
        help="Overwrite existing fit results even if cached",
    )
    parser.add_argument(
        "--dry-run", action="store_true", default=False,
        help="Inspect parameters, output UUID, and cache status without sampling",
    )

    return parser.parse_args()


def run_single_storage_fit(
    cfg: Any,
    model_key: str = "tech_note_200",
    co_cut: Optional[str] = None,
    ndbd_cutoff: Optional[str] = None,
    co_shift: float = 0.0,
    ndbd_shift: float = 0.0,
    blinded: bool = False,
    skip_roi: Optional[Tuple[float, float]] = None,
    prior_width: str = "5sig",
    mkchain_samples: int = 100,
    chains: int = 8,
    draws: int = 20000,
    tune: int = 1000,
    single_core: bool = False,
    save_trace: bool = False,
    study_dir: Optional[Union[str, Path]] = None,
    force: bool = False,
    dry_run: bool = False,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Execute or reuse a single Bayesian fit within the QPyM2 storage architecture.

    Returns:
        Dictionary with fit setup, execution status, output paths, and posterior results.
    """
    # 1. Setup fit objects, parameter hashes, and target paths
    fit_setup = setup_storage_fit(
        cfg=cfg,
        model_key=model_key,
        co_cut=co_cut,
        ndbd_cutoff=ndbd_cutoff,
        co_shift=co_shift,
        ndbd_shift=ndbd_shift,
        blinded=blinded,
        skip_roi=skip_roi,
        prior_width=prior_width,
        mkchain_samples=mkchain_samples,
        chains=chains,
        draws=draws,
        tune=tune,
        dry_run=dry_run,
    )

    output_uuid = fit_setup["output_uuid"]
    fit_prefix = fit_setup["fit_prefix"]
    h5_path = fit_setup["h5_path"]
    nc_path = fit_setup["nc_path"]
    is_cached = fit_setup["is_cached"]

    if verbose:
        print(f"\n[INFO] Storage Fit Point: model={model_key}, co_shift={co_shift}, co_cut={co_cut}")
        print(f"       Parameter-set UUIDs: base={fit_setup['uuid_base']}, co={fit_setup['uuid_co']}, ndbd={fit_setup['uuid_ndbd']}")
        print(f"       Output UUID: {output_uuid}")
        print(f"       Target: {fit_prefix}.h5")
        print(f"       Cached on disk: {is_cached}")

    # Dry run
    if dry_run:
        return {
            "status": "dry_run",
            "output_uuid": output_uuid,
            "fit_prefix": fit_prefix,
            "h5_path": h5_path,
            "is_cached": is_cached,
            "setup": fit_setup,
        }

    # 2. Check cache
    if is_cached and not force:
        if verbose:
            print(f"[CACHE HIT] Reusing existing fit results from: {h5_path}")
        table, trace, data, metadata = read_fit(fit_prefix)
        return {
            "status": "cached",
            "output_uuid": output_uuid,
            "fit_prefix": fit_prefix,
            "h5_path": h5_path,
            "nc_path": nc_path if nc_path.exists() else None,
            "table": table,
            "trace": trace,
            "data": data,
            "metadata": metadata,
        }

    # 3. Perform Bayesian MCMC sampling
    if verbose:
        print(f"[SAMPLING] Running PyMC sampling: chains={chains}, draws={draws}, single_core={single_core}...")

    trace = fit_module.run_fit(
        fit_setup["pymc_model"],
        single_core=single_core,
        chains=chains,
        draws=draws,
        tune=tune,
    )

    # 4. Save results to disk
    written_h5, written_nc = write_fit(
        path_prefix=fit_prefix,
        table=fit_setup["fit_comps"],
        trace=trace,
        data=fit_setup["data"],
        metadata=fit_setup["fit_params"],
        save_trace=save_trace,
        overwrite=force,
    )

    if verbose:
        print(f"[SUCCESS] Wrote fit components table: {written_h5}")
        if written_nc:
            print(f"[SUCCESS] Wrote ArviZ trace: {written_nc}")

    # 5. Optional SQLite registration in study DB
    if study_dir is not None:
        sdir = Path(study_dir)
        db_path = sdir / f"{sdir.name}.sqlite"
        if db_path.exists():
            point_record = {
                "point_id": f"shift_{co_shift}",
                "model_key": model_key,
                "co_shift": co_shift,
                "ndbd_shift": ndbd_shift,
                "co_cut": str(co_cut),
                "ndbd_cutoff": str(ndbd_cutoff),
                "output_uuid": output_uuid,
                "status": "completed",
            }
            record_study_point(db_path, point_record)
            if verbose:
                print(f"[STUDY] Recorded point {point_record['point_id']} into {db_path}")

    return {
        "status": "completed",
        "output_uuid": output_uuid,
        "fit_prefix": fit_prefix,
        "h5_path": written_h5,
        "nc_path": written_nc,
        "table": fit_setup["fit_comps"],
        "trace": trace,
        "data": fit_setup["data"],
        "metadata": fit_setup["fit_params"],
    }


def main():
    args = parse_args()

    # Load configuration
    cfg = HydraConfigManager.get_config(
        environment=args.environment,
        analysis=args.analysis,
        config_dir=args.config_dir,
    )

    co_cut = None if args.co_cut in (None, "None", "null") else args.co_cut
    ndbd_cutoff = None if args.ndbd_cutoff in (None, "None", "null") else args.ndbd_cutoff
    skip_roi = tuple(args.skip_roi) if args.skip_roi else None

    result = run_single_storage_fit(
        cfg=cfg,
        model_key=args.model,
        co_cut=co_cut,
        ndbd_cutoff=ndbd_cutoff,
        co_shift=args.co_shift,
        ndbd_shift=args.ndbd_shift,
        blinded=args.blinded,
        skip_roi=skip_roi,
        prior_width=args.prior_width,
        mkchain_samples=args.mkchain_samples,
        chains=args.chains,
        draws=args.draws,
        tune=args.tune,
        single_core=args.single_core,
        save_trace=args.save_trace,
        study_dir=args.study_dir,
        force=args.force,
        dry_run=args.dry_run,
        verbose=True,
    )

    if result.get("trace") is not None:
        trace = result["trace"]
        table = result["table"]
        print("\n--- Posterior Means ---")
        for var in table.index:
            if var in trace.posterior:
                pmean = float(trace.posterior[var].values.mean())
                pstd = float(trace.posterior[var].values.std())
                print(f"  {var:25s}  mean = {pmean:8.2f} ± {pstd:6.2f}")
        print("-----------------------\n")


if __name__ == "__main__":
    main()
