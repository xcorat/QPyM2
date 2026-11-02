"""
Core Study orchestrator for QPyM2 (v0.6).

Manages study state, checks required input artifacts across all evaluation points,
coordinates on-demand input generation, and exports study manifests.
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import numpy as np
import pandas as pd
from omegaconf import DictConfig

from qpym2.studies.point import StudyPoint
from qpym2.storage.backend import read_fit
from qpym2.storage.study import init_study_db, record_study_point, write_study_manifest


class Study:
    """
    Orchestrator for a collection of parameter evaluation points.

    Attributes:
        name: Name of the study.
        cfg: Merged Hydra configuration.
        points: List of StudyPoint instances composing the study.
        study_dir: Path to study directory.
        description: Human-readable description of study purpose.
    """

    def __init__(
        self,
        name: str,
        cfg: DictConfig,
        points: List[StudyPoint],
        study_dir: Path,
        description: str = "",
    ):
        self.name = name
        self.cfg = cfg
        self.points = points
        self.study_dir = Path(study_dir)
        self.description = description

    def __repr__(self) -> str:
        return f"<Study '{self.name}': {len(self.points)} points at '{self.study_dir}'>"

    def __len__(self) -> int:
        return len(self.points)

    def check_inputs(self) -> pd.DataFrame:
        """
        Check existence of all required inputs on disk across all study points.

        Returns:
            DataFrame with per-point readiness metrics.
        """
        records = [p.check_inputs(self.cfg) for p in self.points]
        return pd.DataFrame(records)

    def get_missing_inputs(self) -> Dict[str, List[str]]:
        """
        Return a mapping from point_id to list of missing input file paths.

        Returns:
            Dictionary of {point_id: [missing_path_1, missing_path_2, ...]}.
        """
        missing_by_point = {}
        for p in self.points:
            report = p.check_inputs(self.cfg)
            missing = list(report.get("missing_template_paths", []))
            if not report.get("data_cached", False):
                missing.append(report.get("data_path", ""))
            if missing:
                missing_by_point[p.point_id] = [m for m in missing if m]

        return missing_by_point

    def print_input_status(self, df: Optional[pd.DataFrame] = None) -> None:
        """
        Print a human-readable summary of input readiness to stdout.

        Args:
            df: Optional pre-computed check_inputs DataFrame.
        """
        status_df = df if df is not None else self.check_inputs()
        total_points = len(status_df)
        ready_points = sum(status_df["is_ready"]) if "is_ready" in status_df else 0
        missing_points = total_points - ready_points

        print("\n" + "=" * 70)
        print(f"Study: {self.name} ({total_points} total points)")
        print(f"Directory: {self.study_dir}")
        if self.description:
            print(f"Description: {self.description}")
        print("=" * 70)

        # Display concise columns
        display_cols = [
            "task_id",
            "point_id",
            "co_shift",
            "co_cut",
            "cached_templates",
            "missing_templates",
            "data_cached",
            "status",
        ]
        available_cols = [c for c in display_cols if c in status_df.columns]
        print(status_df[available_cols].to_string(index=False))

        print("-" * 70)
        if missing_points == 0:
            print(f"✓ All {total_points} points are READY (all inputs cached).")
        else:
            print(f"✗ {missing_points}/{total_points} points have MISSING inputs.")
            missing_dict = self.get_missing_inputs()
            for pid, paths in missing_dict.items():
                print(f"  [{pid}] {len(paths)} missing file(s):")
                for p in paths[:5]:
                    print(f"    - {Path(p).name}")
                if len(paths) > 5:
                    print(f"    ... and {len(paths) - 5} more.")
        print("=" * 70 + "\n")

    def generate_inputs(
        self,
        overwrite: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Materialize all missing (or all if overwrite=True) inputs across the study.

        Uses a shared in-memory histogram cache across points so invariant background
        sources are read from staging only once.

        Args:
            overwrite: If True, re-read and overwrite existing templates.
            verbose: If True, prints progress for each point.

        Returns:
            Dictionary mapping point_id to generation results.
        """
        hist_cache: Dict[Any, Any] = {}
        results: Dict[str, Any] = {}

        for p in self.points:
            report = p.check_inputs(self.cfg)
            if not overwrite and report.get("is_ready", False):
                if verbose:
                    print(f"[{p.point_id}] Inputs already exist. Skipping.")
                results[p.point_id] = {"status": "SKIPPED", "cached": True}
                continue

            if verbose:
                print(f"[{p.point_id}] Generating missing inputs...")
            res = p.generate_inputs(self.cfg, overwrite=overwrite, hist_cache=hist_cache)
            results[p.point_id] = {"status": "GENERATED", "manifest": res}

        return results

    def to_manifest(self) -> Dict[str, Any]:
        """
        Serialize study specification, parameters, and expected artifacts to a dict.

        Returns:
            JSON-serializable manifest dictionary.
        """
        points_info = []
        for p in self.points:
            req = p.get_required_inputs(self.cfg)
            points_info.append(
                {
                    "point_id": p.point_id,
                    "task_id": p.task_id,
                    "model_key": p.model_key,
                    "co_cut": p.co_cut,
                    "ndbd_cutoff": p.ndbd_cutoff,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                    "blinded": p.blinded,
                    "params": p.params,
                    "data": req.get("data", {}),
                    "templates": req.get("templates", {}),
                }
            )

        manifest_data = {
            "study_name": self.name,
            "description": self.description,
            "created_at": datetime.utcnow().isoformat(),
            "study_dir": str(self.study_dir),
            "n_points": len(self.points),
            "points": points_info,
        }
        return manifest_data

    def save_manifest(self, path: Optional[Union[str, Path]] = None) -> Path:
        """
        Save the study manifest JSON to disk.

        Args:
            path: Optional destination path. Defaults to <study_dir>/<name>.manifest.

        Returns:
            Path to written manifest file.
        """
        manifest_path = Path(path) if path else (self.study_dir / f"{self.name}.manifest")
        manifest_data = self.to_manifest()
        return write_study_manifest(manifest_path, manifest_data)

    def check_fits(
        self,
        prior_width: str = "5sig",
        mkchain_samples: int = 100,
        chains: int = 8,
        draws: int = 20000,
        tune: int = 1000,
    ) -> pd.DataFrame:
        """
        Check fit status across all points in the study.

        Returns:
            DataFrame with fit completion metrics for all study points.
        """
        records = [
            p.check_fit(
                self.cfg,
                prior_width=prior_width,
                mkchain_samples=mkchain_samples,
                chains=chains,
                draws=draws,
                tune=tune,
            )
            for p in self.points
        ]
        return pd.DataFrame(records)

    def print_fit_status(self, df: Optional[pd.DataFrame] = None) -> None:
        """
        Print human-readable summary of fit completion status across study points.
        """
        status_df = df if df is not None else self.check_fits()
        total_points = len(status_df)
        fitted_points = sum(status_df["is_fitted"]) if "is_fitted" in status_df else 0

        print("\n" + "=" * 70)
        print(f"Study Fits: {self.name} ({total_points} total points)")
        print(f"Directory: {self.study_dir}")
        print("=" * 70)

        display_cols = ["task_id", "point_id", "co_shift", "output_uuid", "is_fitted", "has_trace", "status"]
        avail_cols = [c for c in display_cols if c in status_df.columns]
        print(status_df[avail_cols].to_string(index=False))

        print("-" * 70)
        if fitted_points == total_points:
            print(f"✓ All {total_points} fits are COMPLETED.")
        else:
            print(f"⏳ {fitted_points}/{total_points} fits completed ({total_points - fitted_points} pending).")
        print("=" * 70 + "\n")

    def run_fits(
        self,
        save_trace: bool = False,
        force: bool = False,
        chains: int = 8,
        draws: int = 20000,
        tune: int = 1000,
        prior_width: str = "5sig",
        mkchain_samples: int = 100,
        single_core: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Execute Bayesian fits across all points in the study.

        Args:
            save_trace: Persist .nc trace if True.
            force: Re-run fits even if already cached.
            chains: Number of MCMC chains.
            draws: Draws per chain.
            tune: Tuning steps.
            prior_width: Prior sigma configuration ('5sig', '10sig', 'syst', '2syst').
            mkchain_samples: JAGS chain samples for priors.
            single_core: Use 1 core for sampler.
            verbose: Verbose logging.

        Returns:
            Dictionary mapping point_id to execution results.
        """
        # Ensure SQLite DB exists
        sqlite_path = self.study_dir / f"{self.name}.sqlite"
        init_study_db(sqlite_path)

        results: Dict[str, Any] = {}
        for p in self.points:
            fit_status = p.check_fit(
                self.cfg,
                prior_width=prior_width,
                mkchain_samples=mkchain_samples,
                chains=chains,
                draws=draws,
                tune=tune,
            )
            if not force and fit_status.get("is_fitted", False):
                if verbose:
                    print(f"[{p.point_id}] Fit {fit_status['output_uuid']} already exists. Skipping.")
                # Ensure point is recorded in SQLite
                record_study_point(sqlite_path, {
                    "point_id": p.point_id,
                    "task_id": p.task_id,
                    "model_key": p.model_key,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                    "co_cut": str(p.co_cut),
                    "output_uuid": fit_status["output_uuid"],
                    "status": "completed",
                })
                results[p.point_id] = {"status": "SKIPPED", "output_uuid": fit_status["output_uuid"]}
                continue

            if verbose:
                print(f"[{p.point_id}] Executing fit for output UUID {fit_status['output_uuid']}...")

            fit_res = p.run_fit(
                cfg=self.cfg,
                save_trace=save_trace,
                force=force,
                chains=chains,
                draws=draws,
                tune=tune,
                prior_width=prior_width,
                mkchain_samples=mkchain_samples,
                single_core=single_core,
                study_dir=self.study_dir,
                verbose=verbose,
            )
            results[p.point_id] = fit_res

        return results

    def collect_results(
        self,
        output_csv: Optional[Union[str, Path]] = None,
        prior_width: str = "5sig",
        chains: int = 8,
        draws: int = 20000,
        tune: int = 1000,
    ) -> pd.DataFrame:
        """
        Aggregate completed fit metrics across all study points into results.csv.

        Extracts posterior summary statistics (means, standard deviations,
        quantiles, and background yields), saves results.csv, and synchronizes
        the study SQLite index.

        Returns:
            Aggregated results DataFrame.
        """
        sqlite_path = self.study_dir / f"{self.name}.sqlite"
        init_study_db(sqlite_path)

        records = []
        for p in self.points:
            fit_status = p.check_fit(
                self.cfg,
                prior_width=prior_width,
                chains=chains,
                draws=draws,
                tune=tune,
            )
            if not fit_status.get("is_fitted", False):
                continue

            fit_prefix = fit_status["fit_prefix"]
            table, trace, data, metadata = read_fit(fit_prefix)

            row: Dict[str, Any] = {
                "point_id": p.point_id,
                "task_id": p.task_id,
                "model_key": p.model_key,
                "co_shift": p.co_shift,
                "ndbd_shift": p.ndbd_shift,
                "co_cut": str(p.co_cut),
                "ndbd_cutoff": str(p.ndbd_cutoff),
                "output_uuid": fit_status["output_uuid"],
                "has_trace": fit_status["has_trace"],
            }

            # If trace is present, extract rich posterior statistics
            if trace is not None and hasattr(trace, "posterior"):
                for var in trace.posterior.data_vars:
                    vals = trace.posterior[var].values
                    row[f"{var}_mean"] = float(vals.mean())
                    row[f"{var}_std"] = float(vals.std())
                    row[f"{var}_q05"] = float(np.percentile(vals, 5))
                    row[f"{var}_q50"] = float(np.percentile(vals, 50))
                    row[f"{var}_q90"] = float(np.percentile(vals, 90))
                    row[f"{var}_q95"] = float(np.percentile(vals, 95))
            elif table is not None and hasattr(table, "index"):
                # Use components table summary
                for comp_name in table.index:
                    if "integral" in table.columns:
                        row[f"{comp_name}_integral"] = float(table.loc[comp_name, "integral"])

            records.append(row)

            # Record in SQLite
            sqlite_rec = dict(row)
            sqlite_rec["status"] = "completed"
            record_study_point(sqlite_path, sqlite_rec)

        results_df = pd.DataFrame(records)
        target_csv = Path(output_csv) if output_csv else (self.study_dir / "results.csv")
        results_df.to_csv(target_csv, index=False)
        print(f"[INFO] Collected results for {len(results_df)} points to: {target_csv}")
        return results_df

