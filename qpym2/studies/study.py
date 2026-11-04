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

    def check_bias_toys(self, n_toys: Optional[int] = None) -> pd.DataFrame:
        """
        Check bias toy completion status across all study points.
        """
        target_n = n_toys
        if target_n is None and hasattr(self.cfg, "bias") and hasattr(self.cfg.bias, "n_toys"):
            target_n = int(self.cfg.bias.n_toys)
        if target_n is None:
            target_n = 200

        records = [p.check_bias_toys(self.study_dir, n_toys=target_n) for p in self.points]
        return pd.DataFrame(records)

    def print_bias_status(self, df: Optional[pd.DataFrame] = None, n_toys: Optional[int] = None) -> None:
        """
        Print human-readable summary of bias toy completion across study points.
        """
        status_df = df if df is not None else self.check_bias_toys(n_toys=n_toys)
        total_points = len(status_df)
        completed_points = sum(status_df["is_completed"]) if "is_completed" in status_df else 0

        print("\n" + "=" * 70)
        print(f"Bias Study Toys: {self.name} ({total_points} total points)")
        print(f"Directory: {self.study_dir}")
        print("=" * 70)

        display_cols = [
            "task_id",
            "point_id",
            "co_shift",
            "n_toys_completed",
            "n_toys_target",
            "is_completed",
            "status",
        ]
        avail_cols = [c for c in display_cols if c in status_df.columns]
        print(status_df[avail_cols].to_string(index=False))

        print("-" * 70)
        if completed_points == total_points:
            print(f"✓ All {total_points} bias toy tasks are COMPLETED.")
        else:
            print(f"⏳ {completed_points}/{total_points} bias toy tasks completed.")
        print("=" * 70 + "\n")

    def run_bias_toys(
        self,
        n_toys: Optional[int] = None,
        bias_min: Optional[float] = None,
        bias_max: Optional[float] = None,
        seed: Optional[int] = None,
        chains: Optional[int] = None,
        draws: Optional[int] = None,
        tune: Optional[int] = None,
        target_accept: Optional[float] = None,
        fit_fn: Optional[Any] = None,
        force: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Execute bias toy fits across all points in the study.

        Pulls default parameters from cfg.bias if available.
        """
        bcfg = getattr(self.cfg, "bias", None)

        def _resolve(val, key, default):
            if val is not None:
                return val
            if bcfg and hasattr(bcfg, key):
                return getattr(bcfg, key)
            return default

        _n_toys = int(_resolve(n_toys, "n_toys", 200))
        _bias_min = float(_resolve(bias_min, "bias_min", 5.0))
        _bias_max = float(_resolve(bias_max, "bias_max", 200.0))
        _seed = int(_resolve(seed, "seed", 42))
        _chains = int(_resolve(chains, "chains", 4))
        _draws = int(_resolve(draws, "draws", 2000))
        _tune = int(_resolve(tune, "tune", 500))
        _target_accept = float(_resolve(target_accept, "target_accept", 0.9))

        results = {}
        for p in self.points:
            if not force:
                st = p.check_bias_toys(self.study_dir, n_toys=_n_toys)
                if st["is_completed"]:
                    if verbose:
                        print(f"[{p.point_id}] Bias toys completed ({st['n_toys_completed']} rows). Skipping.")
                    results[p.point_id] = {"status": "SKIPPED", "task_csv": st["task_csv"]}
                    continue

            if verbose:
                print(f"[{p.point_id}] Running {_n_toys} bias toys (injections in [{_bias_min}, {_bias_max}])...")

            csv_path = p.run_bias_toys(
                cfg=self.cfg,
                study_dir=self.study_dir,
                n_toys=_n_toys,
                bias_min=_bias_min,
                bias_max=_bias_max,
                seed=_seed,
                chains=_chains,
                draws=_draws,
                tune=_tune,
                target_accept=_target_accept,
                fit_fn=fit_fn,
                force=force,
                verbose=verbose,
            )
            results[p.point_id] = {"status": "COMPLETED", "task_csv": str(csv_path)}

        return results

    def collect_bias_results(
        self,
        output_csv: Optional[Union[str, Path]] = None,
        xmin: Optional[float] = None,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Merge bias task CSVs into results.csv and compute calibration parameters via get_bias_from_fitdf.

        Saves bias_calibration.json and synchronizes the study SQLite index.
        """
        import json
        from qpym2.storage.study import collect_study_results, record_study_point
        from qpym2.io.toys import _unflatten_columns
        from qpym2.fit.bias import get_bias_from_fitdf

        bcfg = getattr(self.cfg, "bias", None)
        _xmin = xmin if xmin is not None else (float(getattr(bcfg, "xmin", 15.0)) if bcfg else 15.0)

        sqlite_path = self.study_dir / f"{self.name}.sqlite"
        merged_df = collect_study_results(self.study_dir, self.name)
        if merged_df.empty:
            if verbose:
                print(f"[WARN] No task CSV results found in {self.study_dir}/tasks.")
            return {"results_df": merged_df, "calibration": {}}

        calibrations: Dict[str, Any] = {}
        for p in self.points:
            point_df = merged_df[merged_df["point_id"] == p.point_id].copy()
            if point_df.empty or len(point_df) < 5:
                continue

            unflat = _unflatten_columns(point_df)
            try:
                popt, pcov = get_bias_from_fitdf(unflat, xmin=_xmin)
                a_val = float(popt[0])
                b_val = float(popt[1])
                a_err = float(np.sqrt(max(pcov[0, 0], 0.0)))
                b_err = float(np.sqrt(max(pcov[1, 1], 0.0)))

                calib_entry = {
                    "point_id": p.point_id,
                    "task_id": p.task_id,
                    "co_shift": p.co_shift,
                    "co_cut": str(p.co_cut),
                    "a": a_val,
                    "b": b_val,
                    "a_err": a_err,
                    "b_err": b_err,
                    "pcov_00": float(pcov[0, 0]),
                    "pcov_01": float(pcov[0, 1]),
                    "pcov_10": float(pcov[1, 0]),
                    "pcov_11": float(pcov[1, 1]),
                    "xmin": float(_xmin),
                    "n_toys": len(point_df),
                }
                calibrations[p.point_id] = calib_entry

                # Sync to SQLite
                record_study_point(
                    sqlite_path,
                    {
                        "point_id": p.point_id,
                        "study_name": self.name,
                        "task_id": p.task_id,
                        "status": "completed",
                        "summary_json": calib_entry,
                    },
                )
                if verbose:
                    print(f"[{p.point_id}] Bias fit: a = {a_val:.4f} +/- {a_err:.4f}, b = {b_val:.4f} +/- {b_err:.4f} (xmin={_xmin})")
            except Exception as e:
                if verbose:
                    print(f"[WARN] Could not fit bias for point {p.point_id}: {e}")

        calib_file = self.study_dir / "bias_calibration.json"
        with open(calib_file, "w", encoding="utf-8") as f:
            json.dump(calibrations, f, indent=2)

        if verbose:
            print(f"[INFO] Saved bias calibration parameters for {len(calibrations)} points: {calib_file}")

        return {"results_df": merged_df, "calibration": calibrations}

    def check_sens_toys(self, n_toys: Optional[int] = None) -> pd.DataFrame:
        """
        Check sensitivity toy completion status across all study points.
        """
        target_n = n_toys
        if target_n is None and hasattr(self.cfg, "sens") and hasattr(self.cfg.sens, "n_toys"):
            target_n = int(self.cfg.sens.n_toys)
        if target_n is None:
            target_n = 200

        records = [p.check_sens_toys(self.study_dir, n_toys=target_n) for p in self.points]
        return pd.DataFrame(records)

    def print_sens_status(self, df: Optional[pd.DataFrame] = None, n_toys: Optional[int] = None) -> None:
        """
        Print human-readable summary of sensitivity toy completion across study points.
        """
        status_df = df if df is not None else self.check_sens_toys(n_toys=n_toys)
        total_points = len(status_df)
        completed_points = sum(status_df["is_completed"]) if "is_completed" in status_df else 0

        print("\n" + "=" * 70)
        print(f"Sensitivity Study Toys: {self.name} ({total_points} total points)")
        print(f"Directory: {self.study_dir}")
        print("=" * 70)

        display_cols = [
            "task_id",
            "point_id",
            "co_shift",
            "n_toys_completed",
            "n_toys_target",
            "is_completed",
            "status",
        ]
        avail_cols = [c for c in display_cols if c in status_df.columns]
        print(status_df[avail_cols].to_string(index=False))

        print("-" * 70)
        if completed_points == total_points:
            print(f"✓ All {total_points} sensitivity toy tasks are COMPLETED.")
        else:
            print(f"⏳ {completed_points}/{total_points} sensitivity toy tasks completed.")
        print("=" * 70 + "\n")

    def run_sens_toys(
        self,
        n_toys: Optional[int] = None,
        seed: Optional[int] = None,
        chains: Optional[int] = None,
        draws: Optional[int] = None,
        tune: Optional[int] = None,
        target_accept: Optional[float] = None,
        fit_fn: Optional[Any] = None,
        force: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Execute sensitivity toy fits across all points in the study.

        Pulls default parameters from cfg.sens if available.
        """
        scfg = getattr(self.cfg, "sens", None)

        def _resolve(val, key, default):
            if val is not None:
                return val
            if scfg and hasattr(scfg, key):
                return getattr(scfg, key)
            return default

        _n_toys = int(_resolve(n_toys, "n_toys", 200))
        _seed = int(_resolve(seed, "seed", 42))
        _chains = int(_resolve(chains, "chains", 4))
        _draws = int(_resolve(draws, "draws", 2000))
        _tune = int(_resolve(tune, "tune", 500))
        _target_accept = float(_resolve(target_accept, "target_accept", 0.9))

        results = {}
        for p in self.points:
            if not force:
                st = p.check_sens_toys(self.study_dir, n_toys=_n_toys)
                if st["is_completed"]:
                    if verbose:
                        print(f"[{p.point_id}] Sensitivity toys completed ({st['n_toys_completed']} rows). Skipping.")
                    results[p.point_id] = {"status": "SKIPPED", "task_csv": st["task_csv"]}
                    continue

            if verbose:
                print(f"[{p.point_id}] Running {_n_toys} sensitivity toys (H0: ndbd=0)...")

            csv_path = p.run_sens_toys(
                cfg=self.cfg,
                study_dir=self.study_dir,
                n_toys=_n_toys,
                seed=_seed,
                chains=_chains,
                draws=_draws,
                tune=_tune,
                target_accept=_target_accept,
                fit_fn=fit_fn,
                force=force,
                verbose=verbose,
            )
            results[p.point_id] = {"status": "COMPLETED", "task_csv": str(csv_path)}

        return results

    def collect_sens_results(
        self,
        output_csv: Optional[Union[str, Path]] = None,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Merge sensitivity task CSVs into results.csv and compute summary quantiles per point.

        Calculates expected sensitivity (median ndbd_p90), 1-sigma [p16, p84],
        and 2-sigma [p025, p975] containment intervals. Optionally applies
        bias correction from bias_calibration, saves sens_summary.json, and
        synchronizes the study SQLite index.
        """
        import json
        from qpym2.storage.study import collect_study_results, record_study_point

        sqlite_path = self.study_dir / f"{self.name}.sqlite"
        merged_df = collect_study_results(self.study_dir, self.name)
        if merged_df.empty:
            if verbose:
                print(f"[WARN] No sensitivity task CSV results found in {self.study_dir}/tasks.")
            return {"results_df": merged_df, "summary": {}}

        # Resolve bias calibration if available
        calib_dict: Dict[str, Any] = {}
        if isinstance(bias_calibration, dict):
            calib_dict = bias_calibration
        elif bias_calibration is not None:
            cpath = Path(bias_calibration)
            if not cpath.is_absolute():
                cpath = self.study_dir / cpath if (self.study_dir / cpath).exists() else Path(cpath)
            if cpath.exists():
                with open(cpath, "r", encoding="utf-8") as f:
                    calib_dict = json.load(f)
        else:
            scfg = getattr(self.cfg, "sens", None)
            cfg_calib = getattr(scfg, "bias_calibration", None) if scfg else None
            if cfg_calib:
                cpath = Path(str(cfg_calib))
                if not cpath.is_absolute():
                    cpath = self.study_dir / cpath if (self.study_dir / cpath).exists() else Path(cpath)
                if cpath.exists():
                    with open(cpath, "r", encoding="utf-8") as f:
                        calib_dict = json.load(f)
            elif (self.study_dir / "bias_calibration.json").exists():
                with open(self.study_dir / "bias_calibration.json", "r", encoding="utf-8") as f:
                    calib_dict = json.load(f)

        summaries: Dict[str, Any] = {}
        for p in self.points:
            point_df = merged_df[merged_df["point_id"] == p.point_id].copy()
            if point_df.empty:
                continue

            p90_vals = point_df["ndbd_p90"].dropna().values if "ndbd_p90" in point_df.columns else np.array([])
            mode_vals = point_df["ndbd_mode"].dropna().values if "ndbd_mode" in point_df.columns else np.array([])

            if len(p90_vals) == 0:
                continue

            summary_entry: Dict[str, Any] = {
                "point_id": p.point_id,
                "task_id": p.task_id,
                "co_shift": p.co_shift,
                "ndbd_shift": p.ndbd_shift,
                "co_cut": str(p.co_cut),
                "ndbd_cutoff": str(p.ndbd_cutoff),
                "n_toys": len(point_df),
                "sens_p90_median": float(np.median(p90_vals)),
                "sens_p90_p16": float(np.percentile(p90_vals, 16)),
                "sens_p90_p84": float(np.percentile(p90_vals, 84)),
                "sens_p90_p025": float(np.percentile(p90_vals, 2.5)),
                "sens_p90_p975": float(np.percentile(p90_vals, 97.5)),
                "sens_mode_median": float(np.median(mode_vals)) if len(mode_vals) > 0 else 0.0,
                "sens_mode_p16": float(np.percentile(mode_vals, 16)) if len(mode_vals) > 0 else 0.0,
                "sens_mode_p84": float(np.percentile(mode_vals, 84)) if len(mode_vals) > 0 else 0.0,
            }

            # Check if matching bias calibration entry exists
            p_calib = calib_dict.get(p.point_id)
            if not p_calib:
                # Try finding by matching co_shift
                for ck, cv in calib_dict.items():
                    if isinstance(cv, dict) and cv.get("co_shift") == p.co_shift:
                        p_calib = cv
                        break

            if p_calib and "a" in p_calib and "b" in p_calib:
                a_val = float(p_calib["a"])
                b_val = float(p_calib["b"])
                # Bias-corrected upper limits: (N - b)/a, clipped at 0
                p90_bc = np.clip((p90_vals - b_val) / a_val, 0.0, None)
                summary_entry["bc_a"] = a_val
                summary_entry["bc_b"] = b_val
                summary_entry["sens_bc_p90_median"] = float(np.median(p90_bc))
                summary_entry["sens_bc_p90_p16"] = float(np.percentile(p90_bc, 16))
                summary_entry["sens_bc_p90_p84"] = float(np.percentile(p90_bc, 84))
                summary_entry["sens_bc_p90_p025"] = float(np.percentile(p90_bc, 2.5))
                summary_entry["sens_bc_p90_p975"] = float(np.percentile(p90_bc, 97.5))

            summaries[p.point_id] = summary_entry

            # Sync to SQLite
            record_study_point(
                sqlite_path,
                {
                    "point_id": p.point_id,
                    "study_name": self.name,
                    "task_id": p.task_id,
                    "status": "completed",
                    "summary_json": summary_entry,
                },
            )
            if verbose:
                bc_str = (
                    f" (bias-corrected median: {summary_entry['sens_bc_p90_median']:.2f})"
                    if "sens_bc_p90_median" in summary_entry
                    else ""
                )
                print(
                    f"[{p.point_id}] Expected 90% CL sensitivity: "
                    f"median = {summary_entry['sens_p90_median']:.2f} "
                    f"[-1σ: {summary_entry['sens_p90_p16']:.2f}, +1σ: {summary_entry['sens_p90_p84']:.2f}]{bc_str}"
                )

        summary_file = self.study_dir / "sens_summary.json"
        with open(summary_file, "w", encoding="utf-8") as f:
            json.dump(summaries, f, indent=2)

        if verbose:
            print(f"[INFO] Saved sensitivity summary metrics for {len(summaries)} points: {summary_file}")

        return {"results_df": merged_df, "summary": summaries}



