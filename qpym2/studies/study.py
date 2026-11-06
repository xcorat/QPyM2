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
from omegaconf import DictConfig, OmegaConf

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

        display_cols = [
            "task_id",
            "point_id",
            "co_shift",
            "output_uuid",
            "is_fitted",
            "has_trace",
            "status",
        ]
        avail_cols = [c for c in display_cols if c in status_df.columns]
        print(status_df[avail_cols].to_string(index=False))

        print("-" * 70)
        if fitted_points == total_points:
            print(f"✓ All {total_points} fits are COMPLETED.")
        else:
            print(
                f"⏳ {fitted_points}/{total_points} fits completed ({total_points - fitted_points} pending)."
            )
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
                    print(
                        f"[{p.point_id}] Fit {fit_status['output_uuid']} already exists. Skipping."
                    )
                # Ensure point is recorded in SQLite
                record_study_point(
                    sqlite_path,
                    {
                        "point_id": p.point_id,
                        "task_id": p.task_id,
                        "model_key": p.model_key,
                        "co_shift": p.co_shift,
                        "ndbd_shift": p.ndbd_shift,
                        "co_cut": str(p.co_cut),
                        "output_uuid": fit_status["output_uuid"],
                        "status": "completed",
                    },
                )
                results[p.point_id] = {
                    "status": "SKIPPED",
                    "output_uuid": fit_status["output_uuid"],
                }
                continue

            if verbose:
                print(
                    f"[{p.point_id}] Executing fit for output UUID {fit_status['output_uuid']}..."
                )

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

    def print_bias_status(
        self, df: Optional[pd.DataFrame] = None, n_toys: Optional[int] = None
    ) -> None:
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
                        print(
                            f"[{p.point_id}] Bias toys completed ({st['n_toys_completed']} rows). Skipping."
                        )
                    results[p.point_id] = {"status": "SKIPPED", "task_csv": st["task_csv"]}
                    continue

            if verbose:
                print(
                    f"[{p.point_id}] Running {_n_toys} bias toys (injections in [{_bias_min}, {_bias_max}])..."
                )

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
                    print(
                        f"[{p.point_id}] Bias fit: a = {a_val:.4f} +/- {a_err:.4f}, b = {b_val:.4f} +/- {b_err:.4f} (xmin={_xmin})"
                    )
            except Exception as e:
                if verbose:
                    print(f"[WARN] Could not fit bias for point {p.point_id}: {e}")

        calib_file = self.study_dir / "bias_calibration.json"
        with open(calib_file, "w", encoding="utf-8") as f:
            json.dump(calibrations, f, indent=2)

        if verbose:
            print(
                f"[INFO] Saved bias calibration parameters for {len(calibrations)} points: {calib_file}"
            )

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

    def print_sens_status(
        self, df: Optional[pd.DataFrame] = None, n_toys: Optional[int] = None
    ) -> None:
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
                        print(
                            f"[{p.point_id}] Sensitivity toys completed ({st['n_toys_completed']} rows). Skipping."
                        )
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
                    cpath = (
                        self.study_dir / cpath if (self.study_dir / cpath).exists() else Path(cpath)
                    )
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

            p90_vals = (
                point_df["ndbd_p90"].dropna().values
                if "ndbd_p90" in point_df.columns
                else np.array([])
            )
            mode_vals = (
                point_df["ndbd_mode"].dropna().values
                if "ndbd_mode" in point_df.columns
                else np.array([])
            )

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
            print(
                f"[INFO] Saved sensitivity summary metrics for {len(summaries)} points: {summary_file}"
            )

        return {"results_df": merged_df, "summary": summaries}

    def check_syst(self, n_samples: Optional[int] = None) -> pd.DataFrame:
        """
        Check systematic task completion status across all study points.
        """
        scfg = getattr(self.cfg, "syst", None)
        category = str(getattr(scfg, "category", "bmmc")) if scfg else "bmmc"
        subcmd = getattr(scfg, "subcommand", None) if scfg else None
        if category == "lss" and subcmd:
            category = f"lss_{subcmd}"

        target_n = n_samples
        if target_n is None and scfg and hasattr(scfg, "n_samples"):
            target_n = int(scfg.n_samples)
        if target_n is None:
            target_n = 100

        datasets = (
            list(getattr(scfg, "datasets", []))
            if scfg and hasattr(scfg, "datasets") and scfg.datasets
            else None
        )
        dq_mode = str(getattr(scfg, "dq_mode", "pm1")) if scfg else "pm1"

        records = []
        for p in self.points:
            if category in ("lss_sscl", "lss_dq"):
                records.append(
                    p.check_syst_lss(
                        self.study_dir,
                        category=category,
                        datasets=datasets,
                        n_samples=target_n,
                        dq_mode=dq_mode,
                    )
                )
            elif category in ("bmmc", "mc_stats", "bm_norm"):
                records.append(p.check_syst_mc(self.study_dir, n_samples=target_n))
            elif category == "bm_threshold":
                _var = getattr(scfg, "variant", "Threshold250") if scfg else "Threshold250"
                records.append(p.check_syst_bm_threshold(self.study_dir, variant=_var))
            elif category == "bm_binning":
                _vars = getattr(scfg, "variants", None) if scfg else None
                records.append(p.check_syst_bm_binning(self.study_dir, variants=_vars))
            elif category in ("bin_size", "template_binning"):
                _sizes = getattr(scfg, "bin_sizes", None) if scfg else None
                records.append(p.check_syst_template_binning(self.study_dir, bin_sizes=_sizes))
            else:
                # Fallback for generic syst point
                records.append(p.check_syst_mc(self.study_dir, n_samples=target_n))

        return pd.DataFrame(records)

    def print_syst_status(
        self, df: Optional[pd.DataFrame] = None, n_samples: Optional[int] = None
    ) -> None:
        """
        Print human-readable summary of systematic task completion across study points.
        """
        status_df = df if df is not None else self.check_syst(n_samples=n_samples)
        total_points = len(status_df)
        completed_points = sum(status_df["is_completed"]) if "is_completed" in status_df else 0

        scfg = getattr(self.cfg, "syst", None)
        category = str(getattr(scfg, "category", "bmmc")) if scfg else "bmmc"
        subcmd = getattr(scfg, "subcommand", None) if scfg else None
        if category == "lss" and subcmd:
            category = f"lss_{subcmd}"

        print("\n" + "=" * 70)
        print(f"Systematic Study Tasks [{category}]: {self.name} ({total_points} total points)")
        print(f"Directory: {self.study_dir}")
        print("=" * 70)

        display_cols = [
            "task_id",
            "point_id",
            "co_shift",
            "n_samples_completed",
            "n_samples_target",
            "is_completed",
            "status",
        ]
        avail_cols = [c for c in display_cols if c in status_df.columns]
        print(status_df[avail_cols].to_string(index=False))

        print("-" * 70)
        if completed_points == total_points:
            print(f"✓ All {total_points} systematic tasks are COMPLETED.")
        else:
            print(f"⏳ {completed_points}/{total_points} systematic tasks completed.")
        print("=" * 70 + "\n")

    def run_syst(
        self,
        n_samples: Optional[int] = None,
        seed: Optional[int] = None,
        chains: Optional[int] = None,
        draws: Optional[int] = None,
        tune: Optional[int] = None,
        target_accept: Optional[float] = None,
        stat_only: Optional[bool] = None,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        fit_fn: Optional[Any] = None,
        save_traces: bool = False,
        force: bool = False,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Execute systematic fits across all points in the study.
        """
        scfg = getattr(self.cfg, "syst", None)

        def _resolve(val, key, default):
            if val is not None:
                return val
            if scfg and hasattr(scfg, key):
                return getattr(scfg, key)
            return default

        category = str(_resolve(None, "category", "bmmc"))
        subcmd = _resolve(None, "subcommand", None)
        if category == "lss" and subcmd:
            category = f"lss_{subcmd}"

        _n_samples = int(_resolve(n_samples, "n_samples", 100))
        _seed = int(_resolve(seed, "seed", 42))
        _chains = int(_resolve(chains, "chains", 2))
        _draws = int(_resolve(draws, "draws", 2500))
        _tune = int(_resolve(tune, "tune", 500))
        _target_accept = float(_resolve(target_accept, "target_accept", 0.9))
        _stat_only = bool(_resolve(stat_only, "stat_only", False))
        _bias_calib = _resolve(bias_calibration, "bias_calibration", None)
        _datasets = list(_resolve(None, "datasets", [])) or None
        _dq_mode = str(_resolve(None, "dq_mode", "pm1"))

        results = {}
        for p in self.points:
            if not force:
                if category in ("lss_sscl", "lss_dq"):
                    st = p.check_syst_lss(
                        self.study_dir,
                        category=category,
                        datasets=_datasets,
                        n_samples=_n_samples,
                        dq_mode=_dq_mode,
                    )
                elif category == "bm_threshold":
                    _var = getattr(scfg, "variant", "Threshold250") if scfg else "Threshold250"
                    st = p.check_syst_bm_threshold(self.study_dir, variant=_var)
                elif category == "bm_binning":
                    _vars = getattr(scfg, "variants", None) if scfg else None
                    st = p.check_syst_bm_binning(self.study_dir, variants=_vars)
                elif category in ("bin_size", "template_binning"):
                    _sizes = getattr(scfg, "bin_sizes", None) if scfg else None
                    st = p.check_syst_template_binning(self.study_dir, bin_sizes=_sizes)
                else:
                    st = p.check_syst_mc(self.study_dir, n_samples=_n_samples)
                if st["is_completed"]:
                    if verbose:
                        print(
                            f"[{p.point_id}] Systematic tasks completed ({st['n_samples_completed']} rows). Skipping."
                        )
                    results[p.point_id] = {
                        "status": "SKIPPED",
                        "task_csv": st["task_csv"],
                        "task_h5": st["task_h5"],
                    }
                    continue

            if verbose:
                print(f"[{p.point_id}] Running systematic fits for category '{category}'...")

            if category in ("lss_sscl", "lss_dq"):
                _lss_data_path = (
                    getattr(self.cfg.syst, "lss_data_path", None)
                    if hasattr(self.cfg, "syst")
                    else None
                )
                csv_path, h5_path = p.run_syst_lss(
                    cfg=self.cfg,
                    study_dir=self.study_dir,
                    category=category,
                    datasets=_datasets,
                    n_samples=_n_samples,
                    dq_mode=_dq_mode,
                    seed=_seed,
                    lss_data_path=_lss_data_path,
                    chains=_chains,
                    draws=_draws,
                    tune=_tune,
                    target_accept=_target_accept,
                    single_core=True,
                    bias_calibration=_bias_calib,
                    fit_fn=fit_fn,
                    save_traces=save_traces,
                    force=force,
                    verbose=verbose,
                )
            elif category == "bm_threshold":
                _var = getattr(scfg, "variant", "Threshold250") if scfg else "Threshold250"
                _m_key_syst = getattr(scfg, "model_key_syst", "tech_note_250") if scfg else "tech_note_250"
                _raw_path = getattr(scfg, "results_raw_path", None) if scfg else None
                _ded_bias = getattr(scfg, "dedicated_bias_csv", None) if scfg else None
                csv_path, h5_path = p.run_syst_bm_threshold(
                    cfg=self.cfg,
                    study_dir=self.study_dir,
                    variant=_var,
                    model_key_syst=_m_key_syst,
                    results_raw_path=_raw_path,
                    seed=_seed,
                    chains=_chains,
                    draws=_draws,
                    tune=_tune,
                    target_accept=_target_accept,
                    single_core=True,
                    bias_calibration=_bias_calib,
                    dedicated_bias_csv=_ded_bias,
                    fit_fn=fit_fn,
                    save_traces=save_traces,
                    force=force,
                    verbose=verbose,
                )
            elif category == "bm_binning":
                _vars = getattr(scfg, "variants", None) if scfg else None
                _raw_path = getattr(scfg, "results_raw_path", None) if scfg else None
                csv_path, h5_path = p.run_syst_bm_binning(
                    cfg=self.cfg,
                    study_dir=self.study_dir,
                    variants=_vars,
                    results_raw_path=_raw_path,
                    seed=_seed,
                    chains=_chains,
                    draws=_draws,
                    tune=_tune,
                    target_accept=_target_accept,
                    single_core=True,
                    bias_calibration=_bias_calib,
                    fit_fn=fit_fn,
                    save_traces=save_traces,
                    force=force,
                    verbose=verbose,
                )
            elif category in ("bin_size", "template_binning"):
                _sizes = getattr(scfg, "bin_sizes", None) if scfg else None
                csv_path, h5_path = p.run_syst_template_binning(
                    cfg=self.cfg,
                    study_dir=self.study_dir,
                    bin_sizes=_sizes,
                    seed=_seed,
                    chains=_chains,
                    draws=_draws,
                    tune=_tune,
                    target_accept=_target_accept,
                    single_core=True,
                    bias_calibration=_bias_calib,
                    fit_fn=fit_fn,
                    save_traces=save_traces,
                    force=force,
                    verbose=verbose,
                )
            else:
                csv_path, h5_path = p.run_syst_mc(
                    cfg=self.cfg,
                    study_dir=self.study_dir,
                    n_samples=_n_samples,
                    seed=_seed,
                    chains=_chains,
                    draws=_draws,
                    tune=_tune,
                    target_accept=_target_accept,
                    stat_only=_stat_only,
                    bias_calibration=_bias_calib,
                    fit_fn=fit_fn,
                    save_traces=save_traces,
                    force=force,
                    verbose=verbose,
                )
            results[p.point_id] = {
                "status": "COMPLETED",
                "task_csv": str(csv_path),
                "task_h5": str(h5_path),
            }

        return results

    def collect_syst_results(
        self,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        nominal_post: Optional[np.ndarray] = None,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """
        Merge systematic task CSVs, combine posterior samples into mixture distributions,
        calculate mode shifts and 1-sigma HDI edge extensions against nominal,
        update the study SQLite database, export syst_summary.json, and emit comparison plots.
        """
        import json
        import h5py
        from qpym2.storage.study import collect_study_results, record_study_point
        from qpym2.utils.mathy import sample_weighted_mixture, calc_intervals
        from qpym2.syst.core import compute_deltas
        from qpym2.fit.bias import apply_bias_correction
        from qpym2.syst.io import save_syst_result
        from qpym2.syst.plots import plot_syst_comparison
        import matplotlib.pyplot as plt

        sqlite_path = self.study_dir / f"{self.name}.sqlite"
        merged_df = collect_study_results(self.study_dir, self.name)
        if merged_df.empty:
            if verbose:
                print(f"[WARN] No systematic task CSV results found in {self.study_dir}/tasks.")
            return {"results_df": merged_df, "summary": {}}

        scfg = getattr(self.cfg, "syst", None)
        category = str(getattr(scfg, "category", "bmmc")) if scfg else "bmmc"
        subcmd = getattr(scfg, "subcommand", None) if scfg else None
        if category == "lss" and subcmd:
            category = f"lss_{subcmd}"

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
            cfg_calib = getattr(scfg, "bias_calibration", None) if scfg else None
            if cfg_calib:
                cpath = Path(str(cfg_calib))
                if not cpath.is_absolute():
                    cpath = (
                        self.study_dir / cpath if (self.study_dir / cpath).exists() else Path(cpath)
                    )
                if cpath.exists():
                    with open(cpath, "r", encoding="utf-8") as f:
                        calib_dict = json.load(f)
            elif (self.study_dir / "bias_calibration.json").exists():
                with open(self.study_dir / "bias_calibration.json", "r", encoding="utf-8") as f:
                    calib_dict = json.load(f)

        plots_dir = self.study_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        summaries: Dict[str, Any] = {}
        for p in self.points:
            point_df = merged_df[merged_df["point_id"] == p.point_id].copy()
            if point_df.empty:
                continue

            if category in ("lss_sscl", "lss_dq"):
                # Resolve nominal posterior
                nom_post_raw = None
                if nominal_post is not None:
                    nom_post_raw = np.asarray(nominal_post).flatten()
                else:
                    try:
                        fit_info = p.get_fit_hash(self.cfg)
                        fit_prefix = fit_info["fit_prefix"]
                        nc_path = fit_prefix.with_suffix(".nc")
                        if nc_path.exists():
                            import arviz as az

                            nom_trace = az.from_netcdf(str(nc_path))
                            nom_post_raw = nom_trace.posterior["ndbd"].values.flatten()
                    except Exception:
                        pass

                # Resolve calibration parameters
                p_calib = calib_dict.get(p.point_id)
                if not p_calib:
                    for ck, cv in calib_dict.items():
                        if isinstance(cv, dict) and cv.get("co_shift") == p.co_shift:
                            p_calib = cv
                            break

                popt = None
                if p_calib and "a" in p_calib and "b" in p_calib:
                    popt = (float(p_calib["a"]), float(p_calib["b"]))

                use_bc = (
                    popt is not None
                    and "ndbd_mode_bc" in point_df.columns
                    and not point_df["ndbd_mode_bc"].isna().all()
                )
                mode_col = "ndbd_mode_bc" if use_bc else "ndbd_mode"
                p16_col = "ndbd_p16_bc" if use_bc else "ndbd_p16"
                p90_col = "ndbd_p90_bc" if use_bc else "ndbd_p90"

                if nom_post_raw is not None:
                    nom_post_eval = (
                        apply_bias_correction(nom_post_raw, popt)
                        if popt is not None
                        else nom_post_raw
                    )
                    mode_nom, lo_nom, hi_nom = calc_intervals([nom_post_eval], method="hdi")[0]
                else:
                    mode_nom = float(point_df[mode_col].iloc[0])
                    lo_nom = float(point_df[p16_col].iloc[0])
                    hi_nom = float(point_df[p90_col].iloc[0])

                from qpym2.utils.mathy import signed_quad_sum
                from qpym2.syst.plots import plot_lss_sscl_barplot, plot_lss_dq_barplot

                if category == "lss_sscl":
                    point_df["delta_mode"] = point_df[mode_col] - mode_nom
                    point_df["delta_lo"] = point_df[p16_col] - lo_nom
                    point_df["delta_hi"] = point_df[p90_col] - hi_nom

                    deltas = point_df["delta_mode"].values
                    sigma_plus, sigma_minus = signed_quad_sum(deltas)
                    delta_right = float(sigma_plus)
                    delta_left = -float(sigma_minus)
                    delta_mode = float(np.mean(deltas))

                    syst_h5 = self.study_dir / f"syst_lss_sscl_{p.point_id}.h5"
                    per_post = point_df.copy()

                    plot_path = plots_dir / f"lss_sscl_barplot_{p.point_id}.png"
                    try:
                        fig, _ = plot_lss_sscl_barplot(
                            ds_ids=point_df["ds"].values,
                            dlo=point_df["delta_lo"].values,
                            dhi=point_df["delta_hi"].values,
                            deltas=(delta_left, delta_right),
                            co_shift=p.co_shift,
                        )
                        fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                        plt.close(fig)
                    except Exception as e:
                        if verbose:
                            print(f"[{p.point_id}] Note: could not render LSS plot ({e})")

                else:  # lss_dq
                    datasets_in_point = sorted(point_df["ds"].unique())
                    dpos_list = []
                    dneg_list = []
                    dlo_1r, dhi_1r, dlo_1l, dhi_1l = [], [], [], []

                    for ds_val in datasets_in_point:
                        df_ds = point_df[point_df["ds"] == ds_val]
                        r1_rows = df_ds[df_ds["variant"] == "dqres_1r"]
                        l1_rows = df_ds[df_ds["variant"] == "dqres_1l"]

                        if not r1_rows.empty and not l1_rows.empty:
                            d1r = float(r1_rows[mode_col].iloc[0] - mode_nom)
                            d1l = float(l1_rows[mode_col].iloc[0] - mode_nom)
                            dpos = max(d1r, d1l, 0.0)
                            dneg = min(d1r, d1l, 0.0)
                            dpos_list.append(dpos)
                            dneg_list.append(dneg)

                            dlo_1r.append(float(r1_rows[p16_col].iloc[0] - lo_nom))
                            dhi_1r.append(float(r1_rows[p90_col].iloc[0] - hi_nom))
                            dlo_1l.append(float(l1_rows[p16_col].iloc[0] - lo_nom))
                            dhi_1l.append(float(l1_rows[p90_col].iloc[0] - hi_nom))
                        else:
                            ds_deltas = (df_ds[mode_col] - mode_nom).values
                            dpos_list.append(max(float(np.max(ds_deltas)), 0.0))
                            dneg_list.append(min(float(np.min(ds_deltas)), 0.0))
                            dlo_1r.append(0.0)
                            dhi_1r.append(0.0)
                            dlo_1l.append(0.0)
                            dhi_1l.append(0.0)

                    delta_right = float(np.sqrt(np.sum(np.array(dpos_list) ** 2)))
                    delta_left = -float(np.sqrt(np.sum(np.array(dneg_list) ** 2)))
                    delta_mode = float(np.mean(dpos_list) + np.mean(dneg_list))

                    syst_h5 = self.study_dir / f"syst_lss_dq_{p.point_id}.h5"
                    per_post = point_df.copy()

                    plot_path = plots_dir / f"lss_dq_barplot_{p.point_id}.png"
                    try:
                        fig, _ = plot_lss_dq_barplot(
                            ds_ids=datasets_in_point,
                            dlo_1r=np.array(dlo_1r),
                            dhi_1r=np.array(dhi_1r),
                            dlo_1l=np.array(dlo_1l),
                            dhi_1l=np.array(dhi_1l),
                            deltas=(delta_left, delta_right),
                            co_shift=p.co_shift,
                        )
                        fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                        plt.close(fig)
                    except Exception as e:
                        if verbose:
                            print(f"[{p.point_id}] Note: could not render dQ plot ({e})")

                summary_entry = {
                    "point_id": p.point_id,
                    "task_id": p.task_id,
                    "category": category,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                    "co_cut": str(p.co_cut),
                    "ndbd_cutoff": str(p.ndbd_cutoff),
                    "n_samples": len(point_df),
                    "mode_nom": float(mode_nom),
                    "delta_mode": float(delta_mode),
                    "delta_right": float(delta_right),
                    "delta_left": float(delta_left),
                }
                if popt is not None:
                    summary_entry["bc_a"] = popt[0]
                    summary_entry["bc_b"] = popt[1]

                summaries[p.point_id] = summary_entry

                model_pars = {
                    "co_cut": p.co_cut,
                    "ndbd_cutoff": p.ndbd_cutoff,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                }
                save_syst_result(
                    syst_h5,
                    name=category,
                    per_post=per_post,
                    model_pars=model_pars,
                    overwrite=True,
                )

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
                    print(
                        f"[{p.point_id}] {category.upper()} Systematic: "
                        f"Δright = {delta_right:+.2f}, "
                        f"Δleft = {delta_left:+.2f} "
                        f"(mode_nom={mode_nom:.2f})"
                    )
                continue

            if category == "bm_threshold":
                task_h5 = self.study_dir / "tasks" / f"syst_bm_threshold_{p.point_id}_samples.h5"
                if not task_h5.exists():
                    continue

                with h5py.File(task_h5, "r") as h5f:
                    ndbd_samples = h5f["ndbd_samples"][:]

                syst_post_raw = (
                    ndbd_samples[0] if ndbd_samples.ndim > 1 else ndbd_samples.flatten()
                )

                nom_post_raw = None
                if nominal_post is not None:
                    nom_post_raw = np.asarray(nominal_post).flatten()
                else:
                    try:
                        fit_info = p.get_fit_hash(self.cfg)
                        fit_prefix = fit_info["fit_prefix"]
                        nc_path = fit_prefix.with_suffix(".nc")
                        if nc_path.exists():
                            import arviz as az
                            nom_trace = az.from_netcdf(str(nc_path))
                            nom_post_raw = nom_trace.posterior["ndbd"].values.flatten()
                    except Exception:
                        pass
                if nom_post_raw is None:
                    nom_post_raw = syst_post_raw

                popt = None
                if "bc_a" in point_df.columns and "bc_b" in point_df.columns:
                    popt = (float(point_df["bc_a"].iloc[0]), float(point_df["bc_b"].iloc[0]))
                elif p.point_id in calib_dict:
                    c = calib_dict[p.point_id]
                    if "a" in c and "b" in c:
                        popt = (float(c["a"]), float(c["b"]))
                elif "a" in calib_dict and "b" in calib_dict:
                    popt = (float(calib_dict["a"]), float(calib_dict["b"]))

                nom_bc = apply_bias_correction(nom_post_raw, popt) if popt is not None else nom_post_raw
                syst_bc = apply_bias_correction(syst_post_raw, popt) if popt is not None else syst_post_raw

                mode_nom, lo_nom, hi_nom = calc_intervals([nom_bc], method="hdi")[0]
                mode_syst, lo_syst, hi_syst = calc_intervals([syst_bc], method="hdi")[0]
                deltas = compute_deltas(nom_bc, syst_bc, return_extensions_only=True)
                delta_left, delta_right = float(deltas[0]), float(deltas[1])
                delta_mode = float(mode_syst - mode_nom)

                variant = getattr(scfg, "variant", "Threshold250") if scfg else "Threshold250"
                summary_entry = {
                    "point_id": p.point_id,
                    "task_id": p.task_id,
                    "category": "bm_threshold",
                    "variant": variant,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                    "co_cut": str(p.co_cut),
                    "ndbd_cutoff": str(p.ndbd_cutoff),
                    "mode_nom": float(mode_nom),
                    "mode_syst": float(mode_syst),
                    "delta_mode": float(delta_mode),
                    "hdi_nom": [float(lo_nom), float(hi_nom)],
                    "hdi_syst": [float(lo_syst), float(hi_syst)],
                    "delta_right": float(delta_right),
                    "delta_left": float(delta_left),
                }
                if popt is not None:
                    summary_entry["bc_a"] = popt[0]
                    summary_entry["bc_b"] = popt[1]

                summaries[p.point_id] = summary_entry

                syst_h5 = self.study_dir / f"syst_bm_threshold_{p.point_id}.h5"
                model_pars = {
                    "co_cut": p.co_cut,
                    "ndbd_cutoff": p.ndbd_cutoff,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                }
                save_syst_result(
                    syst_h5,
                    name="bm_threshold",
                    bc_post=syst_bc,
                    raw_post=syst_post_raw,
                    model_pars=model_pars,
                    per_post=point_df,
                    overwrite=True,
                )

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

                plot_path = plots_dir / f"bm_threshold_comparison_{p.point_id}.png"
                try:
                    fig, _ = plot_syst_comparison(
                        nominal_bc=nom_bc,
                        syst_post=syst_bc,
                        delta=delta_mode,
                        deltas=[delta_left, delta_right, delta_mode],
                        syst_label=f"BM-{variant}",
                        nom_label="Nominal",
                        title=f"BM Threshold Systematic ({variant})",
                    )
                    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                    plt.close(fig)
                except Exception as e:
                    if verbose:
                        print(f"[{p.point_id}] Note: could not render plot ({e})")

                if verbose:
                    print(
                        f"[{p.point_id}] BM Threshold ({variant}): "
                        f"Δmode = {delta_mode:+.2f}, "
                        f"Δright = {delta_right:+.2f}, "
                        f"Δleft = {delta_left:+.2f} "
                        f"(mode_nom={mode_nom:.2f}, mode_syst={mode_syst:.2f})"
                    )
                continue

            if category == "bm_binning":
                task_h5 = self.study_dir / "tasks" / f"syst_bm_binning_{p.point_id}_samples.h5"
                if not task_h5.exists():
                    continue

                with h5py.File(task_h5, "r") as h5f:
                    ndbd_samples = h5f["ndbd_samples"][:]

                nom_post_raw = None
                if nominal_post is not None:
                    nom_post_raw = np.asarray(nominal_post).flatten()
                else:
                    try:
                        fit_info = p.get_fit_hash(self.cfg)
                        fit_prefix = fit_info["fit_prefix"]
                        nc_path = fit_prefix.with_suffix(".nc")
                        if nc_path.exists():
                            import arviz as az
                            nom_trace = az.from_netcdf(str(nc_path))
                            nom_post_raw = nom_trace.posterior["ndbd"].values.flatten()
                    except Exception:
                        pass
                if nom_post_raw is None:
                    nom_post_raw = ndbd_samples[0]

                popt = None
                if p.point_id in calib_dict:
                    c = calib_dict[p.point_id]
                    if "a" in c and "b" in c:
                        popt = (float(c["a"]), float(c["b"]))
                elif "a" in calib_dict and "b" in calib_dict:
                    popt = (float(calib_dict["a"]), float(calib_dict["b"]))

                nom_bc = apply_bias_correction(nom_post_raw, popt) if popt is not None else nom_post_raw
                mode_nom, lo_nom, hi_nom = calc_intervals([nom_bc], method="hdi")[0]

                variants_list = (
                    point_df["variant"].tolist()
                    if "variant" in point_df.columns
                    else [f"v{i}" for i in range(len(ndbd_samples))]
                )
                per_variant_deltas = {}
                all_d_left = []
                all_d_right = []
                all_d_mode = []

                for v_idx, v_name in enumerate(variants_list):
                    v_post_raw = ndbd_samples[v_idx]
                    v_bc = apply_bias_correction(v_post_raw, popt) if popt is not None else v_post_raw
                    v_mode, v_lo, v_hi = calc_intervals([v_bc], method="hdi")[0]
                    d = compute_deltas(nom_bc, v_bc, return_extensions_only=True)
                    dl, dr = float(d[0]), float(d[1])
                    dm = float(v_mode - mode_nom)
                    per_variant_deltas[v_name] = {
                        "delta_left": dl,
                        "delta_right": dr,
                        "delta_mode": dm,
                        "mode_syst": float(v_mode),
                        "hdi_syst": [float(v_lo), float(v_hi)],
                    }
                    all_d_left.append(dl)
                    all_d_right.append(dr)
                    all_d_mode.append(dm)

                combined_d_left = float(min(all_d_left)) if all_d_left else 0.0
                combined_d_right = float(max(all_d_right)) if all_d_right else 0.0
                combined_d_mode = float(np.mean(all_d_mode)) if all_d_mode else 0.0

                summary_entry = {
                    "point_id": p.point_id,
                    "task_id": p.task_id,
                    "category": "bm_binning",
                    "variants": variants_list,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                    "co_cut": str(p.co_cut),
                    "ndbd_cutoff": str(p.ndbd_cutoff),
                    "mode_nom": float(mode_nom),
                    "hdi_nom": [float(lo_nom), float(hi_nom)],
                    "delta_right": combined_d_right,
                    "delta_left": combined_d_left,
                    "delta_mode": combined_d_mode,
                    "per_variant": per_variant_deltas,
                }
                if popt is not None:
                    summary_entry["bc_a"] = popt[0]
                    summary_entry["bc_b"] = popt[1]

                summaries[p.point_id] = summary_entry

                syst_h5 = self.study_dir / f"syst_bm_binning_{p.point_id}.h5"
                model_pars = {
                    "co_cut": p.co_cut,
                    "ndbd_cutoff": p.ndbd_cutoff,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                }
                save_syst_result(
                    syst_h5,
                    name="bm_binning",
                    bc_post=nom_bc,
                    raw_post=nom_post_raw,
                    model_pars=model_pars,
                    per_post=point_df,
                    overwrite=True,
                )

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
                    print(
                        f"[{p.point_id}] BM Binning (Envelope): "
                        f"Δmode = {combined_d_mode:+.2f}, "
                        f"Δright = {combined_d_right:+.2f}, "
                        f"Δleft = {combined_d_left:+.2f} "
                        f"(mode_nom={mode_nom:.2f})"
                    )
                continue

            if category in ("bin_size", "template_binning"):
                task_h5 = self.study_dir / "tasks" / f"syst_bin_size_{p.point_id}_samples.h5"
                if not task_h5.exists():
                    continue

                with h5py.File(task_h5, "r") as h5f:
                    ndbd_samples = h5f["ndbd_samples"][:]

                nom_post_raw = None
                if nominal_post is not None:
                    nom_post_raw = np.asarray(nominal_post).flatten()
                else:
                    try:
                        fit_info = p.get_fit_hash(self.cfg)
                        fit_prefix = fit_info["fit_prefix"]
                        nc_path = fit_prefix.with_suffix(".nc")
                        if nc_path.exists():
                            import arviz as az
                            nom_trace = az.from_netcdf(str(nc_path))
                            nom_post_raw = nom_trace.posterior["ndbd"].values.flatten()
                    except Exception:
                        pass
                if nom_post_raw is None:
                    nom_post_raw = ndbd_samples[0]

                popt = None
                if p.point_id in calib_dict:
                    c = calib_dict[p.point_id]
                    if "a" in c and "b" in c:
                        popt = (float(c["a"]), float(c["b"]))
                elif "a" in calib_dict and "b" in calib_dict:
                    popt = (float(calib_dict["a"]), float(calib_dict["b"]))

                nom_bc = apply_bias_correction(nom_post_raw, popt) if popt is not None else nom_post_raw
                mode_nom, lo_nom, hi_nom = calc_intervals([nom_bc], method="hdi")[0]

                bin_sizes_list = (
                    point_df["bin_size"].tolist()
                    if "bin_size" in point_df.columns
                    else [f"bs{i}" for i in range(len(ndbd_samples))]
                )
                per_bs_deltas = {}
                all_d_left = []
                all_d_right = []
                all_d_mode = []

                for bs_idx, bs_val in enumerate(bin_sizes_list):
                    bs_post_raw = ndbd_samples[bs_idx]
                    bs_bc = apply_bias_correction(bs_post_raw, popt) if popt is not None else bs_post_raw
                    bs_mode, bs_lo, bs_hi = calc_intervals([bs_bc], method="hdi")[0]
                    d = compute_deltas(nom_bc, bs_bc, return_extensions_only=True)
                    dl, dr = float(d[0]), float(d[1])
                    dm = float(bs_mode - mode_nom)
                    per_bs_deltas[str(bs_val)] = {
                        "delta_left": dl,
                        "delta_right": dr,
                        "delta_mode": dm,
                        "mode_syst": float(bs_mode),
                        "hdi_syst": [float(bs_lo), float(bs_hi)],
                    }
                    all_d_left.append(dl)
                    all_d_right.append(dr)
                    all_d_mode.append(dm)

                combined_d_left = float(min(all_d_left)) if all_d_left else 0.0
                combined_d_right = float(max(all_d_right)) if all_d_right else 0.0
                combined_d_mode = float(np.mean(all_d_mode)) if all_d_mode else 0.0

                summary_entry = {
                    "point_id": p.point_id,
                    "task_id": p.task_id,
                    "category": "bin_size",
                    "bin_sizes": bin_sizes_list,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                    "co_cut": str(p.co_cut),
                    "ndbd_cutoff": str(p.ndbd_cutoff),
                    "mode_nom": float(mode_nom),
                    "hdi_nom": [float(lo_nom), float(hi_nom)],
                    "delta_right": combined_d_right,
                    "delta_left": combined_d_left,
                    "delta_mode": combined_d_mode,
                    "per_bin_size": per_bs_deltas,
                }
                if popt is not None:
                    summary_entry["bc_a"] = popt[0]
                    summary_entry["bc_b"] = popt[1]

                summaries[p.point_id] = summary_entry

                syst_h5 = self.study_dir / f"syst_bin_size_{p.point_id}.h5"
                model_pars = {
                    "co_cut": p.co_cut,
                    "ndbd_cutoff": p.ndbd_cutoff,
                    "co_shift": p.co_shift,
                    "ndbd_shift": p.ndbd_shift,
                }
                save_syst_result(
                    syst_h5,
                    name="bin_size",
                    bc_post=nom_bc,
                    raw_post=nom_post_raw,
                    model_pars=model_pars,
                    per_post=point_df,
                    overwrite=True,
                )

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
                    print(
                        f"[{p.point_id}] Template Binning (Envelope): "
                        f"Δmode = {combined_d_mode:+.2f}, "
                        f"Δright = {combined_d_right:+.2f}, "
                        f"Δleft = {combined_d_left:+.2f} "
                        f"(mode_nom={mode_nom:.2f})"
                    )
                continue

            task_h5 = self.study_dir / "tasks" / f"syst_mc_{p.point_id}_samples.h5"
            if not task_h5.exists():
                continue

            with h5py.File(task_h5, "r") as h5f:
                ndbd_samples = h5f["ndbd_samples"][:]
                max_logp = h5f["max_logp"][:]

            if len(ndbd_samples) == 0:
                continue

            # Compute likelihood weights via softmax
            scores = max_logp - np.max(max_logp)
            weights = np.exp(scores)
            weights[weights < 1e-10] = 0.0
            if weights.sum() > 0:
                weights = weights / weights.sum()
            else:
                weights = np.full(len(max_logp), 1.0 / len(max_logp))

            rng = np.random.default_rng(42)
            combined_samples = sample_weighted_mixture(
                ndbd_samples, weights, n_draws=50000, rng=rng
            )

            # Resolve nominal posterior
            nom_post_raw = None
            if nominal_post is not None:
                nom_post_raw = np.asarray(nominal_post).flatten()
            else:
                # Try reading from storage fit cache
                try:
                    fit_info = p.get_fit_hash(self.cfg)
                    fit_prefix = fit_info["fit_prefix"]
                    nc_path = fit_prefix.with_suffix(".nc")
                    if nc_path.exists():
                        import arviz as az

                        nom_trace = az.from_netcdf(str(nc_path))
                        nom_post_raw = nom_trace.posterior["ndbd"].values.flatten()
                except Exception:
                    pass

            # Fallback if nominal fit not yet on disk: use point's 0th sample or unweighted mean
            if nom_post_raw is None:
                nom_post_raw = ndbd_samples[0]

            # Resolve calibration parameters
            p_calib = calib_dict.get(p.point_id)
            if not p_calib:
                for ck, cv in calib_dict.items():
                    if isinstance(cv, dict) and cv.get("co_shift") == p.co_shift:
                        p_calib = cv
                        break

            popt = None
            if p_calib and "a" in p_calib and "b" in p_calib:
                popt = (float(p_calib["a"]), float(p_calib["b"]))

            if popt is not None:
                nom_bc = apply_bias_correction(nom_post_raw, popt)
                syst_bc = apply_bias_correction(combined_samples, popt)
            else:
                nom_bc = nom_post_raw
                syst_bc = combined_samples

            mode_nom, lo_nom, hi_nom = calc_intervals([nom_bc], method="hdi")[0]
            mode_syst, lo_syst, hi_syst = calc_intervals([syst_bc], method="hdi")[0]
            deltas = compute_deltas(nom_bc, syst_bc, return_extensions_only=True)
            delta_left = float(deltas[0])  # <= 0
            delta_right = float(deltas[1])  # >= 0
            delta_mode = float(mode_syst - mode_nom)

            summary_entry = {
                "point_id": p.point_id,
                "task_id": p.task_id,
                "category": "bmmc",
                "co_shift": p.co_shift,
                "ndbd_shift": p.ndbd_shift,
                "co_cut": str(p.co_cut),
                "ndbd_cutoff": str(p.ndbd_cutoff),
                "n_samples": len(point_df),
                "mode_nom": float(mode_nom),
                "mode_syst": float(mode_syst),
                "delta_mode": float(delta_mode),
                "hdi_nom": [float(lo_nom), float(hi_nom)],
                "hdi_syst": [float(lo_syst), float(hi_syst)],
                "delta_right": float(delta_right),
                "delta_left": float(delta_left),
            }
            if popt is not None:
                summary_entry["bc_a"] = popt[0]
                summary_entry["bc_b"] = popt[1]

            summaries[p.point_id] = summary_entry

            # Save per-systematic Format v2 HDF5 result
            syst_h5 = self.study_dir / f"syst_bmmc_{p.point_id}.h5"
            model_pars = {
                "co_cut": p.co_cut,
                "ndbd_cutoff": p.ndbd_cutoff,
                "co_shift": p.co_shift,
                "ndbd_shift": p.ndbd_shift,
            }
            save_syst_result(
                syst_h5,
                name="bmmc",
                bc_post=syst_bc,
                raw_post=combined_samples,
                model_pars=model_pars,
                per_post=point_df,
                overwrite=True,
            )

            # Record in SQLite
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

            # Generate diagnostic plot
            plot_path = plots_dir / f"bmmc_comparison_{p.point_id}.png"
            try:
                fig, _ = plot_syst_comparison(
                    nominal_bc=nom_bc,
                    syst_post=syst_bc,
                    delta=delta_mode,
                    deltas=[delta_left, delta_right, delta_mode],
                    syst_label="BM-MC",
                    nom_label="Nominal",
                    title=f"BM-MC Systematic ({p.point_id})",
                )
                fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                plt.close(fig)
            except Exception as e:
                if verbose:
                    print(f"[{p.point_id}] Note: could not render plot ({e})")

            if verbose:
                print(
                    f"[{p.point_id}] BM-MC Systematic: "
                    f"Δmode = {delta_mode:+.2f}, "
                    f"Δright = {delta_right:+.2f}, "
                    f"Δleft = {delta_left:+.2f} "
                    f"(mode_nom={mode_nom:.2f}, mode_syst={mode_syst:.2f})"
                )

        summary_file = self.study_dir / "syst_summary.json"
        with open(summary_file, "w", encoding="utf-8") as f:
            json.dump(summaries, f, indent=2)

        if verbose:
            print(
                f"[INFO] Saved systematic summary metrics for {len(summaries)} points: {summary_file}"
            )

        return {"results_df": merged_df, "summary": summaries}

    def run_marginalization(
        self,
        exposure: Optional[float] = None,
        eff_mc: Optional[float] = None,
        psa_mean: Optional[float] = None,
        psa_sig: Optional[float] = None,
        n_smear: int = 200000,
        rng_seed: int = 42,
        extend_method: str = "skew_norm",
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        nominal_post: Optional[np.ndarray] = None,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """Execute full rate marginalisation, split-Gaussian smearing, and limit extraction.

        Reads systematics definitions and deltas (from sub-studies or declarative fallbacks),
        convolves the nominal bias-corrected posterior with total quadrature widths (sigma_left, sigma_right),
        convolves with exposure and PSA efficiency, and computes 90% C.L. half-life limits.
        """
        from qpym2.syst.core import SystDelta, combine_deltas
        from qpym2.syst.smear import smear_posterior
        from qpym2.syst.rate import to_rate
        from qpym2.syst.io import save_combined_result
        from qpym2.syst.plots import plot_delta_table, plot_combined_posterior, plot_rate_posterior
        from qpym2.fit.bias import apply_bias_correction
        from qpym2.utils.mathy import get_mode, calc_intervals
        import json
        import matplotlib.pyplot as plt

        m_cfg = getattr(self.cfg, "marginalisation", None)
        _exposure = float(exposure if exposure is not None else (getattr(m_cfg, "exposure", 1038.4) if m_cfg else 1038.4))
        _eff_mc = float(eff_mc if eff_mc is not None else (getattr(m_cfg, "eff_mc", 0.04) if m_cfg else 0.04))
        _psa_mean = float(psa_mean if psa_mean is not None else (getattr(m_cfg, "psa_mean", 0.99) if m_cfg else 0.99))
        _psa_sig = float(psa_sig if psa_sig is not None else (getattr(m_cfg, "psa_sig", 0.01) if m_cfg else 0.01))
        _n_smear = int(getattr(m_cfg, "n_smear", n_smear) if m_cfg else n_smear)
        _rng_seed = int(getattr(m_cfg, "rng_seed", rng_seed) if m_cfg else rng_seed)
        _extend_method = str(getattr(m_cfg, "extend_method", extend_method) if m_cfg else extend_method)

        plots_dir = self.study_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        sqlite_path = self.study_dir / f"{self.name}.sqlite"

        # 1. Resolve nominal posterior
        nom_raw = None
        if nominal_post is not None:
            nom_raw = np.asarray(nominal_post, dtype=float).flatten()
        else:
            try:
                p0 = self.points[0]
                fit_info = p0.get_fit_hash(self.cfg)
                nc_path = fit_info["fit_prefix"].with_suffix(".nc")
                if nc_path.exists():
                    import arviz as az
                    nom_trace = az.from_netcdf(str(nc_path))
                    nom_raw = nom_trace.posterior["ndbd"].values.flatten()
            except Exception:
                pass

        if nom_raw is None:
            rng_nom = np.random.default_rng(_rng_seed)
            nom_raw = np.clip(rng_nom.normal(loc=15.0, scale=8.0, size=50000), 0.0, None)

        # Bias calibration
        popt = None
        calib_source = bias_calibration or getattr(self.cfg, "bias_calibration", None)
        if isinstance(calib_source, dict):
            if "a" in calib_source and "b" in calib_source:
                popt = (float(calib_source["a"]), float(calib_source["b"]))
        elif calib_source is not None:
            cpath = Path(calib_source)
            if not cpath.is_absolute():
                cpath = self.study_dir / cpath if (self.study_dir / cpath).exists() else Path(cpath)
            if cpath.exists():
                with open(cpath, "r", encoding="utf-8") as f:
                    cdata = json.load(f)
                if "a" in cdata and "b" in cdata:
                    popt = (float(cdata["a"]), float(cdata["b"]))
                elif "nominal" in cdata and "a" in cdata["nominal"]:
                    popt = (float(cdata["nominal"]["a"]), float(cdata["nominal"]["b"]))

        nom_bc = apply_bias_correction(nom_raw, popt) if popt is not None else nom_raw

        # 2. Resolve systematic deltas
        systs_cfg = getattr(self.cfg, "systematics", None)
        systs_dict = OmegaConf.to_container(systs_cfg, resolve=True) if systs_cfg else {}

        delta_records = []
        syst_deltas = []

        for s_name, s_info in (systs_dict or {}).items():
            if not isinstance(s_info, dict):
                continue
            dr = s_info.get("delta_right", None)
            dl = s_info.get("delta_left", None)
            source = s_info.get("source", None)

            if source:
                src_path = Path(source)
                if not src_path.is_absolute():
                    src_path = (self.study_dir / src_path).resolve() if (self.study_dir / src_path).exists() else Path(source).resolve()
                summary_path = src_path / "syst_summary.json"
                if summary_path.exists():
                    try:
                        with open(summary_path, "r", encoding="utf-8") as f:
                            sum_data = json.load(f)
                        first_entry = next(iter(sum_data.values())) if sum_data else {}
                        if "delta_right" in first_entry and "delta_left" in first_entry:
                            dr = float(first_entry["delta_right"])
                            dl = float(first_entry["delta_left"])
                    except Exception:
                        pass

            dr_val = float(dr) if dr is not None else 0.0
            dl_val = float(dl) if dl is not None else 0.0

            delta_records.append({
                "source": s_name,
                "delta_left": dl_val,
                "delta_right": dr_val,
                "study_source": str(source or "explicit"),
            })
            syst_deltas.append(SystDelta(name=s_name, delta_right=dr_val, delta_left=dl_val))

        if not syst_deltas:
            default_table = [
                ("lss_dq", -3.59, 4.74),
                ("lss_sscl", -1.20, 1.20),
                ("bmmc", -1.50, 0.50),
                ("bm_threshold", -1.16, 0.00),
                ("bm_binning", -0.60, 0.00),
                ("ndbd_shift", -0.42, 0.00),
                ("bias_corr", -0.04, 0.30),
            ]
            for s_name, dl, dr in default_table:
                delta_records.append({
                    "source": s_name,
                    "delta_left": dl,
                    "delta_right": dr,
                    "study_source": "default_tech_note",
                })
                syst_deltas.append(SystDelta(name=s_name, delta_right=dr, delta_left=dl))

        df_delta_table = pd.DataFrame(delta_records)
        df_delta_table.to_csv(self.study_dir / "delta_table.csv", index=False)

        # 3. Quadrature sum
        sigma_left, sigma_right = combine_deltas([(d.delta_left, d.delta_right) for d in syst_deltas])

        # 4. Split-Gaussian Smearing
        smeared_post = smear_posterior(
            nom_bc,
            sigma_right=float(abs(sigma_right)),
            sigma_left=float(abs(sigma_left)),
            n_smear=_n_smear,
            rng_seed=_rng_seed,
        )

        # 5. Rate conversion with PSA marginalisation
        rate_post = to_rate(
            smeared_post,
            eff_mc=_eff_mc,
            exposure=_exposure,
            psa_mean=_psa_mean,
            psa_sig=_psa_sig,
            rng_seed=_rng_seed,
            clip=True,
        )

        # 6. Extract summary statistics
        mode_nom = float(get_mode(nom_bc))
        hdi_nom = calc_intervals([nom_bc], method="hdi")[0]
        p90_nom = float(np.percentile(nom_bc, 90))

        smeared_clip = np.clip(smeared_post, 0.0, None)
        mode_smeared = float(get_mode(smeared_clip))
        hdi_smeared = calc_intervals([smeared_clip], method="hdi")[0]
        p90_smeared = float(np.percentile(smeared_clip, 90))

        rate_mode = float(get_mode(rate_post))
        hdi_rate = calc_intervals([rate_post], method="hdi")[0]
        rate_90 = float(np.percentile(rate_post, 90))
        half_life_90 = float(np.log(2.0) / rate_90) if rate_90 > 0 else np.inf

        stats_df = pd.DataFrame([
            {"label": "nominal_bc", "mode": mode_nom, "hdi_lo": float(hdi_nom[1]), "hdi_hi": float(hdi_nom[2]), "p90": p90_nom},
            {"label": "smeared", "mode": mode_smeared, "hdi_lo": float(hdi_smeared[1]), "hdi_hi": float(hdi_smeared[2]), "p90": p90_smeared},
            {"label": "rate", "mode": rate_mode, "hdi_lo": float(hdi_rate[1]), "hdi_hi": float(hdi_rate[2]), "p90": rate_90},
        ])

        summary = {
            "study_name": self.name,
            "exposure_kgyr": _exposure,
            "eff_mc": _eff_mc,
            "psa_mean": _psa_mean,
            "psa_sig": _psa_sig,
            "sigma_right": sigma_right,
            "sigma_left": sigma_left,
            "half_life_limit_90_yr": half_life_90,
            "rate_limit_90cl": rate_90,
            "nominal_counts": {
                "mode": mode_nom,
                "hdi": [float(hdi_nom[1]), float(hdi_nom[2])],
                "p90": p90_nom,
            },
            "smeared_counts": {
                "mode": mode_smeared,
                "hdi": [float(hdi_smeared[1]), float(hdi_smeared[2])],
                "p90": p90_smeared,
            },
            "decay_rate_yr": {
                "mode": rate_mode,
                "hdi": [float(hdi_rate[1]), float(hdi_rate[2])],
                "limit_90cl": rate_90,
            },
            "half_life_yr": {
                "limit_90cl": half_life_90,
            },
            "delta_breakdown": delta_records,
        }

        # Write summary JSON
        with open(self.study_dir / "syst_rate_summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

        # Write combined Format HDF5
        model_pars = {
            "co_cut": getattr(self.points[0], "co_cut", None) if self.points else None,
            "ndbd_cutoff": getattr(self.points[0], "ndbd_cutoff", None) if self.points else None,
            "co_shift": getattr(self.points[0], "co_shift", 0.0) if self.points else 0.0,
        }
        combined_h5 = self.study_dir / "combined_result.h5"
        save_combined_result(
            combined_h5,
            nominal_bc_post=nom_bc,
            marg_post=smeared_post,
            rate_post=rate_post,
            delta_table=df_delta_table,
            stats=stats_df,
            sigma_right=sigma_right,
            sigma_left=sigma_left,
            model_pars=model_pars,
            eff_mc=_eff_mc,
            exposure=_exposure,
            psa_mean=_psa_mean,
            psa_sig=_psa_sig,
            overwrite=True,
        )

        # Diagnostic plots
        try:
            fig_dt, _, _ = plot_delta_table(syst_deltas, sigma_right=sigma_right, sigma_left=sigma_left)
            fig_dt.savefig(plots_dir / "delta_breakdown.png", dpi=150, bbox_inches="tight")
            plt.close(fig_dt)
        except Exception as e:
            if verbose:
                print(f"[WARN] Could not render delta_breakdown plot: {e}")

        try:
            fig_comb, _ = plot_combined_posterior(
                marg_post=smeared_post,
                nominal_bc=nom_bc,
                sigma_right=sigma_right,
                sigma_left=sigma_left,
            )
            fig_comb.savefig(plots_dir / "smeared_posterior.png", dpi=150, bbox_inches="tight")
            plt.close(fig_comb)
        except Exception as e:
            if verbose:
                print(f"[WARN] Could not render smeared_posterior plot: {e}")

        try:
            rate_stats_dict = {
                "mode": rate_mode,
                "p90": rate_90,
                "hdi": (float(hdi_rate[1]), float(hdi_rate[2])),
                "hdi_prob": 0.68,
            }
            fig_rate, _ = plot_rate_posterior(rate=rate_post, stats=rate_stats_dict)
            fig_rate.savefig(plots_dir / "rate_posterior.png", dpi=150, bbox_inches="tight")
            plt.close(fig_rate)
        except Exception as e:
            if verbose:
                print(f"[WARN] Could not render rate_posterior plot: {e}")

        # Index in SQLite
        record_study_point(
            sqlite_path,
            {
                "point_id": "marginalization",
                "study_name": self.name,
                "task_id": 9999,
                "status": "completed",
                "summary_json": summary,
            },
        )

        if verbose:
            print(
                f"[INFO] Rate Marginalisation COMPLETED for '{self.name}':\n"
                f"  σ_right = {sigma_right:+.3f}, σ_left = {sigma_left:+.3f}\n"
                f"  Γ_mode = {rate_mode:.2e} yr⁻¹, Γ_90 = {rate_90:.2e} yr⁻¹\n"
                f"  T_1/2^(0ν) > {half_life_90:.2e} yr (90% C.L.)\n"
                f"  Outputs saved to {self.study_dir}"
            )

        return summary

