"""Systematics campaign runner, status reporting, and result aggregation dispatcher."""

import json
from pathlib import Path
from typing import Any, Dict, Optional, Union

import pandas as pd

from qpym2.storage.study import collect_study_results
from qpym2.studies import syst_bm, syst_lss, syst_mc, syst_template


def _resolve_category(cfg: Any) -> str:
    """Normalize systematic category and optional subcommand from configuration."""
    scfg = getattr(cfg, "syst", None)
    category = str(getattr(scfg, "category", "bmmc")) if scfg else "bmmc"
    subcmd = getattr(scfg, "subcommand", None) if scfg else None
    if category == "lss" and subcmd:
        category = f"lss_{subcmd}"
    return category


def check_syst(study: Any, n_samples: Optional[int] = None) -> pd.DataFrame:
    """Check systematic task completion status across all study points."""
    scfg = getattr(study.cfg, "syst", None)
    category = _resolve_category(study.cfg)

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
    for p in study.points:
        if category in ("lss_sscl", "lss_dq"):
            records.append(
                syst_lss.check_syst_lss(
                    p,
                    study.study_dir,
                    category=category,
                    datasets=datasets,
                    n_samples=target_n,
                    dq_mode=dq_mode,
                )
            )
        elif category in ("bmmc", "mc_stats", "bm_norm"):
            records.append(syst_mc.check_syst_mc(p, study.study_dir, n_samples=target_n))
        elif category == "bm_threshold":
            _var = getattr(scfg, "variant", "Threshold250") if scfg else "Threshold250"
            records.append(syst_bm.check_syst_bm_threshold(p, study.study_dir, variant=_var))
        elif category == "bm_binning":
            _vars = getattr(scfg, "variants", None) if scfg else None
            records.append(syst_bm.check_syst_bm_binning(p, study.study_dir, variants=_vars))
        elif category in ("bin_size", "template_binning"):
            _sizes = getattr(scfg, "bin_sizes", None) if scfg else None
            records.append(
                syst_template.check_syst_template_binning(p, study.study_dir, bin_sizes=_sizes)
            )
        else:
            records.append(syst_mc.check_syst_mc(p, study.study_dir, n_samples=target_n))

    return pd.DataFrame(records)


def print_syst_status(
    study: Any, df: Optional[pd.DataFrame] = None, n_samples: Optional[int] = None
) -> None:
    """Print human-readable summary of systematic task completion across study points."""
    status_df = df if df is not None else check_syst(study, n_samples=n_samples)
    total_points = len(status_df)
    completed_points = sum(status_df["is_completed"]) if "is_completed" in status_df else 0
    category = _resolve_category(study.cfg)

    print("\n" + "=" * 70)
    print(f"Systematic Study Tasks [{category}]: {study.name} ({total_points} total points)")
    print(f"Directory: {study.study_dir}")
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
    study: Any,
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
    """Execute systematic fits across all points in the study."""
    scfg = getattr(study.cfg, "syst", None)

    def _resolve(val, key, default):
        if val is not None:
            return val
        if scfg and hasattr(scfg, key):
            return getattr(scfg, key)
        return default

    category = _resolve_category(study.cfg)
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
    for p in study.points:
        if not force:
            if category in ("lss_sscl", "lss_dq"):
                st = syst_lss.check_syst_lss(
                    p,
                    study.study_dir,
                    category=category,
                    datasets=_datasets,
                    n_samples=_n_samples,
                    dq_mode=_dq_mode,
                )
            elif category == "bm_threshold":
                _var = getattr(scfg, "variant", "Threshold250") if scfg else "Threshold250"
                st = syst_bm.check_syst_bm_threshold(p, study.study_dir, variant=_var)
            elif category == "bm_binning":
                _vars = getattr(scfg, "variants", None) if scfg else None
                st = syst_bm.check_syst_bm_binning(p, study.study_dir, variants=_vars)
            elif category in ("bin_size", "template_binning"):
                _sizes = getattr(scfg, "bin_sizes", None) if scfg else None
                st = syst_template.check_syst_template_binning(p, study.study_dir, bin_sizes=_sizes)
            else:
                st = syst_mc.check_syst_mc(p, study.study_dir, n_samples=_n_samples)

            if st["is_completed"]:
                if verbose:
                    print(
                        f"[{p.point_id}] Systematic fits for '{category}' already completed. Skipping."
                    )
                results[p.point_id] = (Path(st["task_csv"]), Path(st["task_h5"]))
                continue

        if category in ("lss_sscl", "lss_dq"):
            _lss_path = getattr(scfg, "lss_data_path", None) if scfg else None
            res = syst_lss.run_syst_lss(
                p,
                study.cfg,
                study_dir=study.study_dir,
                category=category,
                datasets=_datasets,
                dq_mode=_dq_mode,
                n_samples=_n_samples,
                lss_data_path=_lss_path,
                seed=_seed,
                chains=_chains,
                draws=_draws,
                tune=_tune,
                target_accept=_target_accept,
                bias_calibration=_bias_calib,
                fit_fn=fit_fn,
                save_traces=save_traces,
                force=force,
                verbose=verbose,
            )
        elif category == "bm_threshold":
            _var = getattr(scfg, "variant", "Threshold250") if scfg else "Threshold250"
            _mk_syst = getattr(scfg, "model_key_syst", "tech_note_250") if scfg else "tech_note_250"
            _raw = getattr(scfg, "results_raw_path", None) if scfg else None
            _dbias = getattr(scfg, "dedicated_bias_csv", None) if scfg else None
            res = syst_bm.run_syst_bm_threshold(
                p,
                study.cfg,
                study_dir=study.study_dir,
                variant=_var,
                model_key_syst=_mk_syst,
                results_raw_path=_raw,
                seed=_seed,
                chains=_chains,
                draws=_draws,
                tune=_tune,
                target_accept=_target_accept,
                bias_calibration=_bias_calib,
                dedicated_bias_csv=_dbias,
                fit_fn=fit_fn,
                save_traces=save_traces,
                force=force,
                verbose=verbose,
            )
        elif category == "bm_binning":
            _vars = getattr(scfg, "variants", None) if scfg else None
            _raw = getattr(scfg, "results_raw_path", None) if scfg else None
            res = syst_bm.run_syst_bm_binning(
                p,
                study.cfg,
                study_dir=study.study_dir,
                variants=_vars,
                results_raw_path=_raw,
                seed=_seed,
                chains=_chains,
                draws=_draws,
                tune=_tune,
                target_accept=_target_accept,
                bias_calibration=_bias_calib,
                fit_fn=fit_fn,
                save_traces=save_traces,
                force=force,
                verbose=verbose,
            )
        elif category in ("bin_size", "template_binning"):
            _sizes = getattr(scfg, "bin_sizes", None) if scfg else None
            res = syst_template.run_syst_template_binning(
                p,
                study.cfg,
                study_dir=study.study_dir,
                bin_sizes=_sizes,
                seed=_seed,
                chains=_chains,
                draws=_draws,
                tune=_tune,
                target_accept=_target_accept,
                bias_calibration=_bias_calib,
                fit_fn=fit_fn,
                save_traces=save_traces,
                force=force,
                verbose=verbose,
            )
        else:
            res = syst_mc.run_syst_mc(
                p,
                study.cfg,
                study_dir=study.study_dir,
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
        results[p.point_id] = res

    return results


def collect_syst_results(
    study: Any,
    bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
    nominal_post: Optional[Any] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Merge systematic task CSVs, compute deltas against nominal, and record summary."""
    merged_df = collect_study_results(study.study_dir, study.name)
    if merged_df.empty:
        if verbose:
            print(f"[WARN] No systematic task CSV results found in {study.study_dir}/tasks.")
        return {"results_df": merged_df, "summary": {}}

    scfg = getattr(study.cfg, "syst", None)
    category = _resolve_category(study.cfg)

    calib_dict: Dict[str, Any] = {}
    if isinstance(bias_calibration, dict):
        calib_dict = bias_calibration
    elif bias_calibration is not None:
        cpath = Path(bias_calibration)
        if not cpath.is_absolute():
            cpath = study.study_dir / cpath if (study.study_dir / cpath).exists() else Path(cpath)
        if cpath.exists():
            with open(cpath, "r", encoding="utf-8") as f:
                calib_dict = json.load(f)
    else:
        cfg_calib = getattr(scfg, "bias_calibration", None) if scfg else None
        if cfg_calib:
            cpath = Path(str(cfg_calib))
            if not cpath.is_absolute():
                cpath = (
                    study.study_dir / cpath if (study.study_dir / cpath).exists() else Path(cpath)
                )
            if cpath.exists():
                with open(cpath, "r", encoding="utf-8") as f:
                    calib_dict = json.load(f)
        elif (study.study_dir / "bias_calibration.json").exists():
            with open(study.study_dir / "bias_calibration.json", "r", encoding="utf-8") as f:
                calib_dict = json.load(f)

    summaries: Dict[str, Any] = {}
    for p in study.points:
        point_df = merged_df[merged_df["point_id"] == p.point_id].copy()
        if point_df.empty:
            continue

        if category in ("lss_sscl", "lss_dq"):
            s_entry = syst_lss.collect_syst_lss(
                study, p, point_df, category, calib_dict, nominal_post, verbose=verbose
            )
        elif category == "bm_threshold":
            s_entry = syst_bm.collect_syst_bm_threshold(
                study, p, point_df, calib_dict, nominal_post, verbose=verbose
            )
        elif category == "bm_binning":
            s_entry = syst_bm.collect_syst_bm_binning(
                study, p, point_df, calib_dict, nominal_post, verbose=verbose
            )
        elif category in ("bin_size", "template_binning"):
            s_entry = syst_template.collect_syst_template_binning(
                study, p, point_df, calib_dict, nominal_post, verbose=verbose
            )
        else:
            s_entry = syst_mc.collect_syst_mc(
                study, p, point_df, calib_dict, nominal_post, verbose=verbose
            )

        if s_entry:
            summaries[p.point_id] = s_entry

    summary_file = study.study_dir / "syst_summary.json"
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)

    if verbose:
        print(
            f"[INFO] Saved systematic summary metrics for {len(summaries)} points: {summary_file}"
        )

    return {"results_df": merged_df, "summary": summaries}
