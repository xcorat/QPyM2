"""Bias and sensitivity toy simulation and results collection."""

import json
from pathlib import Path
from typing import Any, Dict, Optional, Union

import numpy as np
import pandas as pd
from omegaconf import DictConfig

from qpym2.fit.bias import get_bias_from_fitdf
from qpym2.fit.results import _get_stats
from qpym2.io.toys import _unflatten_columns
from qpym2.models import m2
from qpym2.storage import fits as storage_fits
from qpym2.storage.study import collect_study_results, record_study_point
from qpym2.studies._common import (
    check_task_status,
    extract_logp_mode,
    get_task_paths,
    patch_pymc_parallel,
    write_task_csv_atomic,
)
from qpym2.utils.hists import get_sum_hist, sample_hist


def check_bias_toys(point: Any, study_dir: Union[str, Path], n_toys: int = 200) -> Dict[str, Any]:
    """Check completion status of bias toy fits for this study point."""
    return check_task_status(
        study_dir=study_dir,
        prefix="toys_bias",
        point=point,
        target_count=n_toys,
        target_key="n_toys",
        with_h5=False,
    )


def check_sens_toys(point: Any, study_dir: Union[str, Path], n_toys: int = 200) -> Dict[str, Any]:
    """Check completion status of sensitivity toy fits for this study point."""
    return check_task_status(
        study_dir=study_dir,
        prefix="toys_sens",
        point=point,
        target_count=n_toys,
        target_key="n_toys",
        with_h5=False,
    )


def run_bias_toys(
    point: Any,
    cfg: DictConfig,
    study_dir: Union[str, Path],
    n_toys: int = 200,
    bias_min: float = 5.0,
    bias_max: float = 200.0,
    seed: int = 42,
    chains: int = 4,
    draws: int = 2000,
    tune: int = 500,
    target_accept: float = 0.9,
    single_core: bool = True,
    prior_width: str = "5sig",
    mkchain_samples: int = 100,
    fit_fn: Optional[Any] = None,
    force: bool = False,
    verbose: bool = True,
) -> Path:
    """Execute a batch of bias toy fits for this study point."""
    task_csv, _ = get_task_paths(study_dir, "toys_bias", point.point_id, with_h5=False)

    if not force:
        status = check_bias_toys(point, study_dir, n_toys=n_toys)
        if status["is_completed"]:
            if verbose:
                print(
                    f"[{point.point_id}] Bias toys already completed ({status['n_toys_completed']} rows). Skipping."
                )
            return task_csv

    fit_setup = storage_fits.setup_storage_fit(
        cfg=cfg,
        model_key=point.model_key,
        co_cut=point.co_cut,
        ndbd_cutoff=point.ndbd_cutoff,
        co_shift=point.co_shift,
        ndbd_shift=point.ndbd_shift,
        prior_width=prior_width,
        mkchain_samples=mkchain_samples,
        chains=chains,
        draws=draws,
        tune=tune,
        target_accept=target_accept,
        dry_run=False,
    )

    fit_comps = fit_setup["fit_comps"]
    xedges, yedges = fit_setup["edges"]
    component_names = list(fit_comps.index)
    if "ndbd" not in component_names:
        raise ValueError(f"Signal component 'ndbd' not found in components: {component_names}")
    isig = component_names.index("ndbd")

    rng = np.random.default_rng(seed + point.task_id * 10000)
    injections = rng.uniform(bias_min, bias_max, size=n_toys)
    base_norms = np.array(
        [
            float(fit_comps.loc[c, "integral"]) if "integral" in fit_comps.columns else 1.0
            for c in component_names
        ]
    )

    patch_pymc_parallel()
    import pymc as pm

    rows = []
    for i, inj in enumerate(injections):
        norms = base_norms.copy()
        norms[isig] = float(inj)
        hsum = get_sum_hist(fit_comps, norms)
        data_i = sample_hist(hsum, xedges, yedges)

        if fit_fn is not None:
            trace = fit_fn(data_i, fit_comps, xedges, yedges)
        else:
            fik = m2.make_fik(fit_comps, data_i, xedges, yedges)
            model = m2.m2_model(fit_comps, fik)
            with model:
                trace = pm.sample(
                    draws=draws,
                    tune=tune,
                    chains=chains,
                    cores=1 if single_core else chains,
                    return_inferencedata=True,
                    target_accept=target_accept,
                    progressbar=False,
                    idata_kwargs={"log_likelihood": True},
                    random_seed=int(seed + point.task_id * 10000 + i + 1),
                )

        ndbd_post = trace.posterior["ndbd"].values.flatten()
        stats = _get_stats(ndbd_post)
        mode_val = float(stats["mode"])
        upp34_val = float(stats["upp34"])
        uerr = upp34_val - mode_val
        pull = (mode_val - float(inj)) / uerr if abs(uerr) > 1e-9 else 0.0

        logp_val = extract_logp_mode(trace)

        row = {
            "toy_id": i,
            "task_id": point.task_id,
            "point_id": point.point_id,
            "model_key": point.model_key,
            "co_shift": point.co_shift,
            "ndbd_shift": point.ndbd_shift,
            "co_cut": str(point.co_cut),
            "ndbd_cutoff": str(point.ndbd_cutoff),
            "type_": "bias",
            "ndbd_inj_": float(inj),
            "logp_mode_": logp_val,
            "pull": float(pull),
        }

        for var in component_names:
            var_post = trace.posterior[var].values.flatten()
            v_stats = _get_stats(var_post)
            for sname, sval in v_stats.items():
                row[f"{var}_{sname}"] = float(sval)

        rows.append(row)
        if verbose and ((i + 1) % 10 == 0 or (i + 1) == n_toys):
            print(f"[{point.point_id}] Completed {i + 1}/{n_toys} bias toys", flush=True)

    df_res = pd.DataFrame(rows)
    write_task_csv_atomic(df_res, task_csv)

    if verbose:
        print(f"[{point.point_id}] Saved {len(df_res)} bias toys: {task_csv}")

    return task_csv


def run_sens_toys(
    point: Any,
    cfg: DictConfig,
    study_dir: Union[str, Path],
    n_toys: int = 200,
    seed: int = 42,
    chains: int = 4,
    draws: int = 2000,
    tune: int = 500,
    target_accept: float = 0.9,
    single_core: bool = True,
    prior_width: str = "5sig",
    mkchain_samples: int = 100,
    fit_fn: Optional[Any] = None,
    force: bool = False,
    verbose: bool = True,
) -> Path:
    """Execute a batch of sensitivity toy fits (background-only hypothesis) for this study point."""
    task_csv, _ = get_task_paths(study_dir, "toys_sens", point.point_id, with_h5=False)

    if not force:
        status = check_sens_toys(point, study_dir, n_toys=n_toys)
        if status["is_completed"]:
            if verbose:
                print(
                    f"[{point.point_id}] Sensitivity toys already completed ({status['n_toys_completed']} rows). Skipping."
                )
            return task_csv

    fit_setup = storage_fits.setup_storage_fit(
        cfg=cfg,
        model_key=point.model_key,
        co_cut=point.co_cut,
        ndbd_cutoff=point.ndbd_cutoff,
        co_shift=point.co_shift,
        ndbd_shift=point.ndbd_shift,
        prior_width=prior_width,
        mkchain_samples=mkchain_samples,
        chains=chains,
        draws=draws,
        tune=tune,
        target_accept=target_accept,
        dry_run=False,
    )

    fit_comps = fit_setup["fit_comps"]
    xedges, yedges = fit_setup["edges"]
    component_names = list(fit_comps.index)
    if "ndbd" not in component_names:
        raise ValueError(f"Signal component 'ndbd' not found in components: {component_names}")
    isig = component_names.index("ndbd")

    base_norms = np.array(
        [
            float(fit_comps.loc[c, "integral"]) if "integral" in fit_comps.columns else 1.0
            for c in component_names
        ]
    )

    norms = base_norms.copy()
    norms[isig] = 0.0

    patch_pymc_parallel()
    import pymc as pm

    hsum = get_sum_hist(fit_comps, norms)

    rows = []
    for i in range(n_toys):
        toy_seed = int(seed + point.task_id * 10000 + i + 1)
        np.random.seed(toy_seed)

        data_i = sample_hist(hsum, xedges, yedges)

        if fit_fn is not None:
            trace = fit_fn(data_i, fit_comps, xedges, yedges)
        else:
            fik = m2.make_fik(fit_comps, data_i, xedges, yedges)
            model = m2.m2_model(fit_comps, fik)
            with model:
                trace = pm.sample(
                    draws=draws,
                    tune=tune,
                    chains=chains,
                    cores=1 if single_core else chains,
                    return_inferencedata=True,
                    target_accept=target_accept,
                    progressbar=False,
                    idata_kwargs={"log_likelihood": True},
                    random_seed=toy_seed,
                )

        ndbd_post = trace.posterior["ndbd"].values.flatten()
        stats = _get_stats(ndbd_post)
        mode_val = float(stats["mode"])
        upp34_val = float(stats["upp34"])
        uerr = upp34_val - mode_val
        pull = (mode_val - 0.0) / uerr if abs(uerr) > 1e-9 else 0.0

        logp_val = extract_logp_mode(trace)

        row = {
            "toy_id": i,
            "task_id": point.task_id,
            "point_id": point.point_id,
            "model_key": point.model_key,
            "co_shift": point.co_shift,
            "ndbd_shift": point.ndbd_shift,
            "co_cut": str(point.co_cut),
            "ndbd_cutoff": str(point.ndbd_cutoff),
            "type_": "sens",
            "ndbd_inj_": 0.0,
            "logp_mode_": logp_val,
            "pull": float(pull),
        }

        for var in component_names:
            var_post = trace.posterior[var].values.flatten()
            v_stats = _get_stats(var_post)
            for sname, sval in v_stats.items():
                row[f"{var}_{sname}"] = float(sval)

        rows.append(row)
        if verbose and ((i + 1) % 10 == 0 or (i + 1) == n_toys):
            print(f"[{point.point_id}] Completed {i + 1}/{n_toys} sensitivity toys", flush=True)

    df_res = pd.DataFrame(rows)
    write_task_csv_atomic(df_res, task_csv)

    if verbose:
        print(f"[{point.point_id}] Saved {len(df_res)} sensitivity toys: {task_csv}")

    return task_csv


def collect_bias_results(
    study: Any,
    output_csv: Optional[Union[str, Path]] = None,
    xmin: Optional[float] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Merge bias task CSVs into results.csv and compute calibration parameters via get_bias_from_fitdf."""
    bcfg = getattr(study.cfg, "bias", None)
    _xmin = xmin if xmin is not None else (float(getattr(bcfg, "xmin", 15.0)) if bcfg else 15.0)

    sqlite_path = study.study_dir / f"{study.name}.sqlite"
    merged_df = collect_study_results(study.study_dir, study.name)
    if merged_df.empty:
        if verbose:
            print(f"[WARN] No task CSV results found in {study.study_dir}/tasks.")
        return {"results_df": merged_df, "calibration": {}}

    calibrations: Dict[str, Any] = {}
    for p in study.points:
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

            record_study_point(
                sqlite_path,
                {
                    "point_id": p.point_id,
                    "study_name": study.name,
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

    calib_file = study.study_dir / "bias_calibration.json"
    with open(calib_file, "w", encoding="utf-8") as f:
        json.dump(calibrations, f, indent=2)

    if verbose:
        print(
            f"[INFO] Saved bias calibration parameters for {len(calibrations)} points: {calib_file}"
        )

    return {"results_df": merged_df, "calibration": calibrations}


def collect_sens_results(
    study: Any,
    output_csv: Optional[Union[str, Path]] = None,
    bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Merge sensitivity task CSVs into results.csv and compute summary quantiles per point."""
    sqlite_path = study.study_dir / f"{study.name}.sqlite"
    merged_df = collect_study_results(study.study_dir, study.name)
    if merged_df.empty:
        if verbose:
            print(f"[WARN] No sensitivity task CSV results found in {study.study_dir}/tasks.")
        return {"results_df": merged_df, "summary": {}}

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
        scfg = getattr(study.cfg, "sens", None)
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

        p90_vals = (
            point_df["ndbd_p90"].dropna().values if "ndbd_p90" in point_df.columns else np.array([])
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

        p_calib = calib_dict.get(p.point_id)
        if not p_calib:
            for ck, cv in calib_dict.items():
                if isinstance(cv, dict) and cv.get("co_shift") == p.co_shift:
                    p_calib = cv
                    break

        if p_calib and "a" in p_calib and "b" in p_calib:
            a_val = float(p_calib["a"])
            b_val = float(p_calib["b"])
            p90_bc = np.clip((p90_vals - b_val) / a_val, 0.0, None)
            summary_entry["bc_a"] = a_val
            summary_entry["bc_b"] = b_val
            summary_entry["sens_bc_p90_median"] = float(np.median(p90_bc))
            summary_entry["sens_bc_p90_p16"] = float(np.percentile(p90_bc, 16))
            summary_entry["sens_bc_p90_p84"] = float(np.percentile(p90_bc, 84))
            summary_entry["sens_bc_p90_p025"] = float(np.percentile(p90_bc, 2.5))
            summary_entry["sens_bc_p90_p975"] = float(np.percentile(p90_bc, 97.5))

        summaries[p.point_id] = summary_entry

        record_study_point(
            sqlite_path,
            {
                "point_id": p.point_id,
                "study_name": study.name,
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

    summary_file = study.study_dir / "sens_summary.json"
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)

    if verbose:
        print(
            f"[INFO] Saved sensitivity summary metrics for {len(summaries)} points: {summary_file}"
        )

    return {"results_df": merged_df, "summary": summaries}
