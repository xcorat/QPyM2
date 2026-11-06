"""Shared helpers and utilities for qpym2.studies."""

import contextlib
import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import h5py
import numpy as np
import pandas as pd
from omegaconf import DictConfig

from qpym2.fit.bias import apply_bias_correction
from qpym2.fit.results import _get_stats
from qpym2.utils.mathy import get_mode


def get_task_paths(
    study_dir: Union[str, Path],
    prefix: str,
    point_id: str,
    with_h5: bool = False,
) -> Tuple[Path, Optional[Path]]:
    """Resolve task CSV and optional HDF5 output paths, ensuring parent directory exists."""
    sdir = Path(study_dir)
    tasks_dir = sdir / "tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    task_csv = tasks_dir / f"{prefix}_{point_id}.csv"
    task_h5 = tasks_dir / f"{prefix}_{point_id}_samples.h5" if with_h5 else None
    return task_csv, task_h5


def check_task_status(
    study_dir: Union[str, Path],
    prefix: str,
    point: Any,
    target_count: int,
    target_key: str = "n_samples",
    with_h5: bool = False,
    extra_fields: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Check completion status of a task for a given study point."""
    task_csv, task_h5 = get_task_paths(study_dir, prefix, point.point_id, with_h5=with_h5)
    row_count = 0
    is_done = False
    if task_csv.exists() and task_csv.stat().st_size > 0:
        try:
            df = pd.read_csv(task_csv)
            row_count = len(df)
            if with_h5 and task_h5 is not None:
                is_done = (
                    row_count >= target_count and task_h5.exists() and task_h5.stat().st_size > 0
                )
            else:
                is_done = row_count >= target_count
        except Exception:
            pass

    status = "COMPLETED" if is_done else ("IN_PROGRESS" if row_count > 0 else "PENDING")
    res = {
        "point_id": point.point_id,
        "task_id": point.task_id,
        "model_key": point.model_key,
        "co_shift": point.co_shift,
        "ndbd_shift": point.ndbd_shift,
        "co_cut": str(point.co_cut),
        "ndbd_cutoff": str(point.ndbd_cutoff),
        f"{target_key}_target": target_count,
        f"{target_key}_completed": row_count,
        "is_completed": is_done,
        "task_csv": str(task_csv),
        "status": status,
    }
    if with_h5 and task_h5 is not None:
        res["task_h5"] = str(task_h5)
    if extra_fields:
        res.update(extra_fields)
    return res


def load_bias_popt(
    bias_calibration: Optional[Union[str, Path, Dict[str, Any]]],
    study_dir: Union[str, Path],
    point_id: Optional[str] = None,
    co_shift: Optional[float] = None,
    dedicated_bias_csv: Optional[Union[str, Path]] = None,
    cfg: Optional[DictConfig] = None,
) -> Optional[Tuple[float, float]]:
    """Resolve linear bias calibration coefficients (a, b) for a study point."""
    sdir = Path(study_dir)

    # 1. Dedicated bias CSV override
    if dedicated_bias_csv:
        dpath = Path(dedicated_bias_csv)
        if not dpath.is_absolute():
            dpath = sdir / dpath if (sdir / dpath).exists() else Path(dpath)
        if dpath.exists():
            from qpym2.fit.bias import get_bias_from_fitdf

            try:
                vbias_df = pd.read_csv(str(dpath))
                popt, _ = get_bias_from_fitdf(vbias_df)
                if popt is not None:
                    return (float(popt[0]), float(popt[1]))
            except Exception:
                pass

    # 2. Provided calibration dictionary or path
    calib_dict: Optional[Dict[str, Any]] = None
    if isinstance(bias_calibration, dict):
        calib_dict = bias_calibration
    elif bias_calibration is not None:
        cpath = Path(bias_calibration)
        if not cpath.is_absolute():
            cpath = sdir / cpath if (sdir / cpath).exists() else Path(cpath)
        if cpath.exists():
            try:
                with open(cpath, "r", encoding="utf-8") as f:
                    calib_dict = json.load(f)
            except Exception:
                pass

    # 3. Fallback to study configuration or local bias_calibration.json
    if calib_dict is None and cfg is not None:
        scfg = getattr(cfg, "syst", None)
        cfg_calib = getattr(scfg, "bias_calibration", None) if scfg else None
        if not cfg_calib:
            cfg_calib = getattr(cfg, "bias_calibration", None)
        if isinstance(cfg_calib, dict):
            calib_dict = cfg_calib
        elif cfg_calib:
            cpath = Path(str(cfg_calib))
            if not cpath.is_absolute():
                cpath = sdir / cpath if (sdir / cpath).exists() else Path(cpath)
            if cpath.exists():
                try:
                    with open(cpath, "r", encoding="utf-8") as f:
                        calib_dict = json.load(f)
                except Exception:
                    pass

    if calib_dict is None and (sdir / "bias_calibration.json").exists():
        try:
            with open(sdir / "bias_calibration.json", "r", encoding="utf-8") as f:
                calib_dict = json.load(f)
        except Exception:
            pass

    if not calib_dict:
        return None

    # Resolve point calibration from dict
    p_calib = None
    if point_id and point_id in calib_dict:
        p_calib = calib_dict[point_id]
    elif co_shift is not None:
        for cv in calib_dict.values():
            if isinstance(cv, dict) and cv.get("co_shift") == co_shift:
                p_calib = cv
                break

    if p_calib is None:
        if "a" in calib_dict and "b" in calib_dict:
            p_calib = calib_dict
        elif "nominal" in calib_dict and isinstance(calib_dict["nominal"], dict):
            p_calib = calib_dict["nominal"]

    if p_calib and "a" in p_calib and "b" in p_calib:
        return (float(p_calib["a"]), float(p_calib["b"]))

    return None


def patch_pymc_parallel() -> None:
    """Disable threadpool limits in pymc.sampling.parallel to avoid fork deadlocks in threads."""
    try:
        import pymc.sampling.parallel as _parallel

        _parallel.threadpool_limits = lambda **kwargs: contextlib.nullcontext()
    except Exception:
        pass


def extract_logp_mode(trace: Any) -> float:
    """Extract mode of log-likelihood (pik) from PyMC trace InferenceData."""
    if not hasattr(trace, "log_likelihood") or trace.log_likelihood is None:
        return 0.0
    ll = trace.log_likelihood
    if hasattr(ll, "pik"):
        arr = np.asarray(ll.pik.values).flatten()
    elif isinstance(ll, dict) and "pik" in ll:
        arr = np.asarray(ll["pik"]).flatten()
    else:
        arr = np.array([])

    if len(arr) > 1:
        try:
            return float(get_mode(arr))
        except Exception:
            return float(np.median(arr))
    elif len(arr) == 1:
        return float(arr[0])
    return 0.0


def extract_ndbd_stats(
    ndbd_post: np.ndarray,
    popt: Optional[Tuple[float, float]] = None,
) -> Dict[str, float]:
    """Calculate summary statistics and optional bias-corrected metrics for ndbd posterior."""
    stats = _get_stats(ndbd_post)
    mode_val = float(stats["mode"])
    p16_val = float(stats.get("p16", np.percentile(ndbd_post, 16)))
    p90_val = float(stats.get("p90", np.percentile(ndbd_post, 90)))
    upp34_val = float(stats.get("upp34", mode_val))

    mode_bc = p16_bc = p90_bc = np.nan
    if popt is not None:
        corrected = apply_bias_correction(ndbd_post, popt)
        mode_bc = float(get_mode(corrected))
        p16_bc = float(np.percentile(corrected, 16))
        p90_bc = float(np.percentile(corrected, 90))

    return {
        "mode": mode_val,
        "p16": p16_val,
        "p90": p90_val,
        "upp34": upp34_val,
        "mode_bc": mode_bc,
        "p16_bc": p16_bc,
        "p90_bc": p90_bc,
    }


def write_task_csv_atomic(df: pd.DataFrame, target_path: Path) -> None:
    """Atomically write DataFrame to CSV via a temporary file."""
    temp_csv = target_path.parent / f".tmp_{target_path.name}"
    df.to_csv(temp_csv, index=False)
    temp_csv.replace(target_path)


def write_task_h5_atomic(
    target_path: Path,
    datasets: Dict[str, np.ndarray],
    attrs: Optional[Dict[str, Any]] = None,
) -> None:
    """Atomically write datasets and attributes to HDF5 via a temporary file."""
    temp_h5 = target_path.parent / f".tmp_{target_path.name}"
    with h5py.File(temp_h5, "w") as h5f:
        for dname, dval in datasets.items():
            if dval.ndim > 1:
                h5f.create_dataset(dname, data=dval, compression="gzip", compression_opts=4)
            else:
                h5f.create_dataset(dname, data=dval)
        if attrs:
            for ak, av in attrs.items():
                if av is not None:
                    h5f.attrs[ak] = av
    temp_h5.replace(target_path)


def resolve_nominal_posterior(
    point: Any,
    cfg: DictConfig,
    nominal_post: Optional[np.ndarray] = None,
) -> Optional[np.ndarray]:
    """Resolve nominal posterior array from override or cached fit .nc file."""
    if nominal_post is not None:
        return np.asarray(nominal_post, dtype=float).flatten()

    try:
        fit_info = point.get_fit_hash(cfg)
        fit_prefix = fit_info[2] if isinstance(fit_info, tuple) else fit_info.get("fit_prefix")
        if fit_prefix is not None:
            nc_path = Path(fit_prefix).with_suffix(".nc")
            if nc_path.exists():
                import arviz as az

                nom_trace = az.from_netcdf(str(nc_path))
                return nom_trace.posterior["ndbd"].values.flatten()
    except Exception:
        pass
    return None
