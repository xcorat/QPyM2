"""Background model variant systematics (threshold and binning) simulation and collection."""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from omegaconf import DictConfig

from qpym2.cfg import config_builders
from qpym2.fit.bias import apply_bias_correction
from qpym2.io import bkg_model
from qpym2.models import m2
from qpym2.storage import fits as storage_fits
from qpym2.storage.study import record_study_point
from qpym2.studies._common import (
    check_task_status,
    extract_logp_mode,
    extract_ndbd_stats,
    get_task_paths,
    load_bias_popt,
    patch_pymc_parallel,
    resolve_nominal_posterior,
    write_task_csv_atomic,
    write_task_h5_atomic,
)
from qpym2.syst.core import compute_deltas
from qpym2.syst.io import save_syst_result
from qpym2.syst.plots import plot_syst_comparison
from qpym2.utils.mathy import calc_intervals, get_mode


def check_syst_bm_threshold(
    point: Any,
    study_dir: Union[str, Path],
    variant: str = "Threshold250",
) -> Dict[str, Any]:
    """Check completion status of BM threshold systematic fit for this study point."""
    return check_task_status(
        study_dir=study_dir,
        prefix="syst_bm_threshold",
        point=point,
        target_count=1,
        target_key="n_samples",
        with_h5=True,
        extra_fields={"category": "bm_threshold", "variant": variant},
    )


def run_syst_bm_threshold(
    point: Any,
    cfg: DictConfig,
    study_dir: Union[str, Path],
    variant: str = "Threshold250",
    model_key_syst: str = "tech_note_250",
    results_raw_path: Optional[Union[str, Path]] = None,
    seed: int = 42,
    chains: int = 2,
    draws: int = 2500,
    tune: int = 500,
    target_accept: float = 0.9,
    single_core: bool = True,
    bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
    dedicated_bias_csv: Optional[Union[str, Path]] = None,
    fit_fn: Optional[Any] = None,
    save_traces: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> Tuple[Path, Path]:
    """Execute BM threshold systematic fit (e.g. Threshold250 vs Threshold200)."""
    task_csv, task_h5 = get_task_paths(study_dir, "syst_bm_threshold", point.point_id, with_h5=True)
    assert task_h5 is not None

    if not force:
        status = check_syst_bm_threshold(point, study_dir, variant=variant)
        if status["is_completed"]:
            if verbose:
                print(f"[{point.point_id}] BM threshold fit already completed. Skipping.")
            return task_csv, task_h5

    raw_path = results_raw_path or config_builders.get_results_raw_path(cfg)
    if not raw_path or not Path(raw_path).exists():
        default_path = Path("/data/cuore/Analysis/BMMC/inputs/ResultsRaw.h5")
        if default_path.exists():
            raw_path = default_path
        else:
            raise FileNotFoundError(
                f"ResultsRaw.h5 not found at '{raw_path}'. Configure infrastructure.data_files.results_raw."
            )

    bm_modes = bkg_model.read_bm_variant_modes(raw_path, variant)

    coords_syst, hm_syst, _, input_dir_syst = storage_fits.resolve_storage_coords(
        cfg, model_key=model_key_syst
    )
    comp_cfg_raw_syst = [
        dict(c) for c in config_builders.get_comp_cfg(cfg, model_key=model_key_syst)
    ]

    fit_setup_syst = storage_fits.setup_storage_fit(
        cfg=cfg,
        model_key=model_key_syst,
        co_cut=point.co_cut,
        ndbd_cutoff=point.ndbd_cutoff,
        co_shift=point.co_shift,
        ndbd_shift=point.ndbd_shift,
        blinded=point.blinded,
        prior_width="5sig",
        mkchain_samples=100,
        chains=chains,
        draws=draws,
        tune=tune,
        target_accept=target_accept,
        dry_run=False,
    )

    data_syst = fit_setup_syst["data"]
    xedges_syst, yedges_syst = fit_setup_syst["edges"]
    nominal_comps_syst = fit_setup_syst["fit_comps"]
    nominal_priors_syst = list(nominal_comps_syst["prior"])

    tpl_uuids, uuid_base, uuid_co, uuid_ndbd = storage_fits.resolve_fit_template_uuids(
        cfg=cfg,
        model_key=model_key_syst,
        co_cut=point.co_cut,
        ndbd_cutoff=point.ndbd_cutoff,
        co_shift=point.co_shift,
        ndbd_shift=point.ndbd_shift,
    )

    mctable_syst, signal_df_syst, _, _ = storage_fits.assemble_fit_components(
        input_dir=input_dir_syst,
        jags_fpath=config_builders.get_jags_path(cfg),
        comp_cfg=comp_cfg_raw_syst,
        hm=hm_syst,
        uuid_base=uuid_base,
        uuid_co=uuid_co,
        uuid_ndbd=uuid_ndbd,
    )

    mt_syst = mctable_syst.copy()
    mt_syst["norm"] = bm_modes.reindex(mt_syst.index).fillna(0.0)

    bkg_comps_syst = m2.make_bkgcomps_table(comp_cfg_raw_syst, mt_syst, smooth=hm_syst.smooth)
    comps_syst = pd.concat([bkg_comps_syst, signal_df_syst], join="inner")
    comps_syst["prior"] = nominal_priors_syst
    comps_syst["varname"] = comps_syst.index

    popt = load_bias_popt(
        bias_calibration=bias_calibration,
        study_dir=study_dir,
        point_id=point.point_id,
        co_shift=point.co_shift,
        dedicated_bias_csv=dedicated_bias_csv,
        cfg=cfg,
    )

    traces_dir = Path(study_dir) / "traces"
    if save_traces:
        traces_dir.mkdir(parents=True, exist_ok=True)

    patch_pymc_parallel()
    import pymc as pm

    fit_seed = int(seed + point.task_id * 10000 + 1)
    np.random.seed(fit_seed)

    if fit_fn is not None:
        trace = fit_fn(data_syst, comps_syst, xedges_syst, yedges_syst)
    else:
        fik = m2.make_fik(comps_syst, data_syst, xedges_syst, yedges_syst)
        model = m2.m2_model(comps_syst, fik)
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
                random_seed=fit_seed,
            )

    if save_traces and hasattr(trace, "to_netcdf"):
        trace_path = traces_dir / f"trace_bm_threshold_{variant}_{point.point_id}.nc"
        trace.to_netcdf(str(trace_path))

    ndbd_post = trace.posterior["ndbd"].values.flatten()
    ndbd_stats = extract_ndbd_stats(ndbd_post, popt)
    logp_val = extract_logp_mode(trace)

    row = {
        "variant": variant,
        "task_id": point.task_id,
        "point_id": point.point_id,
        "model_key": point.model_key,
        "model_key_syst": model_key_syst,
        "co_shift": point.co_shift,
        "ndbd_shift": point.ndbd_shift,
        "co_cut": str(point.co_cut),
        "ndbd_cutoff": str(point.ndbd_cutoff),
        "ndbd_mode": ndbd_stats["mode"],
        "ndbd_p16": ndbd_stats["p16"],
        "ndbd_p90": ndbd_stats["p90"],
        "max_logp": logp_val,
        "ndbd_mode_bc": ndbd_stats["mode_bc"],
        "ndbd_p16_bc": ndbd_stats["p16_bc"],
        "ndbd_p90_bc": ndbd_stats["p90_bc"],
    }

    for cname in comps_syst.index:
        if cname != "ndbd" and cname in trace.posterior:
            c_post = trace.posterior[cname].values.flatten()
            row[f"{cname}_mode"] = float(get_mode(c_post))

    write_task_h5_atomic(
        task_h5,
        datasets={
            "ndbd_samples": np.array([ndbd_post], dtype=float),
            "max_logp": np.array([logp_val], dtype=float),
        },
        attrs={
            "point_id": point.point_id,
            "variant": variant,
            "model_key_syst": model_key_syst,
            "category": "bm_threshold",
        },
    )

    df_res = pd.DataFrame([row])
    write_task_csv_atomic(df_res, task_csv)

    if verbose:
        print(f"[{point.point_id}] Saved BM threshold fit ({variant}): {task_csv} and {task_h5}")

    return task_csv, task_h5


def collect_syst_bm_threshold(
    study: Any,
    p: Any,
    point_df: pd.DataFrame,
    calib_dict: Dict[str, Any],
    nominal_post: Optional[np.ndarray],
    verbose: bool = True,
) -> Optional[Dict[str, Any]]:
    """Process and collect BM threshold systematic results for a single StudyPoint."""
    _, task_h5 = get_task_paths(study.study_dir, "syst_bm_threshold", p.point_id, with_h5=True)
    if task_h5 is None or not task_h5.exists():
        return None

    with h5py.File(task_h5, "r") as h5f:
        ndbd_samples = h5f["ndbd_samples"][:]

    syst_post_raw = ndbd_samples[0] if ndbd_samples.ndim > 1 else ndbd_samples.flatten()

    nom_post_raw = resolve_nominal_posterior(p, study.cfg, nominal_post)
    if nom_post_raw is None:
        nom_post_raw = syst_post_raw

    popt = None
    if "bc_a" in point_df.columns and "bc_b" in point_df.columns:
        popt = (float(point_df["bc_a"].iloc[0]), float(point_df["bc_b"].iloc[0]))
    else:
        popt = load_bias_popt(calib_dict, study.study_dir, point_id=p.point_id, co_shift=p.co_shift)

    nom_bc = apply_bias_correction(nom_post_raw, popt) if popt is not None else nom_post_raw
    syst_bc = apply_bias_correction(syst_post_raw, popt) if popt is not None else syst_post_raw

    mode_nom, lo_nom, hi_nom = calc_intervals([nom_bc], method="hdi")[0]
    mode_syst, lo_syst, hi_syst = calc_intervals([syst_bc], method="hdi")[0]
    deltas = compute_deltas(nom_bc, syst_bc, return_extensions_only=True)
    delta_left, delta_right = float(deltas[0]), float(deltas[1])
    delta_mode = float(mode_syst - mode_nom)

    scfg = getattr(study.cfg, "syst", None)
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

    syst_h5 = study.study_dir / f"syst_bm_threshold_{p.point_id}.h5"
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

    sqlite_path = study.study_dir / f"{study.name}.sqlite"
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

    plots_dir = study.study_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
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

    return summary_entry


def check_syst_bm_binning(
    point: Any,
    study_dir: Union[str, Path],
    variants: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Check completion status of BM binning systematic fits for this study point."""
    _vars = (
        variants
        if variants is not None
        else ["MinBinSize_5_20", "MinBinSize_25_60", "BinningIntraAlphaBands"]
    )
    return check_task_status(
        study_dir=study_dir,
        prefix="syst_bm_binning",
        point=point,
        target_count=len(_vars),
        target_key="n_samples",
        with_h5=True,
        extra_fields={"category": "bm_binning", "variants": _vars},
    )


def run_syst_bm_binning(
    point: Any,
    cfg: DictConfig,
    study_dir: Union[str, Path],
    variants: Optional[List[str]] = None,
    results_raw_path: Optional[Union[str, Path]] = None,
    seed: int = 42,
    chains: int = 2,
    draws: int = 2500,
    tune: int = 500,
    target_accept: float = 0.9,
    single_core: bool = True,
    bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
    fit_fn: Optional[Any] = None,
    save_traces: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> Tuple[Path, Path]:
    """Execute BM binning systematic fits across variants (MinBinSize, IntraAlpha)."""
    _vars = (
        variants
        if variants is not None
        else ["MinBinSize_5_20", "MinBinSize_25_60", "BinningIntraAlphaBands"]
    )

    task_csv, task_h5 = get_task_paths(study_dir, "syst_bm_binning", point.point_id, with_h5=True)
    assert task_h5 is not None

    if not force:
        status = check_syst_bm_binning(point, study_dir, variants=_vars)
        if status["is_completed"]:
            if verbose:
                print(
                    f"[{point.point_id}] BM binning fits already completed ({status['n_samples_completed']} variants). Skipping."
                )
            return task_csv, task_h5

    raw_path = results_raw_path or config_builders.get_results_raw_path(cfg)
    if not raw_path or not Path(raw_path).exists():
        default_path = Path("/data/cuore/Analysis/BMMC/inputs/ResultsRaw.h5")
        if default_path.exists():
            raw_path = default_path
        else:
            raise FileNotFoundError(
                f"ResultsRaw.h5 not found at '{raw_path}'. Configure infrastructure.data_files.results_raw."
            )

    coords, hm, _, input_dir = storage_fits.resolve_storage_coords(cfg, model_key=point.model_key)
    comp_cfg_raw = [dict(c) for c in config_builders.get_comp_cfg(cfg, model_key=point.model_key)]

    fit_setup = storage_fits.setup_storage_fit(
        cfg=cfg,
        model_key=point.model_key,
        co_cut=point.co_cut,
        ndbd_cutoff=point.ndbd_cutoff,
        co_shift=point.co_shift,
        ndbd_shift=point.ndbd_shift,
        blinded=point.blinded,
        prior_width="5sig",
        mkchain_samples=100,
        chains=chains,
        draws=draws,
        tune=tune,
        target_accept=target_accept,
        dry_run=False,
    )

    data = fit_setup["data"]
    xedges, yedges = fit_setup["edges"]
    nominal_comps = fit_setup["fit_comps"]
    nominal_priors = list(nominal_comps["prior"])

    tpl_uuids, uuid_base, uuid_co, uuid_ndbd = storage_fits.resolve_fit_template_uuids(
        cfg=cfg,
        model_key=point.model_key,
        co_cut=point.co_cut,
        ndbd_cutoff=point.ndbd_cutoff,
        co_shift=point.co_shift,
        ndbd_shift=point.ndbd_shift,
    )

    mctable, signal_df, _, _ = storage_fits.assemble_fit_components(
        input_dir=input_dir,
        jags_fpath=config_builders.get_jags_path(cfg),
        comp_cfg=comp_cfg_raw,
        hm=hm,
        uuid_base=uuid_base,
        uuid_co=uuid_co,
        uuid_ndbd=uuid_ndbd,
    )

    popt = load_bias_popt(
        bias_calibration=bias_calibration,
        study_dir=study_dir,
        point_id=point.point_id,
        co_shift=point.co_shift,
        cfg=cfg,
    )

    traces_dir = Path(study_dir) / "traces"
    if save_traces:
        traces_dir.mkdir(parents=True, exist_ok=True)

    patch_pymc_parallel()
    import pymc as pm

    rows = []
    all_ndbd_samples = []
    all_max_logp = []

    for idx, v in enumerate(_vars):
        v_seed = int(seed + point.task_id * 10000 + idx + 1)
        np.random.seed(v_seed)

        bm_modes_v = bkg_model.read_bm_variant_modes(raw_path, v)
        mt_v = mctable.copy()
        mt_v["norm"] = bm_modes_v.reindex(mt_v.index).fillna(0.0)

        bkg_comps_v = m2.make_bkgcomps_table(comp_cfg_raw, mt_v, smooth=hm.smooth)
        comps_v = pd.concat([bkg_comps_v, signal_df], join="inner")
        comps_v["prior"] = nominal_priors
        comps_v["varname"] = comps_v.index

        if fit_fn is not None:
            trace = fit_fn(data, comps_v, xedges, yedges)
        else:
            fik_v = m2.make_fik(comps_v, data, xedges, yedges)
            model_v = m2.m2_model(comps_v, fik_v)
            with model_v:
                trace = pm.sample(
                    draws=draws,
                    tune=tune,
                    chains=chains,
                    cores=1 if single_core else chains,
                    return_inferencedata=True,
                    target_accept=target_accept,
                    progressbar=False,
                    idata_kwargs={"log_likelihood": True},
                    random_seed=v_seed,
                )

        if save_traces and hasattr(trace, "to_netcdf"):
            trace_path = traces_dir / f"trace_bm_binning_{v}_{point.point_id}.nc"
            trace.to_netcdf(str(trace_path))

        ndbd_post = trace.posterior["ndbd"].values.flatten()
        all_ndbd_samples.append(ndbd_post)

        ndbd_stats = extract_ndbd_stats(ndbd_post, popt)
        logp_val = extract_logp_mode(trace)
        all_max_logp.append(logp_val)

        row = {
            "variant": v,
            "task_id": point.task_id,
            "point_id": point.point_id,
            "model_key": point.model_key,
            "co_shift": point.co_shift,
            "ndbd_shift": point.ndbd_shift,
            "co_cut": str(point.co_cut),
            "ndbd_cutoff": str(point.ndbd_cutoff),
            "ndbd_mode": ndbd_stats["mode"],
            "ndbd_p16": ndbd_stats["p16"],
            "ndbd_p90": ndbd_stats["p90"],
            "max_logp": logp_val,
            "ndbd_mode_bc": ndbd_stats["mode_bc"],
            "ndbd_p16_bc": ndbd_stats["p16_bc"],
            "ndbd_p90_bc": ndbd_stats["p90_bc"],
        }

        for cname in comps_v.index:
            if cname != "ndbd" and cname in trace.posterior:
                c_post = trace.posterior[cname].values.flatten()
                row[f"{cname}_mode"] = float(get_mode(c_post))

        rows.append(row)
        if verbose:
            print(
                f"[{point.point_id}] Completed BM binning fit ({v}): mode={ndbd_stats['mode']:.2f}, mode_bc={ndbd_stats['mode_bc']:.2f}"
            )

    write_task_h5_atomic(
        task_h5,
        datasets={
            "ndbd_samples": np.array(all_ndbd_samples, dtype=float),
            "max_logp": np.array(all_max_logp, dtype=float),
        },
        attrs={
            "point_id": point.point_id,
            "variants": _vars,
            "category": "bm_binning",
        },
    )

    df_res = pd.DataFrame(rows)
    write_task_csv_atomic(df_res, task_csv)

    if verbose:
        print(f"[{point.point_id}] Saved {len(df_res)} BM binning fits: {task_csv} and {task_h5}")

    return task_csv, task_h5


def collect_syst_bm_binning(
    study: Any,
    p: Any,
    point_df: pd.DataFrame,
    calib_dict: Dict[str, Any],
    nominal_post: Optional[np.ndarray],
    verbose: bool = True,
) -> Optional[Dict[str, Any]]:
    """Process and collect BM binning systematic results for a single StudyPoint."""
    _, task_h5 = get_task_paths(study.study_dir, "syst_bm_binning", p.point_id, with_h5=True)
    if task_h5 is None or not task_h5.exists():
        return None

    with h5py.File(task_h5, "r") as h5f:
        ndbd_samples = h5f["ndbd_samples"][:]

    nom_post_raw = resolve_nominal_posterior(p, study.cfg, nominal_post)
    if nom_post_raw is None:
        nom_post_raw = ndbd_samples[0]

    popt = load_bias_popt(calib_dict, study.study_dir, point_id=p.point_id, co_shift=p.co_shift)

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

    syst_h5 = study.study_dir / f"syst_bm_binning_{p.point_id}.h5"
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

    sqlite_path = study.study_dir / f"{study.name}.sqlite"
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
        print(
            f"[{p.point_id}] BM Binning (Envelope): "
            f"Δmode = {combined_d_mode:+.2f}, "
            f"Δright = {combined_d_right:+.2f}, "
            f"Δleft = {combined_d_left:+.2f} "
            f"(mode_nom={mode_nom:.2f})"
        )

    return summary_entry
