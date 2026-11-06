"""Template binning (bin_size) systematics simulation and collection."""

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import h5py
import numpy as np
import pandas as pd
from omegaconf import DictConfig

from qpym2.cfg import config_builders
from qpym2.fit.bias import apply_bias_correction
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
from qpym2.utils import TempConfig
from qpym2.utils.hists import get_empty_hist
from qpym2.utils.mathy import calc_intervals, get_mode


def check_syst_template_binning(
    point: Any,
    study_dir: Union[str, Path],
    bin_sizes: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Check completion status of template binning (bin_size) systematic fits."""
    _sizes = bin_sizes if bin_sizes is not None else [2.0]
    return check_task_status(
        study_dir=study_dir,
        prefix="syst_bin_size",
        point=point,
        target_count=len(_sizes),
        target_key="n_samples",
        with_h5=True,
        extra_fields={"category": "bin_size", "bin_sizes": _sizes},
    )


def run_syst_template_binning(
    point: Any,
    cfg: DictConfig,
    study_dir: Union[str, Path],
    bin_sizes: Optional[List[float]] = None,
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
    """Execute template likelihood binning fits for coarser bin sizes (e.g. 2.0 or 5.0 keV)."""
    _sizes = [float(b) for b in (bin_sizes if bin_sizes is not None else [2.0])]

    task_csv, task_h5 = get_task_paths(study_dir, "syst_bin_size", point.point_id, with_h5=True)
    assert task_h5 is not None

    if not force:
        status = check_syst_template_binning(point, study_dir, bin_sizes=_sizes)
        if status["is_completed"]:
            if verbose:
                print(
                    f"[{point.point_id}] Template binning fits already completed ({status['n_samples_completed']} sizes). Skipping."
                )
            return task_csv, task_h5

    coords, hm_base, _, input_dir = storage_fits.resolve_storage_coords(
        cfg, model_key=point.model_key
    )
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
        hm=hm_base,
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

    for idx, bs in enumerate(_sizes):
        bs_seed = int(seed + point.task_id * 10000 + idx + 1)
        np.random.seed(bs_seed)

        hm_bs = TempConfig(**{k: getattr(hm_base, k) for k in vars(hm_base)})
        hm_bs.bin_size = bs
        hm_bs.binsize = bs / math.sqrt(2)
        _, xedges_bs, yedges_bs = get_empty_hist(hm_bs)

        bkg_comps_bs = m2.make_bkgcomps_table(comp_cfg_raw, mctable, smooth=hm_base.smooth)
        comps_bs = pd.concat([bkg_comps_bs, signal_df], join="inner")
        comps_bs["prior"] = nominal_priors
        comps_bs["varname"] = comps_bs.index

        if fit_fn is not None:
            trace = fit_fn(data, comps_bs, xedges_bs, yedges_bs)
        else:
            fik_bs = m2.make_fik(comps_bs, data, xedges_bs, yedges_bs)
            model_bs = m2.m2_model(comps_bs, fik_bs)
            with model_bs:
                trace = pm.sample(
                    draws=draws,
                    tune=tune,
                    chains=chains,
                    cores=1 if single_core else chains,
                    return_inferencedata=True,
                    target_accept=target_accept,
                    progressbar=False,
                    idata_kwargs={"log_likelihood": True},
                    random_seed=bs_seed,
                )

        if save_traces and hasattr(trace, "to_netcdf"):
            trace_path = traces_dir / f"trace_bin_size_{bs:.1f}_{point.point_id}.nc"
            trace.to_netcdf(str(trace_path))

        ndbd_post = trace.posterior["ndbd"].values.flatten()
        all_ndbd_samples.append(ndbd_post)

        ndbd_stats = extract_ndbd_stats(ndbd_post, popt)
        logp_val = extract_logp_mode(trace)
        all_max_logp.append(logp_val)

        row = {
            "bin_size": bs,
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

        for cname in comps_bs.index:
            if cname != "ndbd" and cname in trace.posterior:
                c_post = trace.posterior[cname].values.flatten()
                row[f"{cname}_mode"] = float(get_mode(c_post))

        rows.append(row)
        if verbose:
            print(
                f"[{point.point_id}] Completed template binning fit ({bs} keV): mode={ndbd_stats['mode']:.2f}, mode_bc={ndbd_stats['mode_bc']:.2f}"
            )

    write_task_h5_atomic(
        task_h5,
        datasets={
            "ndbd_samples": np.array(all_ndbd_samples, dtype=float),
            "max_logp": np.array(all_max_logp, dtype=float),
        },
        attrs={
            "point_id": point.point_id,
            "bin_sizes": _sizes,
            "category": "bin_size",
        },
    )

    df_res = pd.DataFrame(rows)
    write_task_csv_atomic(df_res, task_csv)

    if verbose:
        print(
            f"[{point.point_id}] Saved {len(df_res)} template binning fits: {task_csv} and {task_h5}"
        )

    return task_csv, task_h5


def collect_syst_template_binning(
    study: Any,
    p: Any,
    point_df: pd.DataFrame,
    calib_dict: Dict[str, Any],
    nominal_post: Optional[np.ndarray],
    verbose: bool = True,
) -> Optional[Dict[str, Any]]:
    """Process and collect template binning systematic results for a single StudyPoint."""
    _, task_h5 = get_task_paths(study.study_dir, "syst_bin_size", p.point_id, with_h5=True)
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

    syst_h5 = study.study_dir / f"syst_bin_size_{p.point_id}.h5"
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
            f"[{p.point_id}] Template Binning (Envelope): "
            f"Δmode = {combined_d_mode:+.2f}, "
            f"Δright = {combined_d_right:+.2f}, "
            f"Δleft = {combined_d_left:+.2f} "
            f"(mode_nom={mode_nom:.2f})"
        )

    return summary_entry
