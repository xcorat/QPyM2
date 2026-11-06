"""Background model Monte Carlo (bmmc / mc_stats) systematic simulation and collection."""

from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

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
from qpym2.utils.mathy import calc_intervals, get_mode, sample_weighted_mixture


def check_syst_mc(
    point: Any,
    study_dir: Union[str, Path],
    n_samples: int = 100,
) -> Dict[str, Any]:
    """Check completion status of MC-stats (bmmc) systematic fits for this study point."""
    return check_task_status(
        study_dir=study_dir,
        prefix="syst_mc",
        point=point,
        target_count=n_samples,
        target_key="n_samples",
        with_h5=True,
    )


def run_syst_mc(
    point: Any,
    cfg: DictConfig,
    study_dir: Union[str, Path],
    n_samples: int = 100,
    seed: int = 42,
    chains: int = 2,
    draws: int = 2500,
    tune: int = 500,
    target_accept: float = 0.9,
    stat_only: bool = False,
    single_core: bool = True,
    bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
    fit_fn: Optional[Any] = None,
    save_traces: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> Tuple[Path, Path]:
    """Execute a batch of MC-stats (bmmc) systematic fits for this study point."""
    task_csv, task_h5 = get_task_paths(study_dir, "syst_mc", point.point_id, with_h5=True)
    assert task_h5 is not None

    if not force:
        status = check_syst_mc(point, study_dir, n_samples=n_samples)
        if status["is_completed"]:
            if verbose:
                print(
                    f"[{point.point_id}] MC-stats systematic fits already completed ({status['n_samples_completed']} rows). Skipping."
                )
            return task_csv, task_h5

    coords, hm, param_dir, input_dir = storage_fits.resolve_storage_coords(
        cfg, model_key=point.model_key
    )
    jags_fpath = config_builders.get_jags_path(cfg)
    corrmat_fpath = config_builders.get_corrmat_path(cfg) if not stat_only else None
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
        jags_fpath=jags_fpath,
        comp_cfg=comp_cfg_raw,
        hm=hm,
        uuid_base=uuid_base,
        uuid_co=uuid_co,
        uuid_ndbd=uuid_ndbd,
    )

    if verbose:
        mode_str = "stat-only" if stat_only else "covariance-perturbed"
        print(
            f"[{point.point_id}] Generating {n_samples} {mode_str} normalisation vectors (seed={seed})..."
        )

    mc_samples = bkg_model.get_mcsyst_samples(
        chain=jags_fpath,
        corrmat_fpath=corrmat_fpath,
        n=n_samples,
        seed=seed,
        stat_only=stat_only,
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

    for i in range(n_samples):
        sample_seed = int(seed + point.task_id * 10000 + i + 1)
        np.random.seed(sample_seed)

        mt_i = mctable.copy()
        mt_i["norm"] = pd.Series(mc_samples.iloc[i].to_dict()).reindex(mt_i.index)

        bkg_comps_i = m2.make_bkgcomps_table(comp_cfg_raw, mt_i, smooth=hm.smooth)
        comps_i = pd.concat([bkg_comps_i, signal_df], join="inner")
        comps_i["prior"] = nominal_priors
        comps_i["varname"] = comps_i.index

        if fit_fn is not None:
            trace = fit_fn(data, comps_i, xedges, yedges)
        else:
            fik_i = m2.make_fik(comps_i, data, xedges, yedges)
            model_i = m2.m2_model(comps_i, fik_i)
            with model_i:
                trace = pm.sample(
                    draws=draws,
                    tune=tune,
                    chains=chains,
                    cores=1 if single_core else chains,
                    return_inferencedata=True,
                    target_accept=target_accept,
                    progressbar=False,
                    idata_kwargs={"log_likelihood": True},
                    random_seed=sample_seed,
                )

        if save_traces and hasattr(trace, "to_netcdf"):
            trace_path = traces_dir / f"trace_mc_{point.point_id}_s{i:04d}.nc"
            trace.to_netcdf(str(trace_path))

        ndbd_post = trace.posterior["ndbd"].values.flatten()
        all_ndbd_samples.append(ndbd_post)

        ndbd_stats = extract_ndbd_stats(ndbd_post, popt)
        logp_val = extract_logp_mode(trace)
        all_max_logp.append(logp_val)

        row = {
            "sample_idx": i,
            "task_id": point.task_id,
            "point_id": point.point_id,
            "model_key": point.model_key,
            "co_shift": point.co_shift,
            "ndbd_shift": point.ndbd_shift,
            "co_cut": str(point.co_cut),
            "ndbd_cutoff": str(point.ndbd_cutoff),
            "stat_only": stat_only,
            "ndbd_mode": ndbd_stats["mode"],
            "ndbd_p16": ndbd_stats["p16"],
            "ndbd_p90": ndbd_stats["p90"],
            "max_logp": logp_val,
            "ndbd_mode_bc": ndbd_stats["mode_bc"],
            "ndbd_p16_bc": ndbd_stats["p16_bc"],
            "ndbd_p90_bc": ndbd_stats["p90_bc"],
        }

        for cname in comps_i.index:
            if cname != "ndbd" and cname in trace.posterior:
                c_post = trace.posterior[cname].values.flatten()
                row[f"{cname}_mode"] = float(get_mode(c_post))

        rows.append(row)
        if verbose and ((i + 1) % 10 == 0 or (i + 1) == n_samples):
            print(
                f"[{point.point_id}] Completed {i + 1}/{n_samples} MC-stats systematic fits",
                flush=True,
            )

    write_task_h5_atomic(
        task_h5,
        datasets={
            "ndbd_samples": np.array(all_ndbd_samples, dtype=float),
            "max_logp": np.array(all_max_logp, dtype=float),
        },
        attrs={
            "point_id": point.point_id,
            "n_samples": n_samples,
            "stat_only": stat_only,
        },
    )

    df_res = pd.DataFrame(rows)
    write_task_csv_atomic(df_res, task_csv)

    if verbose:
        print(f"[{point.point_id}] Saved {len(df_res)} MC-stats fits: {task_csv} and {task_h5}")

    return task_csv, task_h5


def collect_syst_mc(
    study: Any,
    p: Any,
    point_df: pd.DataFrame,
    calib_dict: Dict[str, Any],
    nominal_post: Optional[np.ndarray],
    verbose: bool = True,
) -> Optional[Dict[str, Any]]:
    """Process and collect MC-stats systematic results for a single StudyPoint."""
    _, task_h5 = get_task_paths(study.study_dir, "syst_mc", p.point_id, with_h5=True)
    if task_h5 is None or not task_h5.exists():
        return None

    with h5py.File(task_h5, "r") as h5f:
        ndbd_samples = h5f["ndbd_samples"][:]
        max_logp = h5f["max_logp"][:]

    if len(ndbd_samples) == 0:
        return None

    scores = max_logp - np.max(max_logp)
    weights = np.exp(scores)
    weights[weights < 1e-10] = 0.0
    if weights.sum() > 0:
        weights = weights / weights.sum()
    else:
        weights = np.full(len(max_logp), 1.0 / len(max_logp))

    rng = np.random.default_rng(42)
    combined_samples = sample_weighted_mixture(ndbd_samples, weights, n_draws=50000, rng=rng)

    nom_post_raw = resolve_nominal_posterior(p, study.cfg, nominal_post)
    if nom_post_raw is None:
        nom_post_raw = ndbd_samples[0]

    popt = load_bias_popt(calib_dict, study.study_dir, point_id=p.point_id, co_shift=p.co_shift)

    if popt is not None:
        nom_bc = apply_bias_correction(nom_post_raw, popt)
        syst_bc = apply_bias_correction(combined_samples, popt)
    else:
        nom_bc = nom_post_raw
        syst_bc = combined_samples

    mode_nom, lo_nom, hi_nom = calc_intervals([nom_bc], method="hdi")[0]
    mode_syst, lo_syst, hi_syst = calc_intervals([syst_bc], method="hdi")[0]
    deltas = compute_deltas(nom_bc, syst_bc, return_extensions_only=True)
    delta_left = float(deltas[0])
    delta_right = float(deltas[1])
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

    syst_h5 = study.study_dir / f"syst_bmmc_{p.point_id}.h5"
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

    return summary_entry
