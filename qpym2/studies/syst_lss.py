"""Lineshape systematics (LSS: SSCL and dQ) simulation and results collection."""

import gc
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from omegaconf import DictConfig

from qpym2.cfg import config_builders
from qpym2.fit.bias import apply_bias_correction
from qpym2.fit.create_evtrees import create_evtrees
from qpym2.fit.params import model_cfg
from qpym2.fit.syst_lss import (
    build_dq_mctable,
    build_sscl_mctable,
    extract_posterior_stats,
    run_lss_fit,
)
from qpym2.io import bkg_model
from qpym2.io.qls_scaling import load_lss_params, sample_dq_random
from qpym2.storage import fits as storage_fits
from qpym2.storage.study import record_study_point
from qpym2.studies._common import (
    get_task_paths,
    load_bias_popt,
    patch_pymc_parallel,
    resolve_nominal_posterior,
    write_task_csv_atomic,
    write_task_h5_atomic,
)
from qpym2.syst.io import save_syst_result
from qpym2.syst.plots import plot_lss_dq_barplot, plot_lss_sscl_barplot
from qpym2.utils.mathy import calc_intervals, get_mode, signed_quad_sum


def check_syst_lss(
    point: Any,
    study_dir: Union[str, Path],
    category: str = "lss_sscl",
    datasets: Optional[List[int]] = None,
    n_samples: Optional[int] = None,
    dq_mode: str = "pm1",
) -> Dict[str, Any]:
    """Check completion status of Lineshape Systematic (LSS: sscl or dq) fits for this study point."""
    cat_norm = "lss_dq" if category in ("dq", "lss_dq") else "lss_sscl"
    task_csv, task_h5 = get_task_paths(study_dir, f"syst_{cat_norm}", point.point_id, with_h5=True)

    _ds = datasets if datasets is not None else list(range(3601, 3616))
    if cat_norm == "lss_sscl":
        target_n = len(_ds)
    elif dq_mode == "pm1":
        target_n = 2 * len(_ds)
    else:
        target_n = (n_samples or 100) * len(_ds)

    is_done = False
    row_count = 0
    if task_csv.exists() and task_csv.stat().st_size > 0:
        try:
            df = pd.read_csv(task_csv)
            row_count = len(df)
            is_done = row_count >= target_n
        except Exception:
            pass

    return {
        "point_id": point.point_id,
        "task_id": point.task_id,
        "model_key": point.model_key,
        "co_shift": point.co_shift,
        "ndbd_shift": point.ndbd_shift,
        "co_cut": point.co_cut,
        "ndbd_cutoff": point.ndbd_cutoff,
        "category": cat_norm,
        "n_samples_target": target_n,
        "n_samples_completed": row_count,
        "is_completed": is_done,
        "task_csv": str(task_csv),
        "task_h5": str(task_h5) if task_h5 else "",
        "status": "COMPLETED" if is_done else ("IN_PROGRESS" if row_count > 0 else "PENDING"),
    }


def run_syst_lss(
    point: Any,
    cfg: DictConfig,
    study_dir: Union[str, Path],
    category: str = "lss_sscl",
    datasets: Optional[List[int]] = None,
    dq_mode: str = "pm1",
    n_samples: Optional[int] = None,
    lss_data_path: Optional[Union[str, Path]] = None,
    lss_params: Optional[Dict[str, Any]] = None,
    evdata: Optional[Dict[str, Any]] = None,
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
    """Execute Lineshape Systematic (LSS: SSCL or dQ) fits for this study point."""
    cat_norm = "lss_dq" if category in ("dq", "lss_dq") else "lss_sscl"
    task_csv, task_h5 = get_task_paths(study_dir, f"syst_{cat_norm}", point.point_id, with_h5=True)
    assert task_h5 is not None

    _ds_list = datasets if datasets is not None else list(range(3601, 3616))

    if not force:
        status = check_syst_lss(
            point=point,
            study_dir=study_dir,
            category=cat_norm,
            datasets=_ds_list,
            n_samples=n_samples,
            dq_mode=dq_mode,
        )
        if status["is_completed"]:
            if verbose:
                print(
                    f"[{point.point_id}] LSS systematic fits ({cat_norm}) already completed ({status['n_samples_completed']} rows). Skipping."
                )
            return task_csv, task_h5

    coords, hm, param_dir, input_dir = storage_fits.resolve_storage_coords(
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
        chains=chains,
        draws=draws,
        tune=tune,
        target_accept=target_accept,
        dry_run=False,
    )

    data = fit_setup["data"]
    nominal_comps = fit_setup["fit_comps"]
    nominal_priors = list(nominal_comps["prior"])

    if lss_params is None:
        resolved_lss_path = lss_data_path
        if not resolved_lss_path and hasattr(cfg, "syst"):
            resolved_lss_path = getattr(cfg.syst, "lss_data_path", None)
        lss_params = load_lss_params(cfg, n_ds=len(_ds_list), lss_data_path=resolved_lss_path)

    ds_start = _ds_list[0]

    if evdata is None:
        jags_fpath = config_builders.get_jags_path(cfg)
        input_cfg = config_builders.build_input_config(cfg, model_key=point.model_key)
        mcfit_table, _ = bkg_model.read_jagsh5(jags_fpath)
        comp_fnames = dict(zip(mcfit_table.index, mcfit_table.fname))
        comp_fnames["ndbd"] = "ndbd.root"
        evtrees = create_evtrees(input_cfg, comp_fnames)

        _NP_COLS = ["u0", "v0", "ch1", "ch2", "ds"]
        evdata = {}
        for comp, rdf in evtrees.items():
            if rdf is None:
                evdata[comp] = None
            else:
                evdata[comp] = {c: np.array(a) for c, a in rdf.AsNumpy(columns=_NP_COLS).items()}
        del evtrees
        gc.collect()

    mcfg = model_cfg(
        hist_model=hm,
        co_cut=point.co_cut,
        ndbd_cutoff=point.ndbd_cutoff,
        co_shift=point.co_shift,
        ndbd_shift=point.ndbd_shift,
        co_cuts_dict=dict(getattr(cfg.cuts, "co_cuts", {})) if hasattr(cfg, "cuts") else None,
        ndbd_cuts_dict=dict(getattr(cfg.cuts, "ndbd_cutoff", {})) if hasattr(cfg, "cuts") else None,
    )
    cuts = mcfg.get("cuts", {})
    cuts_mc = [f for f in [cuts.get("common"), cuts.get("mc_only")] if f]
    cuts_ndbd = [f for f in [cuts.get("common"), cuts.get("ndbd_only")] if f]

    tpl_uuids, uuid_base, uuid_co, uuid_ndbd = storage_fits.resolve_fit_template_uuids(
        cfg=cfg,
        model_key=point.model_key,
        co_cut=point.co_cut,
        ndbd_cutoff=point.ndbd_cutoff,
        co_shift=point.co_shift,
        ndbd_shift=point.ndbd_shift,
    )

    jags_fpath = config_builders.get_jags_path(cfg)
    mctable, _, _, _ = storage_fits.assemble_fit_components(
        input_dir=input_dir,
        jags_fpath=jags_fpath,
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

    jobs = []
    if cat_norm == "lss_sscl":
        for ds in _ds_list:
            ids = (
                ds - ds_start
                if (ds - ds_start) < len(lss_params["delta_p1_1r"])
                else _ds_list.index(ds)
            )
            ch_sigma0 = lss_params["sigma0_2d"][ids]
            delta_p1 = lss_params["delta_p1_1r"][ids]
            jobs.append((ds, "sscale_1r", {"ch_sigma0": ch_sigma0, "delta_p1": delta_p1}))
    elif dq_mode == "pm1":
        for ds in _ds_list:
            ids = (
                ds - ds_start
                if (ds - ds_start) < len(lss_params["q_pars_arr"])
                else _ds_list.index(ds)
            )
            q_mean = lss_params["q_pars_arr"][ids]
            q_1r = lss_params["q_pars_sample_1r"][ids]
            q_1l = lss_params["q_pars_sample_1l"][ids]
            jobs.append((ds, "dqres_1r", {"q_pars_mean": q_mean, "q_pars_sample": q_1r}))
            jobs.append((ds, "dqres_1l", {"q_pars_mean": q_mean, "q_pars_sample": q_1l}))
    else:
        n_dq = n_samples or 100
        dq_random_samples = sample_dq_random(cfg, n_samples=n_dq, seed=seed)
        for si in range(n_dq):
            sname = f"dq_r{si:04d}"
            qps = dq_random_samples[si]
            for ds in _ds_list:
                ids = (
                    ds - ds_start
                    if (ds - ds_start) < len(lss_params["q_pars_arr"])
                    else _ds_list.index(ds)
                )
                q_mean = lss_params["q_pars_arr"][ids]
                jobs.append((ds, sname, {"q_pars_mean": q_mean, "q_pars_sample": qps[ids]}))

    if verbose:
        print(
            f"[{point.point_id}] Running {len(jobs)} {cat_norm} fits across {len(_ds_list)} datasets..."
        )

    rows = []
    all_ndbd_samples = []

    for idx, (ds, variant, par_dict) in enumerate(jobs):
        fit_seed = int(seed + point.task_id * 10000 + idx + 1)
        rng = np.random.default_rng(fit_seed)

        if cat_norm == "lss_sscl":
            mt_ds, sig_ds = build_sscl_mctable(
                mc_events_dict=evdata,
                mcfit_table=mctable,
                hm=hm,
                target_ds=ds,
                ch_sigma0=par_dict["ch_sigma0"],
                delta_p1=par_dict["delta_p1"],
                cuts_mc=cuts_mc,
                cuts_ndbd=cuts_ndbd,
                rng=rng,
            )
        else:
            mt_ds, sig_ds = build_dq_mctable(
                mc_events_dict=evdata,
                mcfit_table=mctable,
                hm=hm,
                target_ds=ds,
                q_pars_mean=par_dict["q_pars_mean"],
                q_pars_sample=par_dict["q_pars_sample"],
                cuts_mc=cuts_mc,
                cuts_ndbd=cuts_ndbd,
            )

        if fit_fn is not None:
            trace = fit_fn(mt_ds, sig_ds, data, hm, comp_cfg_raw, nominal_priors)
        else:
            trace = run_lss_fit(
                mctable=mt_ds,
                signal_df=sig_ds,
                data=data,
                hm=hm,
                comp_cfg=comp_cfg_raw,
                nominal_priors=nominal_priors,
                chains=chains,
                draws=draws,
            )

        if save_traces and hasattr(trace, "to_netcdf"):
            trace_path = traces_dir / f"trace_{cat_norm}_{variant}_ds{ds}_{point.point_id}.nc"
            trace.to_netcdf(str(trace_path))

        stats = extract_posterior_stats(trace)
        ndbd_post = stats["ndbd_samples"]
        all_ndbd_samples.append(ndbd_post)

        mode_val = float(stats["component_modes"].get("ndbd", get_mode(ndbd_post)))
        p16_val = float(stats["p16"])
        p90_val = float(stats["p90"])

        mode_bc = p16_bc = p90_bc = np.nan
        if popt is not None:
            corrected = apply_bias_correction(ndbd_post, popt)
            mode_bc = float(get_mode(corrected))
            p16_bc = float(np.percentile(corrected, 16))
            p90_bc = float(np.percentile(corrected, 90))

        row = {
            "ds": ds,
            "variant": variant,
            "category": cat_norm,
            "point_id": point.point_id,
            "task_id": point.task_id,
            "model_key": point.model_key,
            "co_shift": point.co_shift,
            "ndbd_shift": point.ndbd_shift,
            "co_cut": str(point.co_cut),
            "ndbd_cutoff": str(point.ndbd_cutoff),
            "ndbd_mode": mode_val,
            "ndbd_p16": p16_val,
            "ndbd_p90": p90_val,
            "ndbd_mode_bc": mode_bc,
            "ndbd_p16_bc": p16_bc,
            "ndbd_p90_bc": p90_bc,
        }

        for cname, cmode in stats["component_modes"].items():
            if cname != "ndbd":
                row[f"{cname}_mode"] = cmode

        rows.append(row)
        if verbose and ((idx + 1) % 5 == 0 or (idx + 1) == len(jobs)):
            print(f"[{point.point_id}] Completed {idx + 1}/{len(jobs)} {cat_norm} fits", flush=True)

    write_task_h5_atomic(
        task_h5,
        datasets={"ndbd_samples": np.array(all_ndbd_samples, dtype=float)},
        attrs={
            "point_id": point.point_id,
            "category": cat_norm,
            "n_fits": len(rows),
        },
    )

    df_res = pd.DataFrame(rows)
    write_task_csv_atomic(df_res, task_csv)

    if verbose:
        print(f"[{point.point_id}] Saved {len(df_res)} {cat_norm} fits: {task_csv} and {task_h5}")

    return task_csv, task_h5


def collect_syst_lss(
    study: Any,
    p: Any,
    point_df: pd.DataFrame,
    category: str,
    calib_dict: Dict[str, Any],
    nominal_post: Optional[np.ndarray],
    verbose: bool = True,
) -> Optional[Dict[str, Any]]:
    """Process and collect LSS (SSCL or dQ) systematic results for a single StudyPoint."""
    nom_post_raw = resolve_nominal_posterior(p, study.cfg, nominal_post)

    popt = load_bias_popt(calib_dict, study.study_dir, point_id=p.point_id, co_shift=p.co_shift)

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
            apply_bias_correction(nom_post_raw, popt) if popt is not None else nom_post_raw
        )
        mode_nom, lo_nom, hi_nom = calc_intervals([nom_post_eval], method="hdi")[0]
    else:
        mode_nom = float(point_df[mode_col].iloc[0])
        lo_nom = float(point_df[p16_col].iloc[0])
        hi_nom = float(point_df[p90_col].iloc[0])

    plots_dir = study.study_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    if category == "lss_sscl":
        point_df["delta_mode"] = point_df[mode_col] - mode_nom
        point_df["delta_lo"] = point_df[p16_col] - lo_nom
        point_df["delta_hi"] = point_df[p90_col] - hi_nom

        deltas = point_df["delta_mode"].values
        sigma_plus, sigma_minus = signed_quad_sum(deltas)
        delta_right = float(sigma_plus)
        delta_left = -float(sigma_minus)
        delta_mode = float(np.mean(deltas))

        syst_h5 = study.study_dir / f"syst_lss_sscl_{p.point_id}.h5"
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

        syst_h5 = study.study_dir / f"syst_lss_dq_{p.point_id}.h5"
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
            f"[{p.point_id}] {category.upper()} Systematic: "
            f"Δright = {delta_right:+.2f}, "
            f"Δleft = {delta_left:+.2f} "
            f"(mode_nom={mode_nom:.2f})"
        )

    return summary_entry
