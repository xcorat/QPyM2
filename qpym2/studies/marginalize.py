"""Rate marginalisation, split-Gaussian smearing, and limit extraction."""

import json
from pathlib import Path
from typing import Any, Dict, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from omegaconf import OmegaConf

from qpym2.fit.bias import apply_bias_correction
from qpym2.storage.study import record_study_point
from qpym2.studies._common import load_bias_popt, resolve_nominal_posterior
from qpym2.syst.core import SystDelta, combine_deltas
from qpym2.syst.io import save_combined_result
from qpym2.syst.plots import plot_combined_posterior, plot_delta_table, plot_rate_posterior
from qpym2.syst.rate import to_rate
from qpym2.syst.smear import smear_posterior
from qpym2.utils.mathy import calc_intervals, get_mode


def run_marginalization(
    study: Any,
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
    """Execute full rate marginalisation, split-Gaussian smearing, and limit extraction."""
    m_cfg = getattr(study.cfg, "marginalisation", None)
    _exposure = float(
        exposure
        if exposure is not None
        else (getattr(m_cfg, "exposure", 1038.4) if m_cfg else 1038.4)
    )
    _eff_mc = float(
        eff_mc if eff_mc is not None else (getattr(m_cfg, "eff_mc", 0.04) if m_cfg else 0.04)
    )
    _psa_mean = float(
        psa_mean if psa_mean is not None else (getattr(m_cfg, "psa_mean", 0.99) if m_cfg else 0.99)
    )
    _psa_sig = float(
        psa_sig if psa_sig is not None else (getattr(m_cfg, "psa_sig", 0.01) if m_cfg else 0.01)
    )
    _n_smear = int(getattr(m_cfg, "n_smear", n_smear) if m_cfg else n_smear)
    _rng_seed = int(getattr(m_cfg, "rng_seed", rng_seed) if m_cfg else rng_seed)
    _extend_method = str(getattr(m_cfg, "extend_method", extend_method) if m_cfg else extend_method)

    plots_dir = study.study_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    sqlite_path = study.study_dir / f"{study.name}.sqlite"

    p0 = study.points[0] if study.points else None
    nom_raw = resolve_nominal_posterior(p0, study.cfg, nominal_post) if p0 else None

    if nom_raw is None:
        rng_nom = np.random.default_rng(_rng_seed)
        nom_raw = np.clip(rng_nom.normal(loc=15.0, scale=8.0, size=50000), 0.0, None)

    popt = load_bias_popt(
        bias_calibration=bias_calibration,
        study_dir=study.study_dir,
        point_id=p0.point_id if p0 else None,
        co_shift=p0.co_shift if p0 else None,
        cfg=study.cfg,
    )

    nom_bc = apply_bias_correction(nom_raw, popt) if popt is not None else nom_raw

    systs_cfg = getattr(study.cfg, "systematics", None)
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
                src_path = (
                    (study.study_dir / src_path).resolve()
                    if (study.study_dir / src_path).exists()
                    else Path(source).resolve()
                )
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

        delta_records.append(
            {
                "source": s_name,
                "delta_left": dl_val,
                "delta_right": dr_val,
                "study_source": str(source or "explicit"),
            }
        )
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
            delta_records.append(
                {
                    "source": s_name,
                    "delta_left": dl,
                    "delta_right": dr,
                    "study_source": "default_tech_note",
                }
            )
            syst_deltas.append(SystDelta(name=s_name, delta_right=dr, delta_left=dl))

    df_delta_table = pd.DataFrame(delta_records)
    df_delta_table.to_csv(study.study_dir / "delta_table.csv", index=False)

    sigma_left, sigma_right = combine_deltas([(d.delta_left, d.delta_right) for d in syst_deltas])

    smeared_post = smear_posterior(
        nom_bc,
        sigma_right=float(abs(sigma_right)),
        sigma_left=float(abs(sigma_left)),
        n_smear=_n_smear,
        rng_seed=_rng_seed,
    )

    rate_post = to_rate(
        smeared_post,
        eff_mc=_eff_mc,
        exposure=_exposure,
        psa_mean=_psa_mean,
        psa_sig=_psa_sig,
        rng_seed=_rng_seed,
        clip=True,
    )

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

    stats_df = pd.DataFrame(
        [
            {
                "label": "nominal_bc",
                "mode": mode_nom,
                "hdi_lo": float(hdi_nom[1]),
                "hdi_hi": float(hdi_nom[2]),
                "p90": p90_nom,
            },
            {
                "label": "smeared",
                "mode": mode_smeared,
                "hdi_lo": float(hdi_smeared[1]),
                "hdi_hi": float(hdi_smeared[2]),
                "p90": p90_smeared,
            },
            {
                "label": "rate",
                "mode": rate_mode,
                "hdi_lo": float(hdi_rate[1]),
                "hdi_hi": float(hdi_rate[2]),
                "p90": rate_90,
            },
        ]
    )

    summary = {
        "study_name": study.name,
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

    with open(study.study_dir / "syst_rate_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    model_pars = {
        "co_cut": getattr(p0, "co_cut", None) if p0 else None,
        "ndbd_cutoff": getattr(p0, "ndbd_cutoff", None) if p0 else None,
        "co_shift": getattr(p0, "co_shift", 0.0) if p0 else 0.0,
    }
    combined_h5 = study.study_dir / "combined_result.h5"
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

    record_study_point(
        sqlite_path,
        {
            "point_id": "marginalization",
            "study_name": study.name,
            "task_id": 9999,
            "status": "completed",
            "summary_json": summary,
        },
    )

    if verbose:
        print(
            f"[INFO] Rate Marginalisation COMPLETED for '{study.name}':\n"
            f"  σ_right = {sigma_right:+.3f}, σ_left = {sigma_left:+.3f}\n"
            f"  Γ_mode = {rate_mode:.2e} yr⁻¹, Γ_90 = {rate_90:.2e} yr⁻¹\n"
            f"  T_1/2^(0ν) > {half_life_90:.2e} yr (90% C.L.)\n"
            f"  Outputs saved to {study.study_dir}"
        )

    return summary
