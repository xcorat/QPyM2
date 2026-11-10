"""Co-60 energy shift (co_shift) marginalization and delta combination."""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from qpym2.fit.bias import apply_bias_correction
from qpym2.plotting.syst_posterior import plot_multi_posterior_overlay, plot_posterior_comparison
from qpym2.storage import fits as storage_fits
from qpym2.storage import paths as storage_paths
from qpym2.storage.study import record_study_point
from qpym2.studies._common import (
    extract_logp_mode,
    load_bias_popt,
)
from qpym2.syst.core import compute_deltas
from qpym2.syst.io import save_syst_result
from qpym2.utils.mathy import calc_intervals, get_mode, sample_weighted_mixture


def fit_co_shift_parabola(
    shifts: Union[List[float], np.ndarray],
    scores: Union[List[float], np.ndarray],
) -> Dict[str, Any]:
    """Fit a 2nd-degree polynomial S(Δ) = aΔ² + bΔ + c to log-likelihood / log-posterior scores.

    Computes the maximum likelihood vertex Δ* = -b / (2a) and the 1-sigma uncertainty
    where ΔS = -0.5: σ_Δ = 1 / sqrt(-2a) if a < 0.

    Parameters
    ----------
    shifts : array-like
        Array of Co-60 shifts in keV.
    scores : array-like
        Array of log-likelihood / log-posterior values corresponding to shifts.

    Returns
    -------
    dict with keys:
        poly_coeffs: [a, b, c]
        vertex_shift: float
        vertex_score: float
        sigma_shift: float
    """
    x = np.asarray(shifts, dtype=float)
    y = np.asarray(scores, dtype=float)

    if len(x) < 3:
        raise ValueError(f"Parabola fit requires at least 3 points, got {len(x)}")

    a, b, c = np.polyfit(x, y, deg=2)
    vertex_x = float(-b / (2.0 * a)) if a != 0 else float(x[np.argmax(y)])
    vertex_y = float(c - (b**2) / (4.0 * a)) if a != 0 else float(np.max(y))

    sigma_x = float(1.0 / np.sqrt(-2.0 * a)) if a < 0 else 0.0

    return {
        "poly_coeffs": [float(a), float(b), float(c)],
        "vertex_shift": vertex_x,
        "vertex_score": vertex_y,
        "sigma_shift": sigma_x,
    }


def calculate_co_shift_weights(
    scores: Union[List[float], np.ndarray],
    draw_counts: Optional[Union[List[int], np.ndarray]] = None,
) -> np.ndarray:
    """Compute softmax weights from log-likelihood / log-posterior scores.

    Optionally scales by (N_avg / n_i) if draw counts differ across fits.
    """
    s = np.asarray(scores, dtype=float)
    weights = np.exp(s - np.max(s))
    weights[weights < 1e-10] = 0.0

    if draw_counts is not None:
        n_draws = np.asarray(draw_counts, dtype=float)
        if len(n_draws) == len(weights) and np.all(n_draws > 0):
            weights = weights * (n_draws.mean() / n_draws)

    total = np.sum(weights)
    if total > 0:
        return weights / total
    return np.full_like(weights, 1.0 / len(weights))


def plot_co_shift_parabola(
    shifts: Union[List[float], np.ndarray],
    scores: Union[List[float], np.ndarray],
    parabola_fit: Dict[str, Any],
    ax: Optional[plt.Axes] = None,
    title: str = "Log-Likelihood Mode vs Co-60 Shift",
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot score points and 2nd-degree polynomial curve with vertex annotation."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.figure

    x = np.asarray(shifts, dtype=float)
    y = np.asarray(scores, dtype=float)
    a, b, c = parabola_fit["poly_coeffs"]
    vx = parabola_fit["vertex_shift"]
    vy = parabola_fit["vertex_score"]
    sig = parabola_fit["sigma_shift"]

    x_dense = np.linspace(min(x.min(), vx - 0.2), max(x.max(), vx + 0.2), 300)
    y_dense = a * x_dense**2 + b * x_dense + c

    ax.scatter(x, y, color="#2b6ea3", s=60, zorder=5, label="Fit scores")
    ax.plot(x_dense, y_dense, color="#1f4e78", lw=1.8, label="Quadratic fit")

    ax.axvline(vx, color="#c97f39", ls="--", lw=1.2, label=f"Vertex: {vx:.3f} keV")
    if sig > 0:
        ax.axvspan(
            vx - sig, vx + sig, color="#e4a66a", alpha=0.2, label=f"±1σ interval ({sig:.3f})"
        )
        ax.axhline(vy - 0.5, color="gray", ls=":", lw=1.0, label="Δlogp = -0.5")

    ax.set_xlabel(r"Co-60 Energy Shift $\Delta_{\mathrm{Co}}$ [keV]", fontsize=12)
    ax.set_ylabel("Score (mode log-likelihood / logp)", fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=9, loc="best", frameon=True)
    ax.grid(True, alpha=0.3)

    return fig, ax


def check_syst_co_shift(study: Any) -> pd.DataFrame:
    """Check fit completion and trace availability across all study points for Co-shift marginalization."""
    records = []
    output_dir = None

    for p in study.points:
        try:
            coords, _, _, _ = storage_fits.resolve_storage_coords(study.cfg, model_key=p.model_key)
            fit_data_root = study.cfg.infrastructure.output_paths.get(
                "fit_data_dir", study.cfg.infrastructure.data_root + "/fit_data"
            )
            param_dir = storage_paths.get_parameter_dir(
                fit_data_root=fit_data_root,
                roi=coords["roi"],
                threshold=coords["threshold"],
                coord=coords["coord"],
                bin_size=coords["bin_size"],
            )
            output_dir = storage_paths.get_output_dir(param_dir)
        except Exception:
            pass

        fit_info = p.get_fit_hash(study.cfg)
        if isinstance(fit_info, tuple):
            _, output_uuid, fit_prefix, _ = fit_info
        else:
            output_uuid = fit_info["output_uuid"]
            fit_prefix = fit_info["fit_prefix"]

        h5_path = Path(str(fit_prefix) + ".h5")
        nc_path = Path(str(fit_prefix) + ".nc")

        is_fitted = h5_path.exists()
        has_trace = nc_path.exists()

        # Fallback check via results.csv if hash aligned across versions
        if not has_trace and output_dir is not None:
            csv_path = study.study_dir / "results.csv"
            if csv_path.exists():
                try:
                    rdf = pd.read_csv(csv_path)
                    match_row = rdf[rdf["co_shift"] == p.co_shift]
                    if not match_row.empty and "output_uuid" in match_row.columns:
                        old_uuid = str(match_row["output_uuid"].iloc[0])
                        old_nc = output_dir / f"fit_{old_uuid}.nc"
                        if old_nc.exists():
                            has_trace = True
                            is_fitted = True
                            output_uuid = old_uuid
                except Exception:
                    pass

        status = "READY" if has_trace else ("FITTED_NO_TRACE" if is_fitted else "PENDING_FIT")
        records.append(
            {
                "point_id": p.point_id,
                "task_id": p.task_id,
                "model_key": p.model_key,
                "co_shift": p.co_shift,
                "output_uuid": output_uuid,
                "is_fitted": is_fitted,
                "has_trace": has_trace,
                "status": status,
                "is_completed": has_trace,
            }
        )

    return pd.DataFrame(records)


def run_syst_co_shift(
    study: Any,
    chains: Optional[int] = None,
    draws: Optional[int] = None,
    tune: Optional[int] = None,
    target_accept: Optional[float] = None,
    force: bool = False,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Execute any pending MCMC fits across the study points, ensuring traces are saved."""
    results = {}
    for p in study.points:
        fit_info = p.get_fit_hash(study.cfg)
        if isinstance(fit_info, tuple):
            _, output_uuid, fit_prefix, _ = fit_info
        else:
            output_uuid = fit_info["output_uuid"]
            fit_prefix = fit_info["fit_prefix"]
        nc_path = Path(str(fit_prefix) + ".nc")

        if not force and nc_path.exists():
            if verbose:
                print(f"[{p.point_id}] Fit with trace already exists ({output_uuid}). Skipping.")
            results[p.point_id] = {"status": "SKIPPED", "output_uuid": output_uuid}
            continue

        if verbose:
            print(f"[{p.point_id}] Running MCMC fit with trace for co_shift={p.co_shift} keV...")

        res = p.run_storage_fit(
            study.cfg,
            chains=chains or 8,
            draws=draws or 20000,
            tune=tune or 1000,
            target_accept=target_accept or 0.9,
            save_trace=True,
            force=force,
            verbose=verbose,
        )
        results[p.point_id] = res

    return results


def collect_syst_co_shift(
    study: Any,
    bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
    reference_shift: Optional[float] = None,
    lp_var: str = "pik",
    estimator: str = "mode",
    n_mixture_draws: int = 160000,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Collect Co-shift posteriors, calculate likelihood parabola and weights, build marginal mixture,

    and extract signed HDI edge deltas against the nominal reference shift.
    """
    scfg = getattr(study.cfg, "syst", None)
    _lp_var = getattr(scfg, "lp_var", lp_var) if scfg else lp_var
    _estimator = getattr(scfg, "estimator", estimator) if scfg else estimator

    _ref_shift = reference_shift
    if _ref_shift is None and scfg and hasattr(scfg, "reference_co_shift"):
        _ref_shift = float(scfg.reference_co_shift)
    if (
        _ref_shift is None
        and hasattr(study.cfg, "reference")
        and hasattr(study.cfg.reference, "co_shift")
    ):
        _ref_shift = float(study.cfg.reference.co_shift)
    if _ref_shift is None:
        _ref_shift = 1.0

    plots_dir = study.study_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    sqlite_path = study.study_dir / f"{study.name}.sqlite"

    # 1. Resolve traces and output UUIDs for each point
    shift_points = sorted(study.points, key=lambda p: float(p.co_shift))
    shifts: List[float] = []
    raw_posts: List[np.ndarray] = []
    scores: List[float] = []
    draw_counts: List[int] = []

    # Map results.csv fallback UUIDs if present
    results_uuid_map: Dict[float, str] = {}
    csv_path = study.study_dir / "results.csv"
    if csv_path.exists():
        try:
            rdf = pd.read_csv(csv_path)
            for _, r in rdf.iterrows():
                if "co_shift" in r and "output_uuid" in r:
                    results_uuid_map[float(r["co_shift"])] = str(r["output_uuid"])
        except Exception:
            pass

    for p in shift_points:
        fit_info = p.get_fit_hash(study.cfg)
        if isinstance(fit_info, tuple):
            _, _, fit_prefix, _ = fit_info
        else:
            fit_prefix = fit_info["fit_prefix"]
        nc_path = Path(str(fit_prefix) + ".nc")

        if not nc_path.exists() and p.co_shift in results_uuid_map:
            alt_uuid = results_uuid_map[p.co_shift]
            alt_nc = fit_prefix.parent / f"fit_{alt_uuid}.nc"
            if alt_nc.exists():
                nc_path = alt_nc

        if not nc_path.exists():
            raise FileNotFoundError(
                f"Missing NetCDF trace for point {p.point_id} (co_shift={p.co_shift}). "
                f"Expected trace at {nc_path}."
            )

        idata = az.from_netcdf(str(nc_path))
        ndbd_arr = idata.posterior["ndbd"].values.flatten()

        # Score extraction from log_likelihood or sample_stats
        score_val = 0.0
        if hasattr(idata, "log_likelihood") and _lp_var in idata.log_likelihood:
            lp_samples = idata.log_likelihood[_lp_var].values.flatten()
            score_val = float(get_mode(lp_samples) if _estimator == "mode" else np.mean(lp_samples))
        elif hasattr(idata, "sample_stats") and "lp" in idata.sample_stats:
            lp_samples = idata.sample_stats["lp"].values.flatten()
            score_val = float(get_mode(lp_samples) if _estimator == "mode" else np.mean(lp_samples))
        else:
            # Fallback to mode of ndbd or 0.0
            score_val = float(extract_logp_mode(idata) or 0.0)

        shifts.append(float(p.co_shift))
        raw_posts.append(ndbd_arr)
        scores.append(score_val)
        draw_counts.append(len(ndbd_arr))

    # 2. Fit 2nd-degree parabola and calculate weights
    parabola = fit_co_shift_parabola(shifts, scores)
    weights = calculate_co_shift_weights(scores, draw_counts=draw_counts)

    # 3. Bias correction for each shift point
    bc_posts: List[np.ndarray] = []
    per_shift_summaries: Dict[str, Any] = {}

    for i, (p, s_val, w_val, raw_arr) in enumerate(zip(shift_points, shifts, weights, raw_posts)):
        popt = load_bias_popt(
            bias_calibration=bias_calibration,
            study_dir=study.study_dir,
            point_id=p.point_id,
            co_shift=s_val,
            cfg=study.cfg,
        )
        bc_arr = apply_bias_correction(raw_arr, popt) if popt is not None else raw_arr
        bc_posts.append(bc_arr)

        m_raw, lo_raw, hi_raw = calc_intervals([raw_arr], method="hdi")[0]
        m_bc, lo_bc, hi_bc = calc_intervals([bc_arr], method="hdi")[0]

        p_info = {
            "point_id": p.point_id,
            "co_shift": s_val,
            "score": scores[i],
            "weight": float(w_val),
            "mode_raw": float(m_raw),
            "hdi_raw": [float(lo_raw), float(hi_raw)],
            "mode_bc": float(m_bc),
            "hdi_bc": [float(lo_bc), float(hi_bc)],
            "bc_a": float(popt[0]) if popt else 1.0,
            "bc_b": float(popt[1]) if popt else 0.0,
        }
        per_shift_summaries[f"shift_{s_val:g}"] = p_info

    # 4. Form weighted mixture distribution
    rng = np.random.default_rng(42)
    marg_mixture = sample_weighted_mixture(
        bc_posts, weights=weights, n_draws=n_mixture_draws, rng=rng
    )

    # 5. Resolve nominal reference posterior
    ref_idx = 0
    min_dist = float("inf")
    for idx, s in enumerate(shifts):
        if abs(s - _ref_shift) < min_dist:
            min_dist = abs(s - _ref_shift)
            ref_idx = idx

    nom_bc = bc_posts[ref_idx]
    actual_ref_shift = shifts[ref_idx]

    nom_mode, nom_lo, nom_hi = calc_intervals([nom_bc], method="hdi")[0]
    marg_mode, marg_lo, marg_hi = calc_intervals([marg_mixture], method="hdi")[0]

    deltas = compute_deltas(nom_bc, marg_mixture, return_extensions_only=True)
    delta_left = float(deltas[0])  # <= 0
    delta_right = float(deltas[1])  # >= 0
    delta_mode = float(marg_mode - nom_mode)

    # 6. Assemble summary payload
    summary_entry = {
        "category": "co_shift",
        "reference_co_shift": actual_ref_shift,
        "vertex_co_shift": float(parabola["vertex_shift"]),
        "vertex_score": float(parabola["vertex_score"]),
        "sigma_co_shift": float(parabola["sigma_shift"]),
        "poly_coeffs": [float(c) for c in parabola["poly_coeffs"]],
        "mode_nom": float(nom_mode),
        "mode_syst": float(marg_mode),
        "delta_mode": float(delta_mode),
        "hdi_nom": [float(nom_lo), float(nom_hi)],
        "hdi_syst": [float(marg_lo), float(marg_hi)],
        "delta_right": float(delta_right),
        "delta_left": float(delta_left),
        "shifts": shifts,
        "weights": {f"{s:g}": float(w) for s, w in zip(shifts, weights)},
        "scores": {f"{s:g}": float(sc) for s, sc in zip(shifts, scores)},
        "per_shift": per_shift_summaries,
    }

    # Format compliant for both direct indexing and marginalize.py first_entry extraction
    summaries = {
        "nominal": summary_entry,
        "co_shift_marginal": summary_entry,
    }

    summary_file = study.study_dir / "syst_summary.json"
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)

    # Save Format v2 HDF5 result
    syst_h5 = study.study_dir / "syst_co_shift_marginal.h5"
    model_pars = {
        "co_cut": getattr(shift_points[ref_idx], "co_cut", None),
        "ndbd_cutoff": getattr(shift_points[ref_idx], "ndbd_cutoff", None),
        "co_shift": actual_ref_shift,
    }
    save_syst_result(
        syst_h5,
        name="co_shift",
        bc_post=marg_mixture,
        raw_post=raw_posts[ref_idx],
        model_pars=model_pars,
        overwrite=True,
    )

    # Update SQLite records
    record_study_point(
        sqlite_path,
        {
            "point_id": "co_shift_marginal",
            "study_name": study.name,
            "task_id": 999,
            "status": "completed",
            "summary_json": summary_entry,
        },
    )

    # 7. Diagnostic Plots
    # A. Parabola fit plot
    try:
        p_fig, p_ax = plot_co_shift_parabola(
            shifts=shifts,
            scores=scores,
            parabola_fit=parabola,
            title=f"Co-60 Log-Likelihood Parabola ({study.name})",
        )
        p_fig.savefig(plots_dir / "co_shift_logp_parabola.png", dpi=150, bbox_inches="tight")
        plt.close(p_fig)
    except Exception as e:
        if verbose:
            print(f"Note: Could not render parabola plot ({e})")

    # B. Multi-posterior overlay
    try:
        o_fig, o_ax = plt.subplots(figsize=(8.5, 4.8))
        labels = [f"cs={s:.1f} keV (w={w:.3f})" for s, w in zip(shifts, weights)]
        plot_multi_posterior_overlay(
            samples_list=bc_posts,
            labels=labels,
            ax=o_ax,
            bins=80,
            highlight_idx=ref_idx,
            xlabel=r"Bias-corrected $n_{0\nu\beta\beta}$ [counts]",
            title="Bias-Corrected Posteriors across Co-60 Shifts",
        )
        o_ax.legend(fontsize=8, loc="upper right")
        o_fig.savefig(plots_dir / "co_shift_posterior_overlay.png", dpi=150, bbox_inches="tight")
        plt.close(o_fig)
    except Exception as e:
        if verbose:
            print(f"Note: Could not render overlay plot ({e})")

    # C. Marginal mixture vs nominal comparison
    try:
        m_fig, m_ax = plt.subplots(figsize=(8.5, 4.8))
        plot_posterior_comparison(
            syst_samples=marg_mixture,
            nom_samples=nom_bc,
            ax=m_ax,
            bins=80,
            syst_label="Co-shift Marginal",
            nom_label=f"Nominal (cs={actual_ref_shift:.1f} keV)",
            xlabel=r"$n_{0\nu\beta\beta}$ [counts]",
            title="Co-60 Energy Shift Systematic: Marginal vs Nominal",
        )
        m_ax.legend(fontsize=8, loc="upper right")
        m_fig.savefig(plots_dir / "co_shift_marginal_comparison.png", dpi=150, bbox_inches="tight")
        plt.close(m_fig)
    except Exception as e:
        if verbose:
            print(f"Note: Could not render comparison plot ({e})")

    if verbose:
        print(
            f"[{study.name}] Co-Shift Marginalization Completed:\n"
            f"  Reference shift:    {actual_ref_shift:.1f} keV\n"
            f"  Optimal vertex:     {parabola['vertex_shift']:.3f} ± {parabola['sigma_shift']:.3f} keV\n"
            f"  Weights:            {', '.join([f'{s:g}: {w:.3f}' for s, w in zip(shifts, weights)])}\n"
            f"  Nominal mode:       {nom_mode:.2f}  HDI: [{nom_lo:.2f}, {nom_hi:.2f}]\n"
            f"  Marginalized mode:  {marg_mode:.2f}  HDI: [{marg_lo:.2f}, {marg_hi:.2f}]\n"
            f"  Deltas:             Δmode = {delta_mode:+.2f}, Δright = {delta_right:+.2f}, Δleft = {delta_left:+.2f}\n"
            f"  Saved summary:      {summary_file}"
        )

    return {"summary": summaries, "mixture_samples": marg_mixture}
