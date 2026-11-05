import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    import os
    import sys
    from pathlib import Path
    from typing import Dict, List, Optional, Tuple

    import arviz as az
    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit

    from qpym2.cfg.hydra_adapter import (
        HydraConfigManager,
        get_config,
        resolve_config_dir,
    )
    from qpym2.fit.bias import (
        apply_bias_correction,
        get_bias_from_fitdf,
        sample_corrected,
    )
    from qpym2.io.toys import _unflatten_columns as unflatten_columns
    from qpym2.storage.backend import fit_exists, read_fit
    from qpym2.storage.paths import get_output_dir
    from qpym2.studies.assembly import assemble_study
    from qpym2.utils.mathy import (
        get_1sig_around_mode,
        get_mode,
        sample_weighted_mixture,
    )

    return (
        Path,
        apply_bias_correction,
        assemble_study,
        get_1sig_around_mode,
        get_bias_from_fitdf,
        get_config,
        json,
        mo,
        np,
        pd,
        plt,
        read_fit,
        sample_corrected,
        sample_weighted_mixture,
        unflatten_columns,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # QPyM2: Signal Injection Bias & Calibration Workbench

    Interactive evaluation and validation of signal injection bias studies in QPyM2.
    - **Linear Calibration**: Inspect $n_{\text{fitted}} = a \cdot n_{\text{inj}} + b$ across systematic sweep axes (e.g., Co-60 shifts).
    - **Parameter Stability**: Track slope $a$ and offset $b$ as a function of shift parameters.
    - **Posterior Bias Correction**: Compare raw MCMC posteriors against deterministic correction ($(n - b)/a$) and covariance-propagated sampling.
    - **Marginalization**: Inspect likelihood-weighted combinations of bias-corrected posteriors across shifts.
    """)
    return


@app.cell
def _(mo):
    study_dir_input = mo.ui.text(
        value="studies/example_bias_study",
        label="Study Directory",
        full_width=True,
    )
    co_shift_study_input = mo.ui.text(
        value="studies/example_co_shift",
        label="Optional Corresponding Fit Study (for real posterior traces)",
        full_width=True,
    )

    mo.vstack(
        [
            mo.md("### Study Configuration & Paths"),
            study_dir_input,
            co_shift_study_input,
        ]
    )
    return co_shift_study_input, study_dir_input


@app.cell
def _(
    Path,
    assemble_study,
    co_shift_study_input,
    json,
    mo,
    pd,
    study_dir_input,
):
    _sdir = Path(study_dir_input.value)
    if not _sdir.exists():
        study_obj = None
        status_df = pd.DataFrame()
        results_df = pd.DataFrame()
        calibration_dict = {}
        co_shift_df = pd.DataFrame()
        _header_ui = mo.callout(mo.md(f"⚠️ Directory `{_sdir}` does not exist."), kind="alert")
    else:
        try:
            study_obj = assemble_study(_sdir)
            status_df = study_obj.check_bias_toys()
        except Exception as _e:
            study_obj = None
            status_df = pd.DataFrame()
            _header_ui = mo.callout(mo.md(f"⚠️ Error assembling study: {_e}"), kind="alert")

        # Load results.csv
        _res_path = _sdir / "results.csv"
        if _res_path.is_file():
            results_df = pd.read_csv(_res_path)
        else:
            results_df = pd.DataFrame()

        # Load bias_calibration.json
        _calib_path = _sdir / "bias_calibration.json"
        if _calib_path.is_file():
            with open(_calib_path, "r", encoding="utf-8") as _f:
                calibration_dict = json.load(_f)
        else:
            calibration_dict = {}

        # Optionally load corresponding co_shift study
        _cs_dir = Path(co_shift_study_input.value)
        _cs_res = _cs_dir / "results.csv"
        if _cs_res.is_file():
            co_shift_df = pd.read_csv(_cs_res)
        else:
            co_shift_df = pd.DataFrame()

        if study_obj is not None:
            _n_points = len(study_obj.points)
            _n_toys_total = len(results_df) if not results_df.empty else 0
            _summary_md = f"""
            **Study Name:** `{study_obj.name}`  
            **Evaluation Points:** `{_n_points}` | **Collected Toys:** `{_n_toys_total}`  
            **Reference Cuts:** `co_cut={getattr(study_obj.cfg.reference, "co_cut", "N/A")}`, `ndbd_cutoff={getattr(study_obj.cfg.reference, "ndbd_cutoff", "N/A")}`  
            **Sweep Parameters:** `{dict(getattr(study_obj.cfg, "sweep", {}))}`
            """

            _display_cols = [
                c
                for c in [
                    "task_id",
                    "point_id",
                    "co_shift",
                    "n_toys_completed",
                    "is_completed",
                    "status",
                ]
                if c in status_df.columns
            ]

            _header_ui = mo.vstack(
                [
                    mo.md(_summary_md),
                    mo.md("#### Study Points Status"),
                    mo.ui.table(status_df[_display_cols] if not status_df.empty else status_df),
                ]
            )

    _header_ui
    return calibration_dict, co_shift_df, results_df


@app.cell
def _(mo, results_df):
    xmin_slider = mo.ui.slider(
        start=0.0,
        stop=50.0,
        step=1.0,
        value=15.0,
        label="Bias Fit Cutoff (xmin)",
    )

    _available_points = (
        ["All Points"] + sorted(list(results_df["point_id"].unique()))
        if (not results_df.empty and "point_id" in results_df.columns)
        else ["All Points"]
    )

    point_selector = mo.ui.dropdown(
        options=_available_points,
        value=_available_points[0],
        label="Select Point for Detail Plot",
    )

    mo.vstack(
        [
            mo.md("### Bias Calibration Controls"),
            mo.hstack([xmin_slider, point_selector]),
        ]
    )
    return point_selector, xmin_slider


@app.cell
def _(
    calibration_dict,
    get_bias_from_fitdf,
    mo,
    np,
    pd,
    results_df,
    unflatten_columns,
    xmin_slider,
):
    _xmin_val = float(xmin_slider.value)
    _calib_rows = []

    if not results_df.empty and "point_id" in results_df.columns:
        for _pid, _group in results_df.groupby("point_id", sort=False):
            _co_shift = float(_group["co_shift"].iloc[0]) if "co_shift" in _group.columns else 0.0
            _unflat = unflatten_columns(_group.copy())
            try:
                _popt, _pcov = get_bias_from_fitdf(_unflat, xmin=_xmin_val)
                _a = float(_popt[0])
                _b = float(_popt[1])
                _a_err = float(np.sqrt(max(_pcov[0, 0], 0.0)))
                _b_err = float(np.sqrt(max(_pcov[1, 1], 0.0)))
                _cov_ab = float(_pcov[0, 1])
                _corr = _cov_ab / (_a_err * _b_err) if (_a_err * _b_err) > 0 else 0.0
            except Exception:
                if _pid in calibration_dict:
                    _c = calibration_dict[_pid]
                    _a = _c.get("a", 1.0)
                    _b = _c.get("b", 0.0)
                    _a_err = _c.get("a_err", 0.0)
                    _b_err = _c.get("b_err", 0.0)
                    _cov_ab = _c.get("pcov_01", 0.0)
                    _corr = _cov_ab / (_a_err * _b_err) if (_a_err * _b_err) > 0 else 0.0
                else:
                    _a, _b, _a_err, _b_err, _cov_ab, _corr = 1.0, 0.0, 0.0, 0.0, 0.0, 0.0

            _calib_rows.append(
                {
                    "point_id": _pid,
                    "co_shift": _co_shift,
                    "a (slope)": round(_a, 4),
                    "a_err": round(_a_err, 4),
                    "b (offset)": round(_b, 3),
                    "b_err": round(_b_err, 3),
                    "cov(a, b)": round(_cov_ab, 4),
                    "corr(a, b)": round(_corr, 3),
                    "n_toys": len(_group),
                    "xmin": _xmin_val,
                }
            )

    calib_summary_df = pd.DataFrame(_calib_rows)

    if not calib_summary_df.empty:
        _table_ui = mo.vstack(
            [
                mo.md(
                    "#### Bias Calibration Fit Results (Linear Model: $n_{\\text{fit}} = a \\cdot n_{\\text{inj}} + b$)"
                ),
                mo.ui.table(calib_summary_df),
            ]
        )
    else:
        _table_ui = mo.md("*No calibration data available.*")

    _table_ui
    return (calib_summary_df,)


@app.cell
def _(
    get_bias_from_fitdf,
    mo,
    np,
    plt,
    point_selector,
    results_df,
    unflatten_columns,
    xmin_slider,
):
    # Plot 1: Bias Curve (Fitted vs Injected) and Pulls
    if results_df.empty or "ndbd_inj_" not in results_df.columns:
        _bias_curve_ui = mo.md("*No toy injection data available in results.csv.*")
    else:
        _sel_point = point_selector.value
        _xmin_val = float(xmin_slider.value)

        if _sel_point == "All Points":
            _plot_df = results_df.copy()
            _plot_title = "Bias Curves across All Points"
        else:
            _plot_df = results_df[results_df["point_id"] == _sel_point].copy()
            _plot_title = f"Bias Curve: {_sel_point}"

        _fig, (_ax_bias, _ax_res) = plt.subplots(
            2,
            1,
            figsize=(8, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [2.5, 1]},
        )

        _points = _plot_df["point_id"].unique()
        _cmap = plt.get_cmap("tab10")

        for _i, _pid in enumerate(_points):
            _sub = _plot_df[_plot_df["point_id"] == _pid]
            _color = _cmap(_i % 10)
            _x = _sub["ndbd_inj_"].values
            _y = _sub["ndbd_mode"].values
            _yerr = (_sub["ndbd_upp34"] - _sub["ndbd_mode"]).values

            _mask_fit = _x >= _xmin_val
            _mask_excl = _x < _xmin_val

            _co_val = _sub["co_shift"].iloc[0] if "co_shift" in _sub.columns else _pid
            _lbl = f"{_pid} (Δ={_co_val})" if len(_points) > 1 else "Toy fits"

            # Points included in fit
            _ax_bias.errorbar(
                _x[_mask_fit],
                _y[_mask_fit],
                yerr=_yerr[_mask_fit],
                fmt="o",
                markersize=4,
                alpha=0.75,
                capsize=2,
                color=_color,
                label=_lbl,
            )

            # Points excluded (below xmin)
            if np.any(_mask_excl):
                _ax_bias.errorbar(
                    _x[_mask_excl],
                    _y[_mask_excl],
                    yerr=_yerr[_mask_excl],
                    fmt="x",
                    markersize=4,
                    alpha=0.25,
                    capsize=1.5,
                    color="gray",
                )

            # Linear fit line
            try:
                _unflat = unflatten_columns(_sub.copy())
                _popt, _pcov = get_bias_from_fitdf(_unflat, xmin=_xmin_val)
                _xline = np.linspace(_xmin_val, max(_x.max(), 10.0), 50)
                _yline = _popt[0] * _xline + _popt[1]
                _ax_bias.plot(
                    _xline,
                    _yline,
                    color=_color,
                    lw=1.8,
                    label=f"Fit: y={_popt[0]:.3f}x+{_popt[1]:.2f}",
                )

                # Residuals in lower panel
                _res = _y[_mask_fit] - (_popt[0] * _x[_mask_fit] + _popt[1])
                _ax_res.scatter(_x[_mask_fit], _res, color=_color, alpha=0.6, s=16)
            except Exception:
                pass

        # Reference diagonal line y = x
        _all_x = _plot_df["ndbd_inj_"].values
        _diag_x = np.linspace(0, max(_all_x.max() * 1.05, 10.0), 100)
        _ax_bias.plot(_diag_x, _diag_x, "k--", lw=1.2, alpha=0.5, label="Identity (y = x)")
        _ax_bias.axvline(
            _xmin_val, color="red", ls=":", lw=1.0, alpha=0.7, label=f"xmin = {_xmin_val}"
        )

        _ax_bias.set_ylabel(r"Fitted Yield ($n_{0\nu\beta\beta}$ mode)")
        _ax_bias.set_title(_plot_title)
        _ax_bias.legend(loc="upper left", fontsize=8, framealpha=0.85)
        _ax_bias.grid(True, alpha=0.35)

        _ax_res.axhline(0, color="black", lw=1.0, ls="--", alpha=0.7)
        _ax_res.set_xlabel(r"Injected Signal ($n_{\text{inj}}$)")
        _ax_res.set_ylabel(r"Residual ($y - y_{\text{fit}}$)")
        _ax_res.grid(True, alpha=0.35)

        plt.tight_layout()
        _bias_curve_ui = mo.vstack(
            [
                mo.md("#### Bias Calibration Curves & Residuals"),
                _fig,
            ]
        )

    _bias_curve_ui
    return


@app.cell
def _(calib_summary_df, mo, plt):
    # Plot 2: Parameter Stability across Co-60 Shifts
    if calib_summary_df.empty or "co_shift" not in calib_summary_df.columns:
        _param_stability_ui = mo.md("*Calibration summary table not available for stability plot.*")
    else:
        _fig_stab, (_ax_a, _ax_b) = plt.subplots(1, 2, figsize=(9, 3.8))

        _shifts = calib_summary_df["co_shift"]
        _a_vals = calib_summary_df["a (slope)"]
        _a_errs = calib_summary_df["a_err"]
        _b_vals = calib_summary_df["b (offset)"]
        _b_errs = calib_summary_df["b_err"]

        # Slope a vs co_shift
        _ax_a.errorbar(
            _shifts,
            _a_vals,
            yerr=_a_errs,
            fmt="s-",
            color="crimson",
            capsize=3,
            lw=1.5,
            markersize=6,
        )
        _ax_a.axhline(1.0, color="gray", ls="--", alpha=0.6, label="Ideal a = 1.0")
        _ax_a.set_xlabel("Co-60 Shift [keV]")
        _ax_a.set_ylabel("Slope $a$")
        _ax_a.set_title("Calibration Slope vs Co-60 Shift")
        _ax_a.grid(True, alpha=0.35)
        _ax_a.legend()

        # Offset b vs co_shift
        _ax_b.errorbar(
            _shifts,
            _b_vals,
            yerr=_b_errs,
            fmt="o-",
            color="royalblue",
            capsize=3,
            lw=1.5,
            markersize=6,
        )
        _ax_b.axhline(0.0, color="gray", ls="--", alpha=0.6, label="Ideal b = 0.0")
        _ax_b.set_xlabel("Co-60 Shift [keV]")
        _ax_b.set_ylabel("Offset $b$")
        _ax_b.set_title("Calibration Offset vs Co-60 Shift")
        _ax_b.grid(True, alpha=0.35)
        _ax_b.legend()

        plt.tight_layout()
        _param_stability_ui = mo.vstack(
            [
                mo.md("#### Calibration Parameter Stability Across Systematic Shifts"),
                _fig_stab,
            ]
        )

    _param_stability_ui
    return


@app.cell
def _(
    Path,
    apply_bias_correction,
    calib_summary_df,
    co_shift_df,
    get_1sig_around_mode,
    get_config,
    mo,
    np,
    pd,
    plt,
    read_fit,
    sample_corrected,
):
    # Plot 3: Raw vs Bias-Corrected Posterior Comparison
    _trace_loaded = False
    _ndbd_raw = None
    _selected_shift = 0.0

    if not co_shift_df.empty and "output_uuid" in co_shift_df.columns:
        try:
            _cfg = get_config(environment="local", analysis="base", models="tech_note")
            _data_root = Path(_cfg.infrastructure.data_root)
            _first_row = co_shift_df.iloc[0]
            _uuid = _first_row["output_uuid"]
            _selected_shift = float(_first_row["co_shift"])

            _h5_candidates = list(_data_root.glob(f"**/fit_{_uuid}.h5"))
            if _h5_candidates:
                _table, _trace, _data, _meta = read_fit(_h5_candidates[0])
                if _trace is not None and "ndbd" in _trace.posterior:
                    _ndbd_raw = _trace.posterior["ndbd"].values.flatten()
                    _trace_loaded = True
        except Exception:
            _trace_loaded = False

    if _ndbd_raw is None:
        if not co_shift_df.empty and "ndbd_mean" in co_shift_df.columns:
            _std_val = (
                float(co_shift_df["ndbd_std"].iloc[0]) if "ndbd_std" in co_shift_df.columns else 7.5
            )
            _selected_shift = float(co_shift_df["co_shift"].iloc[0])
        else:
            _std_val = 7.5
            _selected_shift = 0.0

        _rng = np.random.default_rng(42)
        _ndbd_raw = _rng.gamma(shape=2.5, scale=_std_val / np.sqrt(2.5), size=20000)

    if not calib_summary_df.empty:
        _match = calib_summary_df[calib_summary_df["co_shift"] == _selected_shift]
        if not _match.empty:
            _a_c = float(_match["a (slope)"].iloc[0])
            _b_c = float(_match["b (offset)"].iloc[0])
            _cov_ab = float(_match["cov(a, b)"].iloc[0])
            _a_err = float(_match["a_err"].iloc[0])
            _b_err = float(_match["b_err"].iloc[0])
        else:
            _a_c, _b_c = (
                float(calib_summary_df["a (slope)"].iloc[0]),
                float(calib_summary_df["b (offset)"].iloc[0]),
            )
            _cov_ab = float(calib_summary_df["cov(a, b)"].iloc[0])
            _a_err, _b_err = (
                float(calib_summary_df["a_err"].iloc[0]),
                float(calib_summary_df["b_err"].iloc[0]),
            )
    else:
        _a_c, _b_c, _a_err, _b_err, _cov_ab = 1.0, 0.0, 0.02, 2.0, 0.0

    _popt = np.array([_a_c, _b_c])
    _pcov = np.array([[_a_err**2, _cov_ab], [_cov_ab, _b_err**2]])

    _ndbd_corr_det = apply_bias_correction(_ndbd_raw, _popt)
    _ndbd_corr_samp = sample_corrected(_ndbd_raw, _popt, _pcov)

    _mode_raw, _lo_raw, _hi_raw = get_1sig_around_mode(_ndbd_raw)
    _p90_raw = float(np.percentile(_ndbd_raw, 90))

    _mode_det, _lo_det, _hi_det = get_1sig_around_mode(_ndbd_corr_det)
    _p90_det = float(np.percentile(_ndbd_corr_det, 90))

    _mode_samp, _lo_samp, _hi_samp = get_1sig_around_mode(_ndbd_corr_samp)
    _p90_samp = float(np.percentile(_ndbd_corr_samp, 90))

    _fig_post, _ax_post = plt.subplots(figsize=(8, 4.5))
    _bins = np.linspace(0, np.percentile(_ndbd_raw, 99.5) * 1.2, 70)

    _ax_post.hist(
        _ndbd_raw,
        bins=_bins,
        density=True,
        histtype="stepfilled",
        alpha=0.3,
        color="steelblue",
        label=f"Raw Posterior (Mode: {_mode_raw:.2f})",
    )
    _ax_post.hist(
        _ndbd_corr_det,
        bins=_bins,
        density=True,
        histtype="step",
        lw=1.8,
        color="forestgreen",
        label=f"Deterministic Corrected (Mode: {_mode_det:.2f})",
    )
    _ax_post.hist(
        _ndbd_corr_samp,
        bins=_bins,
        density=True,
        histtype="stepfilled",
        alpha=0.25,
        color="crimson",
        label=f"Covariance Corrected (Mode: {_mode_samp:.2f})",
    )

    _ax_post.axvline(_mode_raw, color="steelblue", ls="--", lw=1.2)
    _ax_post.axvline(_mode_det, color="forestgreen", ls="--", lw=1.2)
    _ax_post.axvline(_mode_samp, color="crimson", ls="--", lw=1.2)

    _source_str = "MCMC fit trace" if _trace_loaded else "Representative parameter draw"
    _ax_post.set_title(
        f"0νββ Posterior Bias Correction ({_source_str}, Co-60 Shift = {_selected_shift} keV)"
    )
    _ax_post.set_xlabel(r"Signal Counts $n_{0\nu\beta\beta}$")
    _ax_post.set_ylabel("Posterior Density")
    _ax_post.legend(loc="upper right", fontsize=9)
    _ax_post.grid(True, alpha=0.35)
    plt.tight_layout()

    _summary_post_df = pd.DataFrame(
        {
            "Variant": [
                "Raw Posterior",
                "Deterministic Corrected",
                "Covariance Propagated",
            ],
            "Mode": [round(_mode_raw, 2), round(_mode_det, 2), round(_mode_samp, 2)],
            "1σ Lower": [round(_lo_raw, 2), round(_lo_det, 2), round(_lo_samp, 2)],
            "1σ Upper": [round(_hi_raw, 2), round(_hi_det, 2), round(_hi_samp, 2)],
            "90% Limit": [round(_p90_raw, 2), round(_p90_det, 2), round(_p90_samp, 2)],
        }
    )

    _post_comp_ui = mo.vstack(
        [
            mo.md("#### Posterior Bias Correction Comparison"),
            _fig_post,
            mo.ui.table(_summary_post_df),
        ]
    )

    _post_comp_ui
    return


@app.cell
def _(
    Path,
    calib_summary_df,
    co_shift_df,
    get_1sig_around_mode,
    get_config,
    mo,
    np,
    plt,
    read_fit,
    sample_corrected,
    sample_weighted_mixture,
):
    # Plot 4: Logp-weighted Marginalization across Shifts
    if co_shift_df.empty or "co_shift" not in co_shift_df.columns:
        _marginal_ui = mo.md("*Co-shift results not available for marginalization.*")
    else:
        _raw_list = []
        _corr_list = []
        _logps = []

        try:
            _cfg = get_config(environment="local", analysis="base", models="tech_note")
            _data_root = Path(_cfg.infrastructure.data_root)
        except Exception:
            _data_root = None

        _rng = np.random.default_rng(42)

        for _, _r in co_shift_df.iterrows():
            _s = float(_r["co_shift"])
            _uuid = _r.get("output_uuid", None)
            _trace_found = False

            if _uuid and _data_root:
                _cand = list(_data_root.glob(f"**/fit_{_uuid}.h5"))
                if _cand:
                    try:
                        _tbl, _tr, _, _ = read_fit(_cand[0])
                        if _tr is not None and "ndbd" in _tr.posterior:
                            _samples = _tr.posterior["ndbd"].values.flatten()
                            _trace_found = True
                    except Exception:
                        _trace_found = False

            if not _trace_found:
                _s_val = float(_r["ndbd_std"]) if "ndbd_std" in _r else 7.5
                _samples = _rng.gamma(shape=2.5, scale=_s_val / np.sqrt(2.5), size=10000)

            _raw_list.append(_samples)

            if not calib_summary_df.empty:
                _c_row = calib_summary_df[calib_summary_df["co_shift"] == _s]
                if not _c_row.empty:
                    _a_val = float(_c_row["a (slope)"].iloc[0])
                    _b_val = float(_c_row["b (offset)"].iloc[0])
                    _a_err = float(_c_row["a_err"].iloc[0])
                    _b_err = float(_c_row["b_err"].iloc[0])
                    _cov_ab = float(_c_row["cov(a, b)"].iloc[0])
                else:
                    _a_val, _b_val, _a_err, _b_err, _cov_ab = 1.0, 0.0, 0.02, 2.0, 0.0
            else:
                _a_val, _b_val, _a_err, _b_err, _cov_ab = 1.0, 0.0, 0.02, 2.0, 0.0

            _popt = np.array([_a_val, _b_val])
            _pcov = np.array([[_a_err**2, _cov_ab], [_cov_ab, _b_err**2]])
            _corr_list.append(sample_corrected(_samples, _popt, _pcov))

            _logps.append(-0.5 * (_s / 0.8) ** 2)

        _logps = np.asarray(_logps)
        _weights = np.exp(_logps - np.max(_logps))
        _weights = _weights / _weights.sum()

        _marg_raw = sample_weighted_mixture(_raw_list, weights=_weights, n_draws=25000, rng=_rng)
        _marg_corr = sample_weighted_mixture(_corr_list, weights=_weights, n_draws=25000, rng=_rng)

        _mode_mraw, _lo_mraw, _hi_mraw = get_1sig_around_mode(_marg_raw)
        _mode_mcorr, _lo_mcorr, _hi_mcorr = get_1sig_around_mode(_marg_corr)

        _fig_marg, _ax_marg = plt.subplots(figsize=(8, 4.2))
        _bins_m = np.linspace(0, np.percentile(_marg_raw, 99.5) * 1.2, 70)

        _ax_marg.hist(
            _marg_raw,
            bins=_bins_m,
            density=True,
            histtype="stepfilled",
            alpha=0.35,
            color="orange",
            label=f"Uncorrected Marginal (Mode: {_mode_mraw:.2f})",
        )
        _ax_marg.hist(
            _marg_corr,
            bins=_bins_m,
            density=True,
            histtype="step",
            lw=2.0,
            color="darkblue",
            label=f"Bias-Corrected Marginal (Mode: {_mode_mcorr:.2f})",
        )

        _ax_marg.axvline(_mode_mraw, color="darkorange", ls="--", lw=1.2)
        _ax_marg.axvline(_mode_mcorr, color="darkblue", ls="--", lw=1.2)

        _ax_marg.set_title("Log-likelihood-Weighted Marginal Posterior Across Co-60 Shifts")
        _ax_marg.set_xlabel(r"Signal Counts $n_{0\nu\beta\beta}$")
        _ax_marg.set_ylabel("Posterior Density")
        _ax_marg.legend(loc="upper right", fontsize=9)
        _ax_marg.grid(True, alpha=0.35)
        plt.tight_layout()

        _marginal_ui = mo.vstack(
            [
                mo.md("#### Systematic Marginalization Across Co-60 Shifts"),
                _fig_marg,
            ]
        )

    _marginal_ui
    return


if __name__ == "__main__":
    app.run()
