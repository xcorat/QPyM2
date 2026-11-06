import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    from pathlib import Path
    from typing import Dict, List, Optional, Tuple, Union

    import arviz as az
    import h5py
    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    import pandas as pd
    from scipy.stats import pearsonr

    from qpym2.fit.bias import apply_bias_correction
    from qpym2.storage.backend import read_fit
    from qpym2.storage.inputs import resolve_storage_coords
    from qpym2.storage.paths import get_output_dir
    from qpym2.studies.assembly import assemble_study
    from qpym2.utils.mathy import calc_intervals, get_mode, sample_weighted_mixture

    return (
        Dict,
        List,
        Optional,
        Path,
        Tuple,
        Union,
        apply_bias_correction,
        assemble_study,
        az,
        calc_intervals,
        get_mode,
        get_output_dir,
        h5py,
        json,
        mo,
        np,
        pd,
        pearsonr,
        plt,
        read_fit,
        resolve_storage_coords,
        sample_weighted_mixture,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # QPyM2: BMMC (Stat-Only) Systematic Study Workbench

    Interactive inspection, diagnostics, and validation of Background Model Monte Carlo (BMMC) systematic fits in QPyM2 v0.6:
    - **Stat-Only Mode**: Background component normalisation vectors are drawn directly and unperturbed from the raw JAGS MCMC chain ($\vec{N}_{\text{sample}} \in \text{mkchain}$), isolating pure background statistical fluctuations without covariance perturbation.
    - **MCMC Fit Diagnostics**: Inspect distributions of $0\nu\beta\beta$ mode (`ndbd_mode`), quantiles (`p16`, `p90`), log-likelihood scores (`max_logp`), and background component modes ($^{60}\text{Co}$, $^{238}\text{U}$, $^{232}\text{Th}$, others).
    - **Weighting Schemes & Collapse Diagnostic**: Compare Equal Weights ($1/N$, recommended for posterior draws) vs Log-Likelihood Softmax Weights, tracking the **Effective Sample Size** ($N_{\text{eff}} = 1/\sum w_i^2$) to detect weight collapse.
    - **Posterior Reconstruction & Nominal Comparison**: Reconstruct the mixture posterior, apply linear bias calibration if available, and extract $\Delta_{\text{mode}}$, $\Delta_{\text{left}}$, and $\Delta_{\text{right}}$ relative to the nominal fit.
    - **Sample-Size Convergence**: Track stability of systematic shifts as a function of sample count $k = 1, \dots, N$.
    - **Sensitivity Analysis**: Correlate sampled background normalizations against $n_{0\nu\beta\beta}$ to identify leading drivers of systematic shifts.
    """)
    return


@app.cell
def _(mo):
    study_dir_input = mo.ui.text(
        value="studies/example_bmmc_study",
        label="BMMC Study Directory",
        full_width=True,
    )
    nominal_study_input = mo.ui.text(
        value="studies/example_co_shift",
        label="Optional Reference / Data Study (for nominal fit trace)",
        full_width=True,
    )
    bias_calib_input = mo.ui.text(
        value="studies/example_bias_study/bias_calibration.json",
        label="Optional Bias Calibration JSON Path",
        full_width=True,
    )

    mo.vstack(
        [
            mo.md("### Study Configuration & Paths"),
            study_dir_input,
            nominal_study_input,
            bias_calib_input,
        ]
    )
    return bias_calib_input, nominal_study_input, study_dir_input


@app.cell
def _(
    Path,
    assemble_study,
    json,
    mo,
    nominal_study_input,
    pd,
    study_dir_input,
):
    _sdir = Path(study_dir_input.value)
    if not _sdir.exists():
        study_obj = None
        status_df = pd.DataFrame()
        results_df = pd.DataFrame()
        syst_summary_dict = {}
        nominal_results_df = pd.DataFrame()
        _header_ui = mo.callout(mo.md(f"⚠️ Study directory `{_sdir}` does not exist."), kind="alert")
    else:
        try:
            study_obj = assemble_study(_sdir)
            status_df = study_obj.check_syst()
        except Exception as _e:
            study_obj = None
            status_df = pd.DataFrame()
            _header_ui = mo.callout(mo.md(f"⚠️ Error assembling study: {_e}"), kind="alert")

        # Load results.csv
        _res_path = _sdir / "results.csv"
        results_df = pd.read_csv(_res_path) if _res_path.is_file() else pd.DataFrame()

        # Load syst_summary.json
        _sum_path = _sdir / "syst_summary.json"
        if _sum_path.is_file():
            with open(_sum_path, "r", encoding="utf-8") as _f:
                syst_summary_dict = json.load(_f)
        else:
            syst_summary_dict = {}

        # Load nominal study results if available
        _nom_dir = Path(nominal_study_input.value)
        _nom_res = _nom_dir / "results.csv"
        nominal_results_df = pd.read_csv(_nom_res) if _nom_res.is_file() else pd.DataFrame()

        if study_obj is not None:
            _scfg = getattr(study_obj.cfg, "syst", None)
            _stat_only_flag = bool(getattr(_scfg, "stat_only", False)) if _scfg else False
            if not _stat_only_flag and not results_df.empty and "stat_only" in results_df.columns:
                _stat_only_flag = bool(results_df["stat_only"].iloc[0])

            _n_points = len(study_obj.points)
            _n_samples_total = len(results_df) if not results_df.empty else 0
            _stat_mode_badge = (
                "🟢 **Stat-Only Mode (True)** — unperturbed JAGS raw chain sampling"
                if _stat_only_flag
                else "🟠 **Covariance Perturbation Mode (False)** — multivariate normal sampling"
            )

            _m_key = (
                study_obj.points[0].model_key
                if (study_obj.points and hasattr(study_obj.points[0], "model_key"))
                else getattr(getattr(study_obj.cfg, "reference", None), "model_key", "default")
            )

            _summary_md = f"""
            **Study Name:** `{study_obj.name}` | **Category:** `{getattr(_scfg, "category", "bmmc")}`  
            **Sampling Mode:** {_stat_mode_badge}  
            **Evaluation Points:** `{_n_points}` | **Collected MCMC Fits:** `{_n_samples_total}`  
            **Reference Model:** `key={_m_key}`, `co_cut={getattr(study_obj.cfg.reference, "co_cut", "N/A")}`, `ndbd_cutoff={getattr(study_obj.cfg.reference, "ndbd_cutoff", "N/A")}`, `co_shift={getattr(study_obj.cfg.reference, "co_shift", 0.0)}`  
            **MCMC Settings:** `chains={getattr(_scfg, "chains", 2)}`, `draws={getattr(_scfg, "draws", 2500)}`, `tune={getattr(_scfg, "tune", 500)}`
            """

            _display_cols = [
                c
                for c in [
                    "task_id",
                    "point_id",
                    "co_shift",
                    "n_samples_completed",
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
    return (
        nominal_results_df,
        results_df,
        status_df,
        study_obj,
        syst_summary_dict,
    )


@app.cell
def _(
    Path,
    bias_calib_input,
    json,
    mo,
    study_dir_input,
    study_obj,
):
    _calib_path = Path(bias_calib_input.value)
    if not _calib_path.is_file():
        # Fallback to study_dir / bias_calibration.json
        _calib_path = Path(study_dir_input.value) / "bias_calibration.json"

    popt = None
    calibration_info = None

    if _calib_path.is_file():
        try:
            with open(_calib_path, "r", encoding="utf-8") as _f:
                _calib_data = json.load(_f)

            # Match reference co_shift or point_id
            _ref_shift = (
                float(getattr(getattr(study_obj, "cfg", None).reference, "co_shift", 1.0))
                if study_obj
                else 1.0
            )
            _matched_entry = None

            for _k, _v in _calib_data.items():
                if isinstance(_v, dict) and float(_v.get("co_shift", -999.0)) == _ref_shift:
                    _matched_entry = _v
                    break

            if _matched_entry is None and "nominal" in _calib_data:
                _matched_entry = _calib_data["nominal"]
            elif _matched_entry is None and len(_calib_data) > 0:
                _first_val = next(iter(_calib_data.values()))
                if isinstance(_first_val, dict) and "a" in _first_val:
                    _matched_entry = _first_val

            if _matched_entry and "a" in _matched_entry and "b" in _matched_entry:
                _a = float(_matched_entry["a"])
                _b = float(_matched_entry["b"])
                popt = (_a, _b)
                calibration_info = _matched_entry
                _calib_ui = mo.callout(
                    mo.md(
                        f"✓ **Linear Bias Calibration Loaded** from `{_calib_path.name}`: "
                        f"**slope $a = {_a:.4f}$**, **offset $b = {_b:.4f}$** "
                        f"(matched for `co_shift = {_ref_shift:.1f}`)."
                    ),
                    kind="success",
                )
            else:
                _calib_ui = mo.callout(
                    mo.md(
                        f"⚠️ Calibration file `{_calib_path}` found, but missing valid `a` and `b` parameters."
                    ),
                    kind="warn",
                )
        except Exception as _e:
            _calib_ui = mo.callout(mo.md(f"⚠️ Error reading calibration JSON: {_e}"), kind="alert")
    else:
        _calib_ui = mo.callout(
            mo.md("ℹ️ **Bias calibration not found.** Displaying raw / uncorrected posteriors."),
            kind="info",
        )

    _calib_ui
    return calibration_info, popt


@app.cell
def _(mo, results_df):
    _points = (
        sorted(list(results_df["point_id"].unique()))
        if (not results_df.empty and "point_id" in results_df.columns)
        else ["nominal"]
    )

    point_selector = mo.ui.dropdown(
        options=_points,
        value=_points[0],
        label="Select Study Point to Inspect",
    )

    point_selector
    return (point_selector,)


@app.cell
def _(
    Path,
    mo,
    np,
    pd,
    point_selector,
    results_df,
    study_dir_input,
):
    _sdir = Path(study_dir_input.value)
    _pid = point_selector.value

    # Load point task CSV
    _task_csv = _sdir / "tasks" / f"syst_mc_{_pid}.csv"
    if _task_csv.is_file():
        point_df = pd.read_csv(_task_csv)
    elif not results_df.empty and "point_id" in results_df.columns:
        point_df = results_df[results_df["point_id"] == _pid].copy().reset_index(drop=True)
    else:
        point_df = pd.DataFrame()

    if not point_df.empty:
        _n_fits = len(point_df)
        _mode_arr = (
            point_df["ndbd_mode"].dropna().values
            if "ndbd_mode" in point_df.columns
            else np.array([])
        )
        _logp_arr = (
            point_df["max_logp"].dropna().values if "max_logp" in point_df.columns else np.array([])
        )

        _stats_summary = f"""
        ### MCMC Fit Diagnostics Table ({_pid})
        **Total Completed Fits:** `{_n_fits}`  
        - **$0\nu\beta\beta$ Mode (`ndbd_mode`):** Mean = `{np.mean(_mode_arr):.2f}` ± `{np.std(_mode_arr):.2f}`, Median = `{np.median(_mode_arr):.2f}` [Min = `{np.min(_mode_arr):.2f}`, Max = `{np.max(_mode_arr):.2f}`]  
        - **Log-Likelihood Mode (`max_logp`):** Mean = `{np.mean(_logp_arr):.2f}`, Max = `{np.max(_logp_arr):.2f}`, Spread = `{np.max(_logp_arr) - np.min(_logp_arr):.2f}`
        """

        _show_cols = [
            c
            for c in [
                "sample_idx",
                "ndbd_mode",
                "ndbd_p16",
                "ndbd_p90",
                "max_logp",
                "co60_mode",
                "u238_mode",
                "th232_mode",
                "others_mode",
            ]
            if c in point_df.columns
        ]

        _table_ui = mo.vstack(
            [
                mo.md(_stats_summary),
                mo.ui.table(point_df[_show_cols]),
            ]
        )
    else:
        _table_ui = mo.callout(mo.md(f"⚠️ No fit tasks found for point `{_pid}`."), kind="warn")

    _table_ui
    return (point_df,)


@app.cell
def _(mo, np, plt, point_df):
    if point_df.empty or "ndbd_mode" not in point_df.columns:
        _dist_ui = mo.md("No fit data available for distribution plots.")
    else:
        _fig, _axes = plt.subplots(1, 3, figsize=(15, 4))

        # 1. ndbd_mode distribution
        _modes = point_df["ndbd_mode"].values
        _axes[0].hist(
            _modes,
            bins=min(max(len(_modes) // 2, 5), 30),
            color="#1f77b4",
            alpha=0.7,
            edgecolor="black",
        )
        _axes[0].axvline(
            np.median(_modes), color="red", linestyle="--", label=f"Median: {np.median(_modes):.2f}"
        )
        _axes[0].set_xlabel(r"$0\nu\beta\beta$ Mode (counts)")
        _axes[0].set_ylabel("Number of Fits")
        _axes[0].set_title(r"Fitted $n_{0\nu\beta\beta}$ Modes Across Samples")
        _axes[0].legend(fontsize=9)
        _axes[0].grid(alpha=0.3)

        # 2. max_logp distribution
        _logps = point_df["max_logp"].values
        _axes[1].hist(
            _logps,
            bins=min(max(len(_logps) // 2, 5), 30),
            color="#2ca02c",
            alpha=0.7,
            edgecolor="black",
        )
        _axes[1].axvline(
            np.max(_logps), color="purple", linestyle="--", label=f"Max: {np.max(_logps):.2f}"
        )
        _axes[1].set_xlabel(r"$\ln \mathcal{L}_{\mathrm{max}}$")
        _axes[1].set_ylabel("Number of Fits")
        _axes[1].set_title("Log-Likelihood Mode Distribution")
        _axes[1].legend(fontsize=9)
        _axes[1].grid(alpha=0.3)

        # 3. Background Component Modes
        _bkg_cols = [
            c
            for c in ["co60_mode", "u238_mode", "th232_mode", "others_mode"]
            if c in point_df.columns
        ]
        if _bkg_cols:
            for _c in _bkg_cols:
                _axes[2].hist(
                    point_df[_c].values, bins=15, alpha=0.5, label=_c.replace("_mode", "")
                )
            _axes[2].set_xlabel("Fitted Component Mode (counts)")
            _axes[2].set_ylabel("Number of Fits")
            _axes[2].set_title("Background Component Modes")
            _axes[2].legend(fontsize=8)
            _axes[2].grid(alpha=0.3)
        else:
            _axes[2].axis("off")

        plt.tight_layout()

        # Outlier Detection Callout
        _outliers = []
        _mean_mode, _std_mode = np.mean(_modes), np.std(_modes)
        if _std_mode > 0:
            _outlier_mask = np.abs(_modes - _mean_mode) > 2.5 * _std_mode
            if np.any(_outlier_mask):
                _outlier_indices = list(point_df.loc[_outlier_mask, "sample_idx"].values)
                _outliers.append(f"Samples with extreme mode (>2.5σ): `{_outlier_indices}`")

        _max_logp_val = np.max(_logps)
        _bad_logp_mask = (_max_logp_val - _logps) > 25.0
        if np.any(_bad_logp_mask):
            _bad_logp_idx = list(point_df.loc[_bad_logp_mask, "sample_idx"].values)
            _outliers.append(
                f"Samples with exceptionally low likelihood (ΔlnL > 25): `{_bad_logp_idx}`"
            )

        if _outliers:
            _alert_ui = mo.callout(
                mo.md("⚠️ **Potential Outlier Fits Detected:**\n- " + "\n- ".join(_outliers)),
                kind="warn",
            )
        else:
            _alert_ui = mo.callout(
                mo.md(
                    "✓ **MCMC Fit Stability Verified:** No severe mode or likelihood outliers detected."
                ),
                kind="success",
            )

        _dist_ui = mo.vstack(
            [
                mo.md("### Sample Mode & Likelihood Distributions"),
                mo.as_html(_fig),
                _alert_ui,
            ]
        )
        plt.close(_fig)

    _dist_ui
    return


@app.cell
def _(mo):
    weighting_selector = mo.ui.radio(
        options=[
            "Equal Weights (1/N, recommended for stat-only)",
            "Softmax Log-Likelihood Mode (pik)",
        ],
        value="Equal Weights (1/N, recommended for stat-only)",
        label="Mixture Weighting Method",
    )

    show_individual_traces = mo.ui.checkbox(
        value=False,
        label="Overlay individual sample posteriors (alpha = 0.15)",
    )

    mo.vstack(
        [
            mo.md("### Posterior Mixture & Weighting Controls"),
            mo.hstack([weighting_selector, show_individual_traces]),
        ]
    )
    return show_individual_traces, weighting_selector


@app.cell
def _(mo, np, plt, point_df, weighting_selector):
    if point_df.empty or "max_logp" not in point_df.columns:
        weights = np.array([])
        _weight_ui = mo.md("No likelihood data available.")
    else:
        _n = len(point_df)
        _logps = point_df["max_logp"].values

        if "Equal" in weighting_selector.value:
            weights = np.full(_n, 1.0 / _n)
            _method_desc = "Equal weights $w_i = 1/N$ (pure posterior MCMC draws from JAGS chain)."
        else:
            _scores = _logps - np.max(_logps)
            _w = np.exp(_scores)
            _w[_w < 1e-10] = 0.0
            weights = _w / _w.sum() if _w.sum() > 0 else np.full(_n, 1.0 / _n)
            _method_desc = r"Softmax mode log-likelihood weights $w_i \propto \exp(\ln \mathcal{L}_i - \ln \mathcal{L}_{\max})$."

        _n_eff = 1.0 / np.sum(weights**2)
        _max_w = np.max(weights)

        # Plot weight bar chart
        _fig, _ax = plt.subplots(figsize=(10, 3))
        _ax.bar(np.arange(_n), weights, color="#ff7f0e", alpha=0.8, edgecolor="black")
        _ax.set_xlabel("Sample Index")
        _ax.set_ylabel("Weight $w_i$")
        _ax.set_title(
            f"Sample Mixture Weights (Effective Sample Size $N_{{\\mathrm{{eff}}}} = {_n_eff:.1f} / {_n}$)"
        )
        _ax.grid(alpha=0.3)
        plt.tight_layout()

        if "Softmax" in weighting_selector.value and (_n_eff < 2.0 or (_n_eff / _n) < 0.25):
            _warning_box = mo.callout(
                mo.md(
                    f"⚠️ **Weight Collapse Detected!**  \n"
                    f"The effective sample size is **$N_{{\\mathrm{{eff}}}} = {_n_eff:.1f}$** out of **{_n}** fits "
                    f"(the single best fit carries **{_max_w * 100:.1f}%** of the entire mixture weight).  \n"
                    f"*Why this happens:* Because individual M2 fits differ by large log-likelihood values, "
                    f"naive softmax exponentiation suppresses almost all MCMC samples.  \n"
                    f"**Recommendation:** Use **Equal Weights ($1/N$)**, which is mathematically proper for "
                    f"unperturbed MCMC draws from the JAGS posterior."
                ),
                kind="alert",
            )
        else:
            _warning_box = mo.callout(
                mo.md(
                    f"✓ **Weighting Health:** Effective sample size $N_{{\\mathrm{{eff}}}} = {_n_eff:.1f}$ fits."
                ),
                kind="success",
            )

        _weight_ui = mo.vstack(
            [
                mo.md(f"**Selected Weighting Scheme:** {_method_desc}"),
                mo.as_html(_fig),
                _warning_box,
            ]
        )
        plt.close(_fig)

    _weight_ui
    return (weights,)


@app.cell
def _(
    Path,
    apply_bias_correction,
    az,
    calc_intervals,
    get_output_dir,
    h5py,
    mo,
    nominal_results_df,
    nominal_study_input,
    np,
    plt,
    point_selector,
    popt,
    resolve_storage_coords,
    sample_weighted_mixture,
    show_individual_traces,
    study_dir_input,
    study_obj,
    weights,
):
    _sdir = Path(study_dir_input.value)
    _pid = point_selector.value
    _h5_path = _sdir / "tasks" / f"syst_mc_{_pid}_samples.h5"

    combined_syst_post = np.array([])
    nom_post = np.array([])

    if not _h5_path.is_file():
        _comp_ui = mo.callout(
            mo.md(f"⚠️ Task sample HDF5 `{_h5_path}` not found. Run fits or collection first."),
            kind="alert",
        )
        delta_mode = 0.0
        delta_right = 0.0
        delta_left = 0.0
    else:
        with h5py.File(_h5_path, "r") as _h5f:
            _ndbd_samples = _h5f["ndbd_samples"][:]

        # Resample mixture
        _rng = np.random.default_rng(42)
        combined_syst_post = sample_weighted_mixture(
            _ndbd_samples, weights, n_draws=50000, rng=_rng
        )

        # Resolve nominal posterior
        # 1. Check storage output directory
        if study_obj is not None:
            try:
                _m_key = (
                    study_obj.points[0].model_key
                    if (study_obj.points and hasattr(study_obj.points[0], "model_key"))
                    else getattr(getattr(study_obj.cfg, "reference", None), "model_key", "default")
                )
                _coords, _hm, _param_dir, _input_dir = resolve_storage_coords(
                    study_obj.cfg, model_key=_m_key
                )
                _out_dir = get_output_dir(_param_dir)
                # Find matching output_uuid from nominal_results_df
                if not nominal_results_df.empty and "output_uuid" in nominal_results_df.columns:
                    _ref_shift = float(
                        getattr(getattr(study_obj, "cfg", None).reference, "co_shift", 1.0)
                    )
                    _match = nominal_results_df[nominal_results_df["co_shift"] == _ref_shift]
                    if not _match.empty:
                        _uuid = _match["output_uuid"].iloc[0]
                        _nc_candidate = _out_dir / f"fit_{_uuid}.nc"
                        if _nc_candidate.is_file():
                            _idata = az.from_netcdf(str(_nc_candidate))
                            nom_post = _idata.posterior["ndbd"].values.flatten()
            except Exception:
                pass

        # 2. Check nominal study dir directly
        if len(nom_post) == 0:
            _nom_sdir = Path(nominal_study_input.value)
            _candidate_nc = list(_nom_sdir.glob("*.nc"))
            if _candidate_nc:
                try:
                    _idata = az.from_netcdf(str(_candidate_nc[0]))
                    nom_post = _idata.posterior["ndbd"].values.flatten()
                except Exception:
                    pass

        # 3. Fallback to sample 0
        if len(nom_post) == 0:
            nom_post = _ndbd_samples[0]
            _nom_source_note = (
                "⚠️ Storage nominal fit not found; using Sample 0 as reference fallback."
            )
        else:
            _nom_source_note = (
                f"✓ Resolved nominal fit trace from storage ({len(nom_post):,} MCMC draws)."
            )

        # Apply bias correction if available
        if popt is not None:
            _nom_plot = apply_bias_correction(nom_post, popt)
            _syst_plot = apply_bias_correction(combined_syst_post, popt)
            _label_suffix = " (Bias-Corrected)"
        else:
            _nom_plot = nom_post
            _syst_plot = combined_syst_post
            _label_suffix = " (Raw)"

        # Intervals & Deltas
        _mode_nom, _lo_nom, _hi_nom = calc_intervals([_nom_plot], method="hdi")[0]
        _mode_syst, _lo_syst, _hi_syst = calc_intervals([_syst_plot], method="hdi")[0]

        delta_mode = float(_mode_syst - _mode_nom)
        delta_right = float(max(0.0, _hi_syst - _hi_nom))
        delta_left = float(min(0.0, _lo_syst - _lo_nom))

        # Render Comparison Plot
        _fig, _ax = plt.subplots(figsize=(10, 5))

        # Individual sample traces if requested
        if show_individual_traces.value:
            for _s_idx in range(len(_ndbd_samples)):
                _s_trace = _ndbd_samples[_s_idx]
                if popt is not None:
                    _s_trace = apply_bias_correction(_s_trace, popt)
                _ax.hist(
                    _s_trace,
                    bins=60,
                    density=True,
                    histtype="step",
                    color="gray",
                    alpha=0.15,
                    linewidth=0.8,
                )

        # Nominal fit
        _ax.hist(
            _nom_plot,
            bins=70,
            density=True,
            histtype="step",
            color="#1f77b4",
            linewidth=2.2,
            label=f"Nominal Fit{_label_suffix}",
        )
        _ax.axvline(
            _mode_nom,
            color="#1f77b4",
            linestyle="--",
            linewidth=1.5,
            label=f"Nominal Mode: {_mode_nom:.2f}",
        )
        _ax.axvspan(
            _lo_nom,
            _hi_nom,
            color="#1f77b4",
            alpha=0.12,
            label=f"Nominal 68% HDI: [{_lo_nom:.2f}, {_hi_nom:.2f}]",
        )

        # Systematic combined mixture
        _ax.hist(
            _syst_plot,
            bins=70,
            density=True,
            histtype="step",
            color="#ff7f0e",
            linewidth=2.2,
            label=f"BMMC Stat-Only Mixture{_label_suffix}",
        )
        _ax.axvline(
            _mode_syst,
            color="#ff7f0e",
            linestyle="--",
            linewidth=1.5,
            label=f"Syst Mode: {_mode_syst:.2f}",
        )
        _ax.axvspan(
            _lo_syst,
            _hi_syst,
            color="#ff7f0e",
            alpha=0.12,
            label=f"Syst 68% HDI: [{_lo_syst:.2f}, {_hi_syst:.2f}]",
        )

        _ax.set_xlabel(r"$0\nu\beta\beta$ Signal Counts ($n_{0\nu\beta\beta}$)")
        _ax.set_ylabel("Posterior Probability Density")
        _ax.set_title(
            f"BMMC Systematic vs Nominal Posterior ({_pid})\n"
            f"Δmode = {delta_mode:+.2f} counts | Δleft = {delta_left:+.2f} | Δright = {delta_right:+.2f}"
        )
        _ax.legend(fontsize=9, loc="upper right")
        _ax.grid(alpha=0.3)
        plt.tight_layout()

        _kpi_md = f"""
        | Metric | Nominal Fit | BMMC (Stat-Only) Mixture | Systematic Shift ($\Delta$) |
        | :--- | :---: | :---: | :---: |
        | **Mode** | `{_mode_nom:.2f}` counts | `{_mode_syst:.2f}` counts | **$\Delta_{{\\mathrm{{mode}}}} = {delta_mode:+.2f}$ counts** |
        | **68.27% HDI Interval** | `[{_lo_nom:.2f}, {_hi_nom:.2f}]` | `[{_lo_syst:.2f}, {_hi_syst:.2f}]` | **$\Delta_{{\\mathrm{{left}}}} = {delta_left:+.2f}$**, **$\Delta_{{\\mathrm{{right}}}} = {delta_right:+.2f}$** |
        """

        _comp_ui = mo.vstack(
            [
                mo.md("### Posterior Mixture vs Nominal Comparison"),
                mo.md(f"*{_nom_source_note}*"),
                mo.as_html(_fig),
                mo.md(_kpi_md),
            ]
        )
        plt.close(_fig)

    _comp_ui
    return (
        combined_syst_post,
        delta_left,
        delta_mode,
        delta_right,
        nom_post,
    )


@app.cell
def _(
    Path,
    apply_bias_correction,
    calc_intervals,
    h5py,
    mo,
    np,
    plt,
    point_df,
    point_selector,
    popt,
    sample_weighted_mixture,
    study_dir_input,
    weights,
):
    _sdir = Path(study_dir_input.value)
    _pid = point_selector.value
    _h5_path = _sdir / "tasks" / f"syst_mc_{_pid}_samples.h5"

    if not _h5_path.is_file() or point_df.empty:
        _conv_ui = mo.md("No sample data for convergence analysis.")
    else:
        with h5py.File(_h5_path, "r") as _h5f:
            _samples_all = _h5f["ndbd_samples"][:]

        _n_total = len(_samples_all)
        if _n_total < 3:
            _conv_ui = mo.callout(
                mo.md(f"Need at least 3 samples to evaluate convergence (currently {_n_total})."),
                kind="info",
            )
        else:
            _rng = np.random.default_rng(123)
            _running_k = []
            _running_mode_syst = []

            # Evaluate progressive running mixture modes
            for _k in range(2, _n_total + 1):
                _sub_samples = _samples_all[:_k]
                _sub_w = weights[:_k]
                _sub_w = _sub_w / _sub_w.sum() if _sub_w.sum() > 0 else np.full(_k, 1.0 / _k)
                _mix_k = sample_weighted_mixture(_sub_samples, _sub_w, n_draws=20000, rng=_rng)
                if popt is not None:
                    _mix_k = apply_bias_correction(_mix_k, popt)
                _mode_k, _, _ = calc_intervals([_mix_k], method="hdi")[0]
                _running_k.append(_k)
                _running_mode_syst.append(_mode_k)

            _fig, _ax = plt.subplots(figsize=(10, 3.5))
            _ax.plot(_running_k, _running_mode_syst, marker="o", color="#d62728", linewidth=2.0)
            _ax.set_xlabel("Number of Accumulated MCMC Sample Fits ($k$)")
            _ax.set_ylabel(r"Mixture Mode $n_{0\nu\beta\beta}$ (counts)")
            _ax.set_title("Systematic Mode Convergence vs Sample Size ($N$-dependence)")
            _ax.grid(alpha=0.3)
            plt.tight_layout()

            _conv_ui = mo.vstack(
                [
                    mo.md("### Sample-Size ($N$) Convergence Diagnostics"),
                    mo.md(
                        "Track whether the mixture mode has stabilized as additional sample fits are accumulated."
                    ),
                    mo.as_html(_fig),
                ]
            )
            plt.close(_fig)

    _conv_ui
    return


@app.cell
def _(mo, np, pearsonr, plt, point_df):
    _bkg_names = [
        c for c in ["co60_mode", "u238_mode", "th232_mode", "others_mode"] if c in point_df.columns
    ]
    if point_df.empty or "ndbd_mode" not in point_df.columns or not _bkg_names:
        _sens_ui = mo.md("No background component data available for sensitivity analysis.")
    else:
        _fig, _axes = plt.subplots(1, len(_bkg_names), figsize=(4 * len(_bkg_names), 3.5))
        if len(_bkg_names) == 1:
            _axes = [_axes]

        _corrs_text = []
        for _idx, _bkg in enumerate(_bkg_names):
            _x = point_df[_bkg].values
            _y = point_df["ndbd_mode"].values
            _r, _pval = pearsonr(_x, _y) if len(_x) > 2 else (0.0, 1.0)
            _bkg_clean = _bkg.replace("_mode", "")
            _corrs_text.append(f"**{_bkg_clean}:** $r = {_r:+.3f}$ ($p = {_pval:.3f}$)")

            _axes[_idx].scatter(_x, _y, color="#1f77b4", alpha=0.8, edgecolors="black")
            # Trendline
            if len(_x) > 1 and np.std(_x) > 0:
                _pfit = np.polyfit(_x, _y, 1)
                _x_line = np.linspace(np.min(_x), np.max(_x), 50)
                _axes[_idx].plot(
                    _x_line, np.polyval(_pfit, _x_line), color="red", linestyle="--", linewidth=1.2
                )

            _axes[_idx].set_xlabel(f"{_bkg_clean} Mode (counts)")
            _axes[_idx].set_ylabel(r"$n_{0\nu\beta\beta}$ Mode (counts)")
            _axes[_idx].set_title(f"{_bkg_clean} vs $0\\nu\\beta\\beta$ ($r={_r:+.2f}$)")
            _axes[_idx].grid(alpha=0.3)

        plt.tight_layout()

        _sens_ui = mo.vstack(
            [
                mo.md("### Sensitivity Analysis: Background Normalization Correlations"),
                mo.md(
                    "Identifies which background model component normalisation most strongly shifts the $0\nu\beta\beta$ signal:  \n"
                    + " | ".join(_corrs_text)
                ),
                mo.as_html(_fig),
            ]
        )
        plt.close(_fig)

    _sens_ui
    return


if __name__ == "__main__":
    app.run()
