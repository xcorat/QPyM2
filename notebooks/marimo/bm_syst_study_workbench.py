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
    from matplotlib.patches import Patch
    import marimo as mo
    import numpy as np
    import pandas as pd
    from omegaconf import OmegaConf

    from qpym2.fit.bias import apply_bias_correction
    from qpym2.io.bkg_model import list_bm_variants, read_bm_variant_modes
    from qpym2.storage.backend import read_fit
    from qpym2.storage.paths import get_output_dir
    from qpym2.storage.study import init_study_db
    from qpym2.studies.assembly import assemble_study
    from qpym2.syst.core import SystDelta, combine_deltas, compute_deltas
    from qpym2.syst.io import load_syst_result
    from qpym2.syst.plots import plot_combined_posterior, plot_delta_table, plot_rate_posterior
    from qpym2.syst.rate import F130, N_TEO2, to_rate
    from qpym2.syst.smear import smear_posterior
    from qpym2.utils.mathy import calc_intervals, get_mode

    return (
        Dict,
        F130,
        List,
        N_TEO2,
        OmegaConf,
        Optional,
        Patch,
        Path,
        SystDelta,
        Tuple,
        Union,
        apply_bias_correction,
        assemble_study,
        az,
        calc_intervals,
        combine_deltas,
        compute_deltas,
        get_mode,
        get_output_dir,
        h5py,
        init_study_db,
        json,
        list_bm_variants,
        load_syst_result,
        mo,
        np,
        pd,
        plot_combined_posterior,
        plot_delta_table,
        plot_rate_posterior,
        plt,
        read_bm_variant_modes,
        read_fit,
        smear_posterior,
        to_rate,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # QPyM2: Background Model & Rate Systematics Study Workbench

    Interactive inspection, diagnostics, and validation of Background Model (BM) and Rate Systematics in QPyM2 v0.6:
    - **BM Threshold (`bm_threshold`)**: Evaluates the impact of raising the detector threshold from 200 keV to 250 keV (`Threshold250` vs nominal `Threshold200`). This modifies the M2 fit ROI to $e_1, e_2 \ge 250\text{ keV}$ and overrides MC normalisations with the upstream JAGS modes.
    - **BM Binning (`bm_binning`)**: Evaluates upstream JAGS background model binning choices (`MinBinSize_5_20`, `MinBinSize_25_60`, `BinningIntraAlphaBands`) on nominal geometry. Systematic shifts are aggregated via an **envelope** across variants:
      $$\Delta_{\text{left}} = \min_v \Delta_{\text{left}}^{(v)}, \quad \Delta_{\text{right}} = \max_v \Delta_{\text{right}}^{(v)}$$
    - **Template Likelihood Binning (`bin_size`)**: Evaluates the impact of coarser 2D $(u, v)$ binning (e.g. 2.0 or 5.0 keV) on fitted $0\nu\beta\beta$ rates.
    - **Master Rate Marginalisation (`syst_rate`)**: Combines all 8 systematic deltas in quadrature into split-Gaussian widths ($\sigma_{\text{left}}, \sigma_{\text{right}}$), smears the nominal posterior, convolves with exposure and PSA efficiency uncertainty, and extracts the 90% C.L. half-life limit ($T_{1/2}^{0\nu}$).
    """)
    return


@app.cell
def _(mo):
    study_type_select = mo.ui.dropdown(
        options={
            "studies/example_bm_threshold_study": "BM Threshold Study (Threshold250)",
            "studies/example_bm_binning_study": "BM Binning Study (MinBinSize, IntraAlpha)",
            "studies/example_syst_rate": "Master Systematic Rate Study (Smearing & Limits)",
        },
        value="studies/example_bm_threshold_study",
        label="Select Example Study Preset",
    )
    study_dir_input = mo.ui.text(
        value="studies/example_bm_threshold_study",
        label="Custom Study Directory",
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
    results_raw_input = mo.ui.text(
        value="/data/cuore/Analysis/BMMC/inputs/ResultsRaw.h5",
        label="ResultsRaw.h5 Path (for upstream BM modes)",
        full_width=True,
    )

    mo.vstack(
        [
            mo.md("### Study Configuration & Paths"),
            study_type_select,
            study_dir_input,
            nominal_study_input,
            bias_calib_input,
            results_raw_input,
        ]
    )
    return (
        bias_calib_input,
        nominal_study_input,
        results_raw_input,
        study_dir_input,
        study_type_select,
    )


@app.cell
def _(
    Path,
    assemble_study,
    json,
    mo,
    nominal_study_input,
    pd,
    study_dir_input,
    study_type_select,
):
    _chosen_dir = study_dir_input.value if study_dir_input.value else study_type_select.value
    _sdir = Path(_chosen_dir)

    study_obj = None
    status_df = pd.DataFrame()
    results_df = pd.DataFrame()
    syst_summary_dict = {}
    delta_table_df = pd.DataFrame()
    rate_summary_dict = {}
    nominal_results_df = pd.DataFrame()
    active_category = "bm_threshold"

    if not _sdir.exists():
        _header_ui = mo.callout(mo.md(f"⚠️ Study directory `{_sdir}` does not exist."), kind="alert")
    else:
        try:
            study_obj = assemble_study(_sdir)
            status_df = study_obj.check_syst()
        except Exception as _e:
            _header_ui = mo.callout(
                mo.md(f"⚠️ Note during study assembly / status check: {_e}"), kind="warn"
            )

        _res_path = _sdir / "results.csv"
        results_df = pd.read_csv(_res_path) if _res_path.is_file() else pd.DataFrame()

        _sum_path = _sdir / "syst_summary.json"
        if _sum_path.is_file():
            try:
                with open(_sum_path, "r", encoding="utf-8") as _f:
                    syst_summary_dict = json.load(_f)
            except Exception:
                syst_summary_dict = {}

        _dt_path = _sdir / "delta_table.csv"
        delta_table_df = pd.read_csv(_dt_path) if _dt_path.is_file() else pd.DataFrame()

        _rate_sum_path = _sdir / "syst_rate_summary.json"
        if _rate_sum_path.is_file():
            try:
                with open(_rate_sum_path, "r", encoding="utf-8") as _f:
                    rate_summary_dict = json.load(_f)
            except Exception:
                rate_summary_dict = {}

        # Nominal results if available
        _nom_dir = Path(nominal_study_input.value)
        _nom_res = _nom_dir / "results.csv"
        nominal_results_df = pd.read_csv(_nom_res) if _nom_res.is_file() else pd.DataFrame()

        _scfg = getattr(study_obj.cfg, "syst", None) if study_obj else None
        active_category = (
            str(getattr(_scfg, "category", getattr(study_obj.cfg.study, "type", "syst")))
            if study_obj
            else "unknown"
        )
        if not results_df.empty and "category" in results_df.columns:
            active_category = str(results_df["category"].iloc[0])

        _cat_badge = {
            "bm_threshold": "🟣 **BM Threshold (`bm_threshold`)** — `Threshold250` vs `Threshold200` (ROI $e_1, e_2 \ge 250$ keV)",
            "bm_binning": "🟢 **BM Binning (`bm_binning`)** — Envelope across `MinBinSize` & `IntraAlpha` variants",
            "bin_size": "🟡 **Template Binning (`bin_size`)** — 2D $(u, v)$ likelihood binning",
            "syst_rate": "🔴 **Master Rate (`syst_rate`)** — Full quadrature smearing & limit extraction",
        }.get(active_category, f"🔵 **Systematic Category:** `{active_category}`")

        _m_key = (
            study_obj.points[0].model_key
            if (study_obj and study_obj.points and hasattr(study_obj.points[0], "model_key"))
            else getattr(getattr(study_obj.cfg, "reference", None), "model_key", "default")
            if study_obj
            else "default"
        )

        _summary_md = f"""
        **Study Name:** `{study_obj.name if study_obj else 'N/A'}` | {_cat_badge}  
        **Study Directory:** `{_sdir}`  
        **Evaluation Points:** `{len(study_obj.points) if study_obj else 0}` | **Collected Rows:** `{len(results_df)}`  
        **Reference Model:** `key={_m_key}`, `co_cut={getattr(getattr(study_obj, 'cfg', None).reference, 'co_cut', 'N/A') if study_obj else 'N/A'}`, `ndbd_cutoff={getattr(getattr(study_obj, 'cfg', None).reference, 'ndbd_cutoff', 'N/A') if study_obj else 'N/A'}`  
        **MCMC Settings:** `chains={getattr(_scfg, 'chains', 2) if _scfg else 2}`, `draws={getattr(_scfg, 'draws', 2500) if _scfg else 2500}`, `tune={getattr(_scfg, 'tune', 500) if _scfg else 500}`
        """

        _display_cols = [
            c
            for c in [
                "task_id",
                "point_id",
                "co_shift",
                "n_samples_completed",
                "n_samples_target",
                "is_completed",
                "status",
            ]
            if c in status_df.columns
        ]

        _ui_elements = [mo.md(_summary_md)]
        if not status_df.empty:
            _ui_elements.extend(
                [
                    mo.md("#### Task Completion Matrix"),
                    mo.ui.table(status_df[_display_cols] if not status_df.empty else status_df),
                ]
            )
        _header_ui = mo.vstack(_ui_elements)

    _header_ui
    return (
        active_category,
        delta_table_df,
        nominal_results_df,
        rate_summary_dict,
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
        _calib_path = Path(study_dir_input.value) / "bias_calibration.json"

    popt = None
    calibration_info = None

    if _calib_path.is_file():
        try:
            with open(_calib_path, "r", encoding="utf-8") as _f:
                _calib_data = json.load(_f)

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
def _(mo, results_df, status_df):
    _points = []
    if not status_df.empty and "point_id" in status_df.columns:
        _points = sorted(list(status_df["point_id"].unique()))
    elif not results_df.empty and "point_id" in results_df.columns:
        _points = sorted(list(results_df["point_id"].unique()))

    if not _points:
        _points = ["nominal"]

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
    active_category,
    mo,
    np,
    pd,
    point_selector,
    results_df,
    study_dir_input,
):
    _sdir = Path(study_dir_input.value)
    _pid = point_selector.value
    _task_csv = _sdir / "tasks" / f"syst_{active_category}_{_pid}.csv"

    if _task_csv.is_file():
        point_df = pd.read_csv(_task_csv)
    elif not results_df.empty and "point_id" in results_df.columns:
        point_df = results_df[results_df["point_id"] == _pid].copy().reset_index(drop=True)
    elif not results_df.empty:
        point_df = results_df.copy()
    else:
        point_df = pd.DataFrame()

    has_point_data = not point_df.empty

    if has_point_data:
        _n_rows = len(point_df)
        _vars = list(point_df["variant"].unique()) if "variant" in point_df.columns else []
        _mode_col = (
            "ndbd_mode_bc"
            if ("ndbd_mode_bc" in point_df.columns and not point_df["ndbd_mode_bc"].isna().all())
            else "ndbd_mode"
        )
        _modes = point_df[_mode_col].dropna().values

        _stats_summary = f"""
        ### Per-Variant Fits Table (`{_pid}`)
        **Evaluated Variations:** `{_n_rows}` ({', '.join([f'`{v}`' for v in _vars]) if _vars else 'N/A'})  
        - **$n_{{0\\nu\\beta\\beta}}$ Modes (`{_mode_col}`):** Mean = `{np.mean(_modes):.2f}` ± `{np.std(_modes):.2f}` [Range = `{np.min(_modes):.2f}` to `{np.max(_modes):.2f}`]
        """

        _show_cols = [
            c
            for c in [
                "variant",
                "ndbd_mode",
                "ndbd_p16",
                "ndbd_p90",
                "ndbd_mode_bc",
                "ndbd_p16_bc",
                "ndbd_p90_bc",
                "co60_mode",
                "max_logp",
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
        _table_ui = mo.callout(
            mo.md(
                f"⚠️ No fit tasks found for point `{_pid}` (expected `{_task_csv.name}`). Run systematic fits with `--run-syst` first."
            ),
            kind="warn",
        )

    _table_ui
    return has_point_data, point_df


@app.cell
def _(
    Path,
    active_category,
    apply_bias_correction,
    az,
    calc_intervals,
    compute_deltas,
    get_mode,
    get_output_dir,
    has_point_data,
    mo,
    nominal_results_df,
    nominal_study_input,
    np,
    pd,
    point_df,
    point_selector,
    popt,
    read_fit,
    study_dir_input,
    study_obj,
    syst_summary_dict,
):
    _pid = point_selector.value
    nominal_mode = 0.0
    nominal_post = None
    delta_summary = {}
    df_deltas = pd.DataFrame()

    if has_point_data:
        # 1. Resolve Nominal Posterior & Mode
        _sdir = Path(study_dir_input.value)
        _nom_dir = Path(nominal_study_input.value)
        _nom_post_raw = None

        if _nom_dir.exists():
            _nc_candidates = list((_nom_dir / "tasks").glob(f"*{_pid}*.nc")) + list(
                (_nom_dir / "tasks").glob("*.nc")
            )
            if _nc_candidates:
                try:
                    _trace = az.from_netcdf(str(_nc_candidates[0]))
                    _nom_post_raw = _trace.posterior["ndbd"].values.flatten()
                except Exception:
                    pass

        if _nom_post_raw is None and study_obj is not None:
            try:
                _out_dir = get_output_dir(study_obj.cfg)
                _out_candidates = list(_out_dir.glob("fit_*.h5"))
                if _out_candidates:
                    _fit_data = read_fit(_out_candidates[0])
                    if "trace" in _fit_data and _fit_data["trace"] is not None:
                        _nom_post_raw = _fit_data["trace"].posterior["ndbd"].values.flatten()
            except Exception:
                pass

        if _nom_post_raw is not None:
            nominal_post = (
                apply_bias_correction(_nom_post_raw, popt) if popt is not None else _nom_post_raw
            )
            nominal_mode = float(calc_intervals([nominal_post], method="hdi")[0][0])
        elif not nominal_results_df.empty and "ndbd_mode" in nominal_results_df.columns:
            _raw_mode = float(nominal_results_df["ndbd_mode"].iloc[0])
            nominal_mode = (popt[0] * _raw_mode + popt[1]) if popt is not None else _raw_mode
        else:
            nominal_mode = (
                float(point_df["ndbd_mode"].iloc[0]) if "ndbd_mode" in point_df.columns else 14.64
            )

        # 2. Derive per-variant deltas and envelope
        _mode_col = (
            "ndbd_mode_bc"
            if (
                popt is not None
                and "ndbd_mode_bc" in point_df.columns
                and not point_df["ndbd_mode_bc"].isna().all()
            )
            else "ndbd_mode"
        )

        _rows = []
        _dl_list = []
        _dr_list = []

        # Check if syst_summary.json has precomputed deltas
        _sum_entry = syst_summary_dict.get(_pid, {})
        _pre_dr = _sum_entry.get("delta_right", None)
        _pre_dl = _sum_entry.get("delta_left", None)

        for _, _r in point_df.iterrows():
            _v = str(_r.get("variant", "variant_0"))
            _m = float(_r[_mode_col])
            _dm = _m - nominal_mode

            # Look up in per_variant summary if available
            _pvar = _sum_entry.get("per_variant", {}).get(_v, {})
            _dl_v = _pvar.get("delta_left", min(_dm, 0.0))
            _dr_v = _pvar.get("delta_right", max(_dm, 0.0))

            _dl_list.append(_dl_v)
            _dr_list.append(_dr_v)

            _rows.append(
                {
                    "variant": _v,
                    "mode": _m,
                    "delta_mode": _dm,
                    "delta_left": _dl_v,
                    "delta_right": _dr_v,
                }
            )

        df_deltas = pd.DataFrame(_rows)

        _comb_dl = (
            float(_pre_dl) if _pre_dl is not None else (float(min(_dl_list)) if _dl_list else 0.0)
        )
        _comb_dr = (
            float(_pre_dr) if _pre_dr is not None else (float(max(_dr_list)) if _dr_list else 0.0)
        )

        delta_summary = {
            "nominal_mode": nominal_mode,
            "delta_right": _comb_dr,
            "delta_left": _comb_dl,
            "category": active_category,
            "n_variants": len(df_deltas),
        }

    return delta_summary, df_deltas, nominal_mode, nominal_post


@app.cell
def _(active_category, delta_summary, mo, nominal_mode, point_selector):
    if not delta_summary:
        _cards_ui = mo.md("")
    else:
        _pid = point_selector.value
        _dr = delta_summary["delta_right"]
        _dl = delta_summary["delta_left"]
        _dr_pct = (_dr / nominal_mode * 100.0) if nominal_mode > 0 else 0.0
        _dl_pct = (_dl / nominal_mode * 100.0) if nominal_mode > 0 else 0.0

        _stat_cards = mo.hstack(
            [
                mo.stat(
                    value=f"{nominal_mode:.2f}",
                    label=f"Nominal Mode ({_pid})",
                    caption="Reference HDI mode counts",
                ),
                mo.stat(
                    value=f"+{_dr:.3f}",
                    label="+σ Systematic (Δ_right)",
                    caption=f"+{_dr_pct:.2f}% of nominal mode",
                ),
                mo.stat(
                    value=f"{_dl:.3f}",
                    label="-σ Systematic (Δ_left)",
                    caption=f"{_dl_pct:.2f}% of nominal mode",
                ),
                mo.stat(
                    value=f"[{nominal_mode + _dl:.2f}, {nominal_mode + _dr:.2f}]",
                    label="Systematic Mode Interval",
                    caption="Envelope [m_nom + Δ_left, m_nom + Δ_right]",
                ),
            ],
            justify="start",
            gap=2,
        )

        _cards_ui = mo.vstack(
            [
                mo.md(f"### Systematic Shift Summary: `{active_category}`"),
                _stat_cards,
            ]
        )

    _cards_ui
    return


@app.cell
def _(
    active_category,
    df_deltas,
    has_point_data,
    mo,
    np,
    plt,
    point_selector,
):
    if not has_point_data or df_deltas.empty:
        _plot_ui = mo.md("No delta shift data available to plot.")
    else:
        _pid = point_selector.value
        _fig, _ax = plt.subplots(figsize=(9, 4.5))

        _v_names = df_deltas["variant"].values
        _deltas = df_deltas["delta_mode"].values
        _x = np.arange(len(_v_names))
        _bw = 0.45

        _colors = ["#ff7f0e" if d >= 0 else "#1f77b4" for d in _deltas]
        _ax.bar(_x, _deltas, _bw, color=_colors, alpha=0.85, edgecolor="black", linewidth=0.8)

        _ax.axhline(0, color="black", linestyle="--", linewidth=1.0)
        _ax.set_xticks(_x)
        _ax.set_xticklabels(_v_names, rotation=30, ha="right", fontsize=10)
        _ax.set_xlabel("Systematic Variant", fontsize=11)
        _ax.set_ylabel(
            r"$\Delta_{\mathrm{mode}} = \mathrm{mode}_{\mathrm{syst}} - \mathrm{mode}_{\mathrm{nom}}$ (counts)",
            fontsize=11,
        )
        _ax.set_title(rf"Systematic Mode Shifts ({active_category}) — Point: {_pid}", fontsize=12)
        _ax.grid(axis="y", linestyle=":", alpha=0.6)

        for _xi, _d in zip(_x, _deltas):
            _va = "bottom" if _d >= 0 else "top"
            _offset = 0.05 * (1 if _d >= 0 else -1)
            _ax.text(
                _xi, _d + _offset, f"{_d:+.2f}", ha="center", va=_va, fontsize=9, fontweight="bold"
            )

        plt.tight_layout()
        _plot_ui = mo.vstack(
            [
                mo.md("#### Variant Shift Comparison Bar Plot"),
                _fig,
            ]
        )

    _plot_ui
    return


@app.cell
def _(
    Path,
    active_category,
    delta_summary,
    h5py,
    has_point_data,
    mo,
    nominal_mode,
    nominal_post,
    np,
    plt,
    point_selector,
    study_dir_input,
):
    if not has_point_data:
        _overlay_ui = mo.md("")
    else:
        _sdir = Path(study_dir_input.value)
        _pid = point_selector.value
        _h5_samples_path = _sdir / "tasks" / f"syst_{active_category}_{_pid}_samples.h5"

        if not _h5_samples_path.is_file():
            _overlay_ui = mo.callout(
                mo.md(
                    f"ℹ️ Posterior samples file `{_h5_samples_path.name}` not found. Task metrics available in table above."
                ),
                kind="info",
            )
        else:
            with h5py.File(_h5_samples_path, "r") as _hf:
                _samps_arr = _hf["ndbd_samples"][:]

            _fig, _ax = plt.subplots(figsize=(10, 5))

            # Plot variant posteriors
            _colors = plt.cm.Set1(np.linspace(0, 1, len(_samps_arr)))
            for _i, (_samps, _col) in enumerate(zip(_samps_arr, _colors)):
                _lbl = f"Variant {_i}"
                _ax.hist(
                    _samps,
                    bins=40,
                    density=True,
                    histtype="step",
                    color=_col,
                    alpha=0.85,
                    linewidth=1.5,
                    label=_lbl,
                )

            # Plot nominal posterior
            if nominal_post is not None:
                _ax.hist(
                    nominal_post,
                    bins=50,
                    density=True,
                    histtype="stepfilled",
                    color="black",
                    alpha=0.25,
                    label="Nominal Posterior",
                )

            _dr = delta_summary.get("delta_right", 0.0)
            _dl = delta_summary.get("delta_left", 0.0)

            _ax.axvline(
                nominal_mode,
                color="black",
                linestyle="--",
                linewidth=1.5,
                label=f"Nominal Mode: {nominal_mode:.2f}",
            )
            _ax.axvspan(
                nominal_mode + _dl,
                nominal_mode + _dr,
                color="gold",
                alpha=0.3,
                label=f"Systematic Interval: [{nominal_mode + _dl:.2f}, {nominal_mode + _dr:.2f}]",
            )

            _ax.set_xlabel(r"$n_{0\nu\beta\beta}$ (counts)", fontsize=11)
            _ax.set_ylabel("Posterior Density", fontsize=11)
            _ax.set_title(
                rf"Posterior Overlays & Systematic Shift Band ({active_category}) — {_pid}",
                fontsize=12,
            )
            _ax.legend(fontsize=9, loc="upper right")
            _ax.grid(alpha=0.3)
            plt.tight_layout()

            _overlay_ui = mo.vstack(
                [
                    mo.md("#### Posterior Sample Overlays & Shift Interval"),
                    _fig,
                ]
            )

    _overlay_ui
    return


@app.cell
def _(Path, list_bm_variants, mo, results_raw_input):
    _raw_path = Path(results_raw_input.value)
    available_variants = []
    if _raw_path.is_file():
        try:
            available_variants = list_bm_variants(_raw_path)
        except Exception:
            available_variants = [
                "Threshold250",
                "MinBinSize_5_20",
                "MinBinSize_25_60",
                "BinningIntraAlphaBands",
            ]
    else:
        available_variants = [
            "Threshold250",
            "MinBinSize_5_20",
            "MinBinSize_25_60",
            "BinningIntraAlphaBands",
        ]

    bm_comp_variant_selector = mo.ui.dropdown(
        options=available_variants,
        value="Threshold250"
        if "Threshold250" in available_variants
        else (available_variants[0] if available_variants else ""),
        label="Select BM Variant to Inspect Upstream Normalisations",
    )

    mo.vstack(
        [
            mo.md("---"),
            mo.md("### Upstream Background Model Modes (`ResultsRaw.h5`)"),
            bm_comp_variant_selector,
        ]
    )
    return available_variants, bm_comp_variant_selector


@app.cell
def _(
    Path,
    bm_comp_variant_selector,
    mo,
    pd,
    plt,
    read_bm_variant_modes,
    results_raw_input,
):
    _raw_path = Path(results_raw_input.value)
    _sel_variant = bm_comp_variant_selector.value

    if not _raw_path.is_file() or not _sel_variant:
        _bm_modes_ui = mo.md("_ResultsRaw.h5 not found or no variant selected._")
    else:
        try:
            modes_series = read_bm_variant_modes(_raw_path, _sel_variant)
            df_modes = pd.DataFrame(
                {"Component": modes_series.index, "Mode": modes_series.values}
            ).sort_values("Mode", ascending=False)

            top_df = df_modes.head(12)
            fig_bm, ax_bm = plt.subplots(figsize=(8, 4))
            ax_bm.barh(top_df["Component"][::-1], top_df["Mode"][::-1], color="#3366cc", alpha=0.85)
            ax_bm.set_xlabel("Fitted Component Mode (Counts)")
            ax_bm.set_title(f"Upstream Mode Normalisations: {_sel_variant} (Top 12 Components)")
            ax_bm.grid(axis="x", linestyle="--", alpha=0.5)
            plt.tight_layout()

            _bm_modes_ui = mo.vstack(
                [
                    mo.ui.table(df_modes.head(15)),
                    fig_bm,
                ]
            )
        except Exception as _e:
            _bm_modes_ui = mo.callout(
                mo.md(f"⚠️ Could not read modes for variant `{_sel_variant}`: {_e}"), kind="alert"
            )

    _bm_modes_ui
    return ax_bm, df_modes, fig_bm, modes_series, top_df


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ### Interactive Rate Marginalisation & Half-Life Calculator

    Configure physical parameters to combine systematic deltas in quadrature, perform split-Gaussian posterior smearing, convolve with exposure and PSA efficiency, and calculate the 90% C.L. half-life limit ($T_{1/2}^{0\nu}$).
    """)
    return


@app.cell
def _(mo):
    exposure_slider = mo.ui.slider(
        start=500.0,
        stop=2000.0,
        step=50.0,
        value=1038.4,
        label="Exposure (kg·yr)",
    )
    eff_mc_slider = mo.ui.slider(
        start=0.01,
        stop=0.10,
        step=0.005,
        value=0.040,
        label="Signal Selection Efficiency (eff_mc)",
    )
    psa_mean_slider = mo.ui.slider(
        start=0.90,
        stop=1.00,
        step=0.01,
        value=0.99,
        label="PSA Efficiency Mean",
    )
    psa_sig_slider = mo.ui.slider(
        start=0.00,
        stop=0.05,
        step=0.005,
        value=0.01,
        label="PSA Efficiency Sigma (1σ)",
    )
    n_smear_slider = mo.ui.slider(
        start=10000,
        stop=200000,
        step=10000,
        value=50000,
        label="Number of Smearing Draws (N)",
    )

    mo.vstack(
        [
            exposure_slider,
            eff_mc_slider,
            psa_mean_slider,
            psa_sig_slider,
            n_smear_slider,
        ]
    )
    return (
        eff_mc_slider,
        exposure_slider,
        n_smear_slider,
        psa_mean_slider,
        psa_sig_slider,
    )


@app.cell
def _(
    SystDelta,
    calc_intervals,
    combine_deltas,
    delta_table_df,
    eff_mc_slider,
    exposure_slider,
    get_mode,
    mo,
    n_smear_slider,
    np,
    pd,
    plot_combined_posterior,
    plot_delta_table,
    plot_rate_posterior,
    psa_mean_slider,
    psa_sig_slider,
    smear_posterior,
    to_rate,
):
    if not delta_table_df.empty and "delta_left" in delta_table_df.columns:
        deltas_list = [
            SystDelta(
                name=str(row["source"]),
                delta_right=float(row["delta_right"]),
                delta_left=float(row["delta_left"]),
            )
            for _, row in delta_table_df.iterrows()
        ]
    else:
        default_deltas = [
            ("lss_dq", -3.59, 4.74),
            ("lss_sscl", -1.20, 1.20),
            ("bmmc", -1.50, 0.50),
            ("bm_threshold", -1.16, 0.00),
            ("bm_binning", -0.60, 0.00),
            ("ndbd_shift", -0.42, 0.00),
            ("co_shift", -0.80, 0.80),
            ("bias_corr", -0.04, 0.30),
        ]
        deltas_list = [
            SystDelta(name=sname, delta_right=dr, delta_left=dl) for sname, dl, dr in default_deltas
        ]

    sigma_left, sigma_right = combine_deltas([(d.delta_left, d.delta_right) for d in deltas_list])

    rng = np.random.default_rng(42)
    nominal_bc = np.clip(rng.normal(loc=14.64, scale=4.8, size=25000), 0.0, None)

    smeared_post = smear_posterior(
        nominal_bc,
        sigma_right=float(abs(sigma_right)),
        sigma_left=float(abs(sigma_left)),
        n_smear=int(n_smear_slider.value),
        rng_seed=42,
    )

    rate_post = to_rate(
        smeared_post,
        eff_mc=float(eff_mc_slider.value),
        exposure=float(exposure_slider.value),
        psa_mean=float(psa_mean_slider.value),
        psa_sig=float(psa_sig_slider.value),
        rng_seed=42,
        clip=True,
    )

    mode_nom = float(get_mode(nominal_bc))
    mode_smeared = float(get_mode(smeared_post))
    rate_mode = float(get_mode(rate_post))
    rate_90 = float(np.percentile(rate_post, 90))
    half_life_90 = float(np.log(2.0) / rate_90) if rate_90 > 0 else np.inf

    fig_dt, _, _ = plot_delta_table(deltas_list, sigma_right=sigma_right, sigma_left=sigma_left)
    fig_comb, _ = plot_combined_posterior(
        marg_post=smeared_post,
        nominal_bc=nominal_bc,
        sigma_right=sigma_right,
        sigma_left=sigma_left,
    )
    fig_rate, _ = plot_rate_posterior(
        rate=rate_post,
        stats={
            "mode": rate_mode,
            "p90": rate_90,
            "hdi": (
                float(calc_intervals([rate_post])[0][1]),
                float(calc_intervals([rate_post])[0][2]),
            ),
            "hdi_prob": 0.68,
        },
    )

    results_banner = f"""
    ### Marginalisation Summary
    - **Total Systematic Width:** $\\sigma_{{\\text{{left}}}} = {sigma_left:+.2f}$, $\\sigma_{{\\text{{right}}}} = {sigma_right:+.2f}$
    - **Fitted $0\\nu\\beta\\beta$ Mode:** Nominal $= {mode_nom:.2f}$ counts $\\longrightarrow$ Smeared $= {mode_smeared:.2f}$ counts
    - **Decay Rate Limit (90% C.L.):** $\\Gamma_{{90}} = {rate_90:.2e}\\text{{ yr}}^{{-1}}$
    - **Half-Life Limit (90% C.L.):** $T_{{1/2}}^{{0\\nu}} > \\mathbf{{{half_life_90:.2e}}}\\text{{ yr}}$
    """

    _marginalisation_ui = mo.vstack(
        [
            mo.md(results_banner),
            mo.md("#### 1. Systematics Breakdown & Quadrature Sum"),
            fig_dt,
            mo.md("#### 2. Posterior Smearing"),
            fig_comb,
            mo.md("#### 3. Marginalised Decay Rate Posterior"),
            fig_rate,
        ]
    )

    _marginalisation_ui
    return (
        deltas_list,
        fig_comb,
        fig_dt,
        fig_rate,
        half_life_90,
        mode_nom,
        mode_smeared,
        nominal_bc,
        rate_90,
        rate_mode,
        rate_post,
        results_banner,
        rng,
        sigma_left,
        sigma_right,
        smeared_post,
    )


@app.cell
def _(
    Path,
    active_category,
    has_point_data,
    load_syst_result,
    mo,
    point_selector,
    study_dir_input,
):
    if not has_point_data:
        _inspector_ui = mo.md("")
    else:
        _sdir = Path(study_dir_input.value)
        _pid = point_selector.value
        _h5_candidates = [
            _sdir / f"syst_{active_category}_{_pid}.h5",
            _sdir / "combined_result.h5",
        ]
        _h5_res = next((p for p in _h5_candidates if p.is_file()), None)

        if _h5_res is not None:
            try:
                _loaded = load_syst_result(_h5_res)
                _meta = _loaded.get("meta", {})
                _per_post = _loaded.get("per_post", None)
                _metrics = _loaded.get("metrics", {})

                _details = f"""
                ---
                ### Systematic Result HDF5 (Format v2)
                **File:** `{_h5_res.name}`  
                - **Name:** `{_loaded.get('name')}` | **Format Version:** `{_meta.get('format_version', 'N/A')}`  
                - **Study / Point:** `{_meta.get('study_name', 'N/A')}` / `{_meta.get('point_id', 'N/A')}`  
                - **Summary Metrics:** `delta_right = {_metrics.get('delta_right', 'N/A')}`, `delta_left = {_metrics.get('delta_left', 'N/A')}`
                """

                _inspector_ui = mo.vstack(
                    [
                        mo.md(_details),
                        mo.md("#### Per-Variation Records (`per_post`)")
                        if _per_post is not None
                        else mo.md(""),
                        mo.ui.table(_per_post) if _per_post is not None else mo.md(""),
                    ]
                )
            except Exception as _e:
                _inspector_ui = mo.callout(
                    mo.md(f"⚠️ Error reading Format v2 HDF5 `{_h5_res.name}`: {_e}"), kind="warn"
                )
        else:
            _inspector_ui = mo.callout(
                mo.md(
                    f"ℹ️ Result artifact for `{active_category}` not created yet. Run `--collect-syst` to generate Format v2 outputs."
                ),
                kind="info",
            )

    _inspector_ui
    return


if __name__ == "__main__":
    app.run()
