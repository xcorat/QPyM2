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

    from qpym2.fit.bias import apply_bias_correction
    from qpym2.storage.backend import read_fit
    from qpym2.storage.inputs import resolve_storage_coords
    from qpym2.storage.paths import get_output_dir
    from qpym2.studies.assembly import assemble_study
    from qpym2.syst.plots import plot_lss_sscl_barplot, plot_lss_dq_barplot
    from qpym2.syst.io import load_syst_result
    from qpym2.utils.mathy import calc_intervals, get_mode, signed_quad_sum

    return (
        Dict,
        List,
        Optional,
        Patch,
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
        load_syst_result,
        mo,
        np,
        pd,
        plot_lss_dq_barplot,
        plot_lss_sscl_barplot,
        plt,
        read_fit,
        resolve_storage_coords,
        signed_quad_sum,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # QPyM2: Lineshape Systematic (LSS) Study Workbench

    Interactive inspection, diagnostics, and validation of Lineshape Systematics (**Energy Scale `lss_sscl`** and **Energy Resolution `lss_dq`**) in QPyM2 v0.6:
    - **Tech Note Method 1 Analysis Logic**: Per-dataset MCMC fits across the 15 CUORE datasets ($3601 \le \text{ds} \le 3615$), bias correction on posterior draws, and signed quadrature summation of mode shifts $\Delta_{ds}$:
      $$\delta_{\text{right}} = +\sqrt{\sum_{ds} (\max(\Delta_{ds}, 0))^2}, \quad \delta_{\text{left}} = -\sqrt{\sum_{ds} (\min(\Delta_{ds}, 0))^2}$$
    - **Resolution Envelopes (`lss_dq`)**: Evaluates asymmetric variations between $+1\sigma$ (`1r`) and $-1\sigma$ (`1l`) calibration curves per dataset, extracting positive and negative extrema envelopes before quadrature summation.
    - **Per-Dataset Shift Diagnostics**: Inspect mode shifts $\Delta_{ds}$, HDI intervals ($p_{16}, p_{90}$), fit likelihoods ($\ln \mathcal{L}_{\text{max}}$), and $^{60}\text{Co}$ normalizations across all datasets.
    - **Dataset Dominance & Pareto Breakdown**: Identify which datasets are the primary drivers of energy scale or resolution systematic uncertainties.
    - **Posterior Overlays**: Direct visualization of shifted dataset posteriors compared to the nominal posterior and the composite uncertainty band $[m_{\text{nom}} + \delta_{\text{left}}, m_{\text{nom}} + \delta_{\text{right}}]$.
    - **Format v2 Artifact Inspector**: Validate checkpointed task CSVs, samples HDF5, Format v2 result files (`per_post`, `per_fit`, `metrics`), and SQLite study records.
    """)
    return


@app.cell
def _(mo):
    study_dir_input = mo.ui.text(
        value="studies/example_lss_sscl_study",
        label="LSS Study Directory",
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
        is_lss_dq = False
        dq_mode_val = "pm1"
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

        _scfg = getattr(study_obj.cfg, "syst", None) if study_obj else None
        _cat = str(getattr(_scfg, "category", "lss_sscl")) if _scfg else "lss_sscl"
        if not results_df.empty and "category" in results_df.columns:
            _cat = str(results_df["category"].iloc[0])

        is_lss_dq = _cat in ("dq", "lss_dq")
        dq_mode_val = str(getattr(_scfg, "dq_mode", "pm1")) if _scfg else "pm1"

        if study_obj is not None:
            _n_points = len(study_obj.points)
            _n_samples_total = len(results_df) if not results_df.empty else 0
            _cat_badge = (
                "🟡 **Energy Resolution (`lss_dq`)** — "
                f"Mode `{dq_mode_val}` (envelope across ±1σ calibration curves)"
                if is_lss_dq
                else "🔵 **Energy Scale (`lss_sscl`)** — Single-variant scaling (`sscale_1r`)"
            )

            _m_key = (
                study_obj.points[0].model_key
                if (study_obj.points and hasattr(study_obj.points[0], "model_key"))
                else getattr(getattr(study_obj.cfg, "reference", None), "model_key", "default")
            )

            _summary_md = f"""
            **Study Name:** `{study_obj.name}` | **Systematic Type:** {_cat_badge}  
            **Evaluation Points:** `{_n_points}` | **Collected Rows:** `{_n_samples_total}`  
            **Reference Coordinates:** `key={_m_key}`, `co_cut={getattr(study_obj.cfg.reference, "co_cut", "N/A")}`, `ndbd_cutoff={getattr(study_obj.cfg.reference, "ndbd_cutoff", "N/A")}`, `co_shift={getattr(study_obj.cfg.reference, "co_shift", 0.0)}`  
            **MCMC Sampler:** `chains={getattr(_scfg, "chains", 2)}`, `draws={getattr(_scfg, "draws", 2500)}`, `tune={getattr(_scfg, "tune", 500)}`
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

            _header_ui = mo.vstack(
                [
                    mo.md(_summary_md),
                    mo.md("#### Study Task Completion Status"),
                    mo.ui.table(status_df[_display_cols] if not status_df.empty else status_df),
                ]
            )

    _header_ui
    return (
        dq_mode_val,
        is_lss_dq,
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
    is_lss_dq,
    mo,
    np,
    pd,
    point_selector,
    results_df,
    study_dir_input,
):
    _sdir = Path(study_dir_input.value)
    _pid = point_selector.value
    _cat_tag = "lss_dq" if is_lss_dq else "lss_sscl"

    # Attempt to load point task CSV
    _task_csv = _sdir / "tasks" / f"syst_{_cat_tag}_{_pid}.csv"
    if _task_csv.is_file():
        point_df = pd.read_csv(_task_csv)
    elif not results_df.empty and "point_id" in results_df.columns:
        point_df = results_df[results_df["point_id"] == _pid].copy().reset_index(drop=True)
    else:
        point_df = pd.DataFrame()

    has_point_data = not point_df.empty

    if has_point_data:
        _n_rows = len(point_df)
        _ds_unique = sorted(point_df["ds"].unique().tolist()) if "ds" in point_df.columns else []
        _mode_col = (
            "ndbd_mode_bc"
            if ("ndbd_mode_bc" in point_df.columns and not point_df["ndbd_mode_bc"].isna().all())
            else "ndbd_mode"
        )
        _modes = point_df[_mode_col].dropna().values

        _stats_summary = f"""
        ### Per-Dataset Fits Table ({_pid})
        **Completed Variations:** `{_n_rows}` across `{len(_ds_unique)}` datasets ({_ds_unique[0] if _ds_unique else 'N/A'}–{_ds_unique[-1] if _ds_unique else 'N/A'})  
        - **$n_{{0\\nu\\beta\\beta}}$ Modes (`{_mode_col}`):** Mean = `{np.mean(_modes):.2f}` ± `{np.std(_modes):.2f}`, Median = `{np.median(_modes):.2f}` [Range = `{np.min(_modes):.2f}` to `{np.max(_modes):.2f}`]
        """

        _show_cols = [
            c
            for c in [
                "ds",
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
    apply_bias_correction,
    az,
    calc_intervals,
    get_mode,
    get_output_dir,
    has_point_data,
    is_lss_dq,
    mo,
    nominal_results_df,
    nominal_study_input,
    np,
    pd,
    point_df,
    point_selector,
    popt,
    read_fit,
    signed_quad_sum,
    study_dir_input,
    study_obj,
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

        # Try nominal trace from reference study or local storage
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

        # Fallback to local storage output
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

        # Fallback to nominal_results_df or synthetic baseline if trace not accessible
        if _nom_post_raw is not None:
            nominal_post = (
                apply_bias_correction(_nom_post_raw, popt) if popt is not None else _nom_post_raw
            )
            nominal_mode = float(calc_intervals([nominal_post], method="hdi")[0][0])
        elif not nominal_results_df.empty and "ndbd_mode" in nominal_results_df.columns:
            _raw_mode = float(nominal_results_df["ndbd_mode"].iloc[0])
            nominal_mode = popt[0] * _raw_mode + popt[1] if popt is not None else _raw_mode
        else:
            # Baseline from first row or median as fallback
            nominal_mode = (
                float(point_df["ndbd_mode"].median()) if "ndbd_mode" in point_df.columns else 10.0
            )

        # 2. Compute per-dataset shifts
        _use_bc = (
            popt is not None
            and "ndbd_mode_bc" in point_df.columns
            and not point_df["ndbd_mode_bc"].isna().all()
        )
        _mode_col = "ndbd_mode_bc" if _use_bc else "ndbd_mode"

        _datasets = sorted(list(point_df["ds"].unique()))
        _rows = []

        if not is_lss_dq:
            # Energy Scale (lss_sscl)
            _dmodes = []
            for _d in _datasets:
                _sub = point_df[point_df["ds"] == _d]
                if not _sub.empty:
                    _m = float(_sub[_mode_col].iloc[0])
                    _delta = _m - nominal_mode
                    _dmodes.append(_delta)
                    _rows.append(
                        {
                            "ds": _d,
                            "variant": "sscale_1r",
                            "mode": _m,
                            "delta_mode": _delta,
                            "delta_right": max(_delta, 0.0),
                            "delta_left": min(_delta, 0.0),
                        }
                    )

            _d_arr = np.array(_dmodes)
            _dr, _dl = signed_quad_sum(_d_arr)
            delta_summary = {
                "nominal_mode": nominal_mode,
                "delta_right": float(_dr),
                "delta_left": float(_dl),
                "category": "lss_sscl",
                "n_datasets": len(_datasets),
            }
        else:
            # Energy Resolution (lss_dq)
            _dpos_list = []
            _dneg_list = []
            for _d in _datasets:
                _sub = point_df[point_df["ds"] == _d]
                _row_1r = _sub[_sub["variant"].str.contains("1r", na=False)]
                _row_1l = _sub[_sub["variant"].str.contains("1l", na=False)]

                _m_1r = float(_row_1r[_mode_col].iloc[0]) if not _row_1r.empty else nominal_mode
                _m_1l = float(_row_1l[_mode_col].iloc[0]) if not _row_1l.empty else nominal_mode

                _d1r = _m_1r - nominal_mode
                _d1l = _m_1l - nominal_mode

                _dpos = max(_d1r, _d1l, 0.0)
                _dneg = min(_d1r, _d1l, 0.0)
                _dpos_list.append(_dpos)
                _dneg_list.append(_dneg)

                _rows.append(
                    {
                        "ds": _d,
                        "variant": "dqres_pm1",
                        "mode_1r": _m_1r,
                        "mode_1l": _m_1l,
                        "delta_1r": _d1r,
                        "delta_1l": _d1l,
                        "delta_right": _dpos,
                        "delta_left": _dneg,
                    }
                )

            _dr = float(np.sqrt(np.sum(np.array(_dpos_list) ** 2)))
            _dl = -float(np.sqrt(np.sum(np.array(_dneg_list) ** 2)))
            delta_summary = {
                "nominal_mode": nominal_mode,
                "delta_right": _dr,
                "delta_left": _dl,
                "category": "lss_dq",
                "n_datasets": len(_datasets),
            }

        df_deltas = pd.DataFrame(_rows)

    return (
        delta_summary,
        df_deltas,
        nominal_mode,
        nominal_post,
    )


@app.cell
def _(delta_summary, is_lss_dq, mo, nominal_mode, point_selector):
    if not delta_summary:
        _cards_ui = mo.md("")
    else:
        _pid = point_selector.value
        _dr = delta_summary["delta_right"]
        _dl = delta_summary["delta_left"]
        _cat_title = "Energy Resolution (lss_dq)" if is_lss_dq else "Energy Scale (lss_sscl)"

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
                    label="+σ Systematic (δ_right)",
                    caption=f"+{_dr_pct:.2f}% of nominal mode",
                ),
                mo.stat(
                    value=f"{_dl:.3f}",
                    label="-σ Systematic (δ_left)",
                    caption=f"{_dl_pct:.2f}% of nominal mode",
                ),
                mo.stat(
                    value=f"[{nominal_mode + _dl:.2f}, {nominal_mode + _dr:.2f}]",
                    label="Systematic Mode Interval",
                    caption="Envelope [m_nom + δ_left, m_nom + δ_right]",
                ),
            ],
            justify="start",
            gap=2,
        )

        _cards_ui = mo.vstack(
            [
                mo.md(f"### Systematic Mode Shift Summary: {_cat_title}"),
                _stat_cards,
            ]
        )

    _cards_ui
    return


@app.cell
def _(
    df_deltas,
    has_point_data,
    is_lss_dq,
    mo,
    np,
    plot_lss_dq_barplot,
    plot_lss_sscl_barplot,
    plt,
    point_selector,
):
    if not has_point_data or df_deltas.empty:
        _plot_ui = mo.md("No delta shift data available to plot.")
    else:
        _pid = point_selector.value
        _fig, _ax = plt.subplots(figsize=(13, 5))

        _ds_ids = df_deltas["ds"].values
        _x = np.arange(len(_ds_ids))

        if not is_lss_dq:
            # Energy scale (single variant)
            _deltas = df_deltas["delta_mode"].values
            _bw = 0.55

            _colors = ["#ff7f0e" if d >= 0 else "#1f77b4" for d in _deltas]
            _bars = _ax.bar(
                _x, _deltas, _bw, color=_colors, alpha=0.85, edgecolor="black", linewidth=0.8
            )

            _ax.axhline(0, color="black", linestyle="--", linewidth=1.0)
            _ax.set_xticks(_x)
            _ax.set_xticklabels([f"ds{d}" for d in _ds_ids], rotation=45, ha="right", fontsize=9)
            _ax.set_xlabel("CUORE Dataset", fontsize=11)
            _ax.set_ylabel(
                r"$\Delta_{\mathrm{mode}} = \mathrm{mode}_{ds} - \mathrm{mode}_{\mathrm{nom}}$ (counts)",
                fontsize=11,
            )
            _ax.set_title(
                rf"Energy Scale Systematic Shifts by Dataset ($\mathtt{{lss\_sscl}}$) — Point: {_pid}",
                fontsize=12,
            )
            _ax.grid(axis="y", linestyle=":", alpha=0.6)

            # Annotate delta values above/below bars
            for _xi, _d in zip(_x, _deltas):
                _va = "bottom" if _d >= 0 else "top"
                _offset = 0.04 * (1 if _d >= 0 else -1)
                _ax.text(
                    _xi,
                    _d + _offset,
                    f"{_d:+.2f}",
                    ha="center",
                    va=_va,
                    fontsize=8,
                    fontweight="bold",
                )
        else:
            # Energy resolution (dqres_1r and dqres_1l side-by-side)
            _d1r = df_deltas["delta_1r"].values
            _d1l = df_deltas["delta_1l"].values
            _bw = 0.35

            _b1 = _ax.bar(
                _x - _bw / 2,
                _d1r,
                _bw,
                label="dQ +1σ (1r)",
                color="#1f77b4",
                alpha=0.8,
                edgecolor="black",
                linewidth=0.8,
            )
            _b2 = _ax.bar(
                _x + _bw / 2,
                _d1l,
                _bw,
                label="dQ -1σ (1l)",
                color="#2ca02c",
                alpha=0.8,
                edgecolor="black",
                linewidth=0.8,
            )

            _ax.axhline(0, color="black", linestyle="--", linewidth=1.0)
            _ax.set_xticks(_x)
            _ax.set_xticklabels([f"ds{d}" for d in _ds_ids], rotation=45, ha="right", fontsize=9)
            _ax.set_xlabel("CUORE Dataset", fontsize=11)
            _ax.set_ylabel(
                r"$\Delta_{\mathrm{mode}} = \mathrm{mode}_{ds} - \mathrm{mode}_{\mathrm{nom}}$ (counts)",
                fontsize=11,
            )
            _ax.set_title(
                rf"Energy Resolution Systematic Shifts ($\mathtt{{lss\_dq}}$) — Point: {_pid}",
                fontsize=12,
            )
            _ax.legend(fontsize=9, loc="upper right")
            _ax.grid(axis="y", linestyle=":", alpha=0.6)

        plt.tight_layout()
        _plot_ui = mo.vstack(
            [
                mo.md("#### Per-Dataset Mode Shifts (Tech Note Method 1)"),
                _fig,
            ]
        )

    _plot_ui
    return


@app.cell
def _(df_deltas, has_point_data, is_lss_dq, mo, np, plt, point_selector):
    if not has_point_data or df_deltas.empty:
        _pareto_ui = mo.md("")
    else:
        _pid = point_selector.value
        _fig, _axes = plt.subplots(1, 2, figsize=(14, 4.5))

        # 1. Positive shifts (contributing to delta_right)
        _df_pos = df_deltas[df_deltas["delta_right"] > 0].copy()
        if not _df_pos.empty:
            _df_pos["var_contrib"] = _df_pos["delta_right"] ** 2
            _tot_var_pos = _df_pos["var_contrib"].sum()
            _df_pos["frac_pos"] = (
                (_df_pos["var_contrib"] / _tot_var_pos * 100.0) if _tot_var_pos > 0 else 0.0
            )
            _df_pos = _df_pos.sort_values(by="frac_pos", ascending=True)

            _axes[0].barh(
                [f"ds{d}" for d in _df_pos["ds"]],
                _df_pos["frac_pos"],
                color="#ff7f0e",
                alpha=0.85,
                edgecolor="black",
            )
            _axes[0].set_xlabel("Fraction of Positive Uncertainty Variance (%)", fontsize=10)
            _axes[0].set_title(
                r"Top Contributors to $+\sigma_{\mathrm{syst}}$ ($\delta_{\mathrm{right}}$)",
                fontsize=11,
            )
            _axes[0].grid(axis="x", linestyle=":", alpha=0.6)
            for _i, _v in enumerate(_df_pos["frac_pos"]):
                _axes[0].text(_v + 0.5, _i, f"{_v:.1f}%", va="center", fontsize=8)
        else:
            _axes[0].text(
                0.5,
                0.5,
                "No positive shifts",
                ha="center",
                va="center",
                transform=_axes[0].transAxes,
            )
            _axes[0].axis("off")

        # 2. Negative shifts (contributing to delta_left)
        _df_neg = df_deltas[df_deltas["delta_left"] < 0].copy()
        if not _df_neg.empty:
            _df_neg["var_contrib"] = _df_neg["delta_left"] ** 2
            _tot_var_neg = _df_neg["var_contrib"].sum()
            _df_neg["frac_neg"] = (
                (_df_neg["var_contrib"] / _tot_var_neg * 100.0) if _tot_var_neg > 0 else 0.0
            )
            _df_neg = _df_neg.sort_values(by="frac_neg", ascending=True)

            _axes[1].barh(
                [f"ds{d}" for d in _df_neg["ds"]],
                _df_neg["frac_neg"],
                color="#1f77b4",
                alpha=0.85,
                edgecolor="black",
            )
            _axes[1].set_xlabel("Fraction of Negative Uncertainty Variance (%)", fontsize=10)
            _axes[1].set_title(
                r"Top Contributors to $-\sigma_{\mathrm{syst}}$ ($\delta_{\mathrm{left}}$)",
                fontsize=11,
            )
            _axes[1].grid(axis="x", linestyle=":", alpha=0.6)
            for _i, _v in enumerate(_df_neg["frac_neg"]):
                _axes[1].text(_v + 0.5, _i, f"{_v:.1f}%", va="center", fontsize=8)
        else:
            _axes[1].text(
                0.5,
                0.5,
                "No negative shifts",
                ha="center",
                va="center",
                transform=_axes[1].transAxes,
            )
            _axes[1].axis("off")

        plt.tight_layout()
        _pareto_ui = mo.vstack(
            [
                mo.md("#### Dataset Dominance & Variance Contribution Breakdown"),
                _fig,
            ]
        )

    _pareto_ui
    return


@app.cell
def _(
    Path,
    delta_summary,
    h5py,
    has_point_data,
    mo,
    nominal_mode,
    nominal_post,
    np,
    plt,
    point_selector,
    popt,
    study_dir_input,
):
    if not has_point_data:
        _overlay_ui = mo.md("")
    else:
        _sdir = Path(study_dir_input.value)
        _pid = point_selector.value
        _cat_tag = delta_summary.get("category", "lss_sscl")
        _h5_samples_path = _sdir / "tasks" / f"syst_{_cat_tag}_{_pid}_samples.h5"

        if not _h5_samples_path.is_file():
            _overlay_ui = mo.callout(
                mo.md(
                    f"ℹ️ Posterior samples file `{_h5_samples_path.name}` not found. "
                    "Task summary metrics are available in the table and barplots above."
                ),
                kind="info",
            )
        else:
            # Load stored MCMC posterior samples
            _dataset_samples = {}
            with h5py.File(_h5_samples_path, "r") as _hf:
                for _k in _hf.keys():
                    _dataset_samples[_k] = _hf[_k][:]

            _fig, _ax = plt.subplots(figsize=(12, 5))

            # Plot each dataset posterior
            _colors = plt.cm.tab20(np.linspace(0, 1, len(_dataset_samples)))
            for (_k, _samps), _col in zip(_dataset_samples.items(), _colors):
                _ax.hist(
                    _samps,
                    bins=40,
                    density=True,
                    histtype="step",
                    color=_col,
                    alpha=0.4,
                    linewidth=1.0,
                    label=_k if len(_dataset_samples) <= 15 else None,
                )

            # Plot nominal posterior if available
            if nominal_post is not None:
                _ax.hist(
                    nominal_post,
                    bins=60,
                    density=True,
                    histtype="stepfilled",
                    color="black",
                    alpha=0.25,
                    label="Nominal Posterior",
                )

            # Mark nominal mode and systematic bounds
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
                rf"Per-Dataset Posterior Distributions & Systematic Interval — {_pid}", fontsize=12
            )
            _ax.legend(fontsize=8, loc="upper right")
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
def _(
    Path,
    has_point_data,
    is_lss_dq,
    load_syst_result,
    mo,
    pd,
    point_selector,
    study_dir_input,
):
    if not has_point_data:
        _inspector_ui = mo.md("")
    else:
        _sdir = Path(study_dir_input.value)
        _pid = point_selector.value
        _cat = "lss_dq" if is_lss_dq else "lss_sscl"
        _h5_res = _sdir / f"syst_{_cat}_{_pid}.h5"

        if _h5_res.is_file():
            try:
                _loaded = load_syst_result(_h5_res)
                _meta = _loaded.get("meta", {})
                _per_post = _loaded.get("per_post", None)
                _metrics = _loaded.get("metrics", {})

                _details = f"""
                ### Systematic Result HDF5 (Format v2)
                **File:** `{_h5_res.name}`  
                - **Name:** `{_loaded.get('name')}` | **Format Version:** `{_meta.get('format_version', 'N/A')}`  
                - **Study / Point:** `{_meta.get('study_name', 'N/A')}` / `{_meta.get('point_id', 'N/A')}`  
                - **Evaluated Datasets:** `{len(_per_post) if _per_post is not None else 'N/A'}` variations  
                - **Summary Metrics:** `delta_right = {_metrics.get('delta_right', 'N/A')}`, `delta_left = {_metrics.get('delta_left', 'N/A')}`
                """

                _inspector_ui = mo.vstack(
                    [
                        mo.md(_details),
                        mo.md("#### Per-Variation Posterior Records (`per_post`)"),
                        mo.ui.table(_per_post)
                        if _per_post is not None
                        else mo.md("No per_post table."),
                    ]
                )
            except Exception as _e:
                _inspector_ui = mo.callout(
                    mo.md(f"⚠️ Error reading Format v2 HDF5 `{_h5_res.name}`: {_e}"), kind="warn"
                )
        else:
            _inspector_ui = mo.callout(
                mo.md(
                    f"ℹ️ Result artifact `{_h5_res.name}` not created yet. Run `--collect-syst` to generate Format v2 outputs."
                ),
                kind="info",
            )

    _inspector_ui
    return


if __name__ == "__main__":
    app.run()
