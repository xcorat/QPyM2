import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    from pathlib import Path

    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    import pandas as pd

    from qpym2.fit.bias import apply_bias_correction
    from qpym2.studies.assembly import assemble_study

    return (
        Path,
        apply_bias_correction,
        assemble_study,
        json,
        mo,
        np,
        pd,
        plt,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # QPyM2: Expected Sensitivity Study Workbench

    Interactive evaluation and validation of Bayesian expected sensitivity ($0\nu\beta\beta$) studies in QPyM2 v0.6:
    - **Background-Only Hypothesis ($H_0$)**: Pseudo-experiments generated with zero injected signal ($N_{\text{inj}} = 0$).
    - **90% C.L. Sensitivity ($N_{90}$)**: Distribution of Bayesian 90% credibility upper limits (`ndbd_p90`).
    - **Expected Sensitivity & Containment**: Median sensitivity with $\pm 1\sigma$ ($[16\%, 84\%]$) and $\pm 2\sigma$ ($[2.5\%, 97.5\%]$) containment bands.
    - **Systematic Scans**: Sensitivity evolution as a function of systematic parameters (e.g., Co-60 energy shifts).
    - **Bias Correction**: Integration with linear calibration ($N_{\text{bc}} = (N_{90} - b)/a$) from bias studies.
    - **Data Comparison**: Overlay observed data limit and evaluate $P(T_{1/2}^{\text{toy}} < T_{1/2}^{\text{obs}})$.
    """)
    return


@app.cell
def _(mo):
    study_dir_input = mo.ui.text(
        value="studies/example_sens_study",
        label="Sensitivity Study Directory",
        full_width=True,
    )
    bias_calib_input = mo.ui.text(
        value="studies/example_bias_study/bias_calibration.json",
        label="Optional Bias Calibration JSON",
        full_width=True,
    )
    data_study_input = mo.ui.text(
        value="studies/example_co_shift",
        label="Optional Data Fit Study (for real observed limits)",
        full_width=True,
    )

    mo.vstack(
        [
            mo.md("### Study Configuration & Paths"),
            study_dir_input,
            bias_calib_input,
            data_study_input,
        ]
    )
    return bias_calib_input, data_study_input, study_dir_input


@app.cell
def _(
    Path,
    assemble_study,
    bias_calib_input,
    data_study_input,
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
        sens_summary_dict = {}
        bias_calib_dict = {}
        data_results_df = pd.DataFrame()
        _header_ui = mo.callout(mo.md(f"⚠️ Directory `{_sdir}` does not exist."), kind="alert")
    else:
        try:
            study_obj = assemble_study(_sdir)
            status_df = study_obj.check_sens_toys()
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

        # Load sens_summary.json
        _sum_path = _sdir / "sens_summary.json"
        if _sum_path.is_file():
            with open(_sum_path, "r", encoding="utf-8") as _f:
                sens_summary_dict = json.load(_f)
        else:
            sens_summary_dict = {}

        # Load bias calibration if given
        _bpath = Path(bias_calib_input.value)
        if not _bpath.is_absolute() and not _bpath.exists():
            _bpath = _sdir / _bpath
        if _bpath.is_file():
            with open(_bpath, "r", encoding="utf-8") as _f:
                bias_calib_dict = json.load(_f)
        else:
            bias_calib_dict = {}

        # Load data fit study if available
        _dpath = Path(data_study_input.value)
        _dres = _dpath / "results.csv"
        if _dres.is_file():
            data_results_df = pd.read_csv(_dres)
        else:
            data_results_df = pd.DataFrame()

        if study_obj is not None:
            _n_points = len(study_obj.points)
            _n_toys_total = len(results_df) if not results_df.empty else 0
            _summary_md = f"""
            **Study Name:** `{study_obj.name}`  
            **Evaluation Points:** `{_n_points}` | **Collected Sensitivity Toys:** `{_n_toys_total}`  
            **Reference Coordinates:** `co_cut={getattr(study_obj.cfg.reference, "co_cut", "N/A")}`, `ndbd_cutoff={getattr(study_obj.cfg.reference, "ndbd_cutoff", "N/A")}`  
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
    return (
        bias_calib_dict,
        data_results_df,
        results_df,
        sens_summary_dict,
        status_df,
        study_obj,
    )


@app.cell
def _(mo, results_df):
    _available_points = (
        ["All Points"] + sorted(list(results_df["point_id"].unique()))
        if (not results_df.empty and "point_id" in results_df.columns)
        else ["All Points"]
    )

    point_dropdown = mo.ui.dropdown(
        options=_available_points,
        value=_available_points[0],
        label="Select Point for 1D Detail",
    )

    apply_bias_cb = mo.ui.checkbox(
        value=True,
        label="Apply Bias Correction (using bias_calibration.json)",
    )

    nbins_slider = mo.ui.slider(
        start=10,
        stop=60,
        step=5,
        value=30,
        label="Histogram Bins",
    )

    mo.vstack(
        [
            mo.md("### Visualization Controls"),
            mo.hstack([point_dropdown, apply_bias_cb, nbins_slider]),
        ]
    )
    return apply_bias_cb, nbins_slider, point_dropdown


@app.cell
def _(
    apply_bias_cb,
    bias_calib_dict,
    data_results_df,
    mo,
    nbins_slider,
    np,
    plt,
    point_dropdown,
    results_df,
):
    if results_df.empty or "ndbd_p90" not in results_df.columns:
        _dist_plot = mo.md("*No sensitivity results available to plot. Run `--run-sens` first.*")
    else:
        _sel_point = point_dropdown.value
        _sub_df = (
            results_df
            if _sel_point == "All Points"
            else results_df[results_df["point_id"] == _sel_point].copy()
        )

        _p90_raw = _sub_df["ndbd_p90"].dropna().values
        _mode_raw = _sub_df["ndbd_mode"].dropna().values if "ndbd_mode" in _sub_df.columns else np.array([])

        _fig, _axes = plt.subplots(1, 2, figsize=(11, 4.2))

        # Check bias calibration for selected point
        _calib = None
        if apply_bias_cb.value and bias_calib_dict:
            if _sel_point in bias_calib_dict:
                _calib = bias_calib_dict[_sel_point]
            elif len(_sub_df["co_shift"].unique()) == 1:
                _cs = float(_sub_df["co_shift"].iloc[0])
                for _k, _v in bias_calib_dict.items():
                    if isinstance(_v, dict) and _v.get("co_shift") == _cs:
                        _calib = _v
                        break

        if _calib and "a" in _calib and "b" in _calib:
            _a, _b = float(_calib["a"]), float(_calib["b"])
            _p90_vals = np.clip((_p90_raw - _b) / _a, 0.0, None)
            _p90_title_suffix = f" (Bias Corrected: a={_a:.3f}, b={_b:.2f})"
        else:
            _p90_vals = _p90_raw
            _p90_title_suffix = " (Raw)"

        # Axis 0: 90% C.L. upper limits
        _ax0 = _axes[0]
        _med = float(np.median(_p90_vals))
        _p16 = float(np.percentile(_p90_vals, 16))
        _p84 = float(np.percentile(_p90_vals, 84))

        _ax0.hist(
            _p90_vals,
            bins=int(nbins_slider.value),
            color="#2a75d3",
            alpha=0.65,
            edgecolor="white",
            density=True,
            label=f"Sensitivity Toys (N={len(_p90_vals)})",
        )
        _ax0.axvline(_med, color="#124182", linestyle="--", linewidth=2.0, label=f"Median = {_med:.2f}")
        _ax0.axvspan(_p16, _p84, color="#2a75d3", alpha=0.2, label=f"±1σ: [{_p16:.2f}, {_p84:.2f}]")

        # Overlay observed data limit if available
        if not data_results_df.empty:
            _dsub = (
                data_results_df
                if _sel_point == "All Points"
                else data_results_df[data_results_df["point_id"] == _sel_point]
            )
            if not _dsub.empty and "ndbd_q90" in _dsub.columns:
                _dval = float(_dsub["ndbd_q90"].iloc[0])
                _ax0.axvline(_dval, color="#d9381e", linewidth=2.2, label=f"Observed = {_dval:.2f}")

        _ax0.set_xlabel("ndbd 90% CL Upper Limit (counts)")
        _ax0.set_ylabel("Probability Density")
        _ax0.set_title(f"Sensitivity Distribution{_p90_title_suffix}")
        _ax0.legend(fontsize=8)

        # Axis 1: Posterior mode
        _ax1 = _axes[1]
        if len(_mode_raw) > 0:
            _mode_med = float(np.median(_mode_raw))
            _ax1.hist(
                _mode_raw,
                bins=int(nbins_slider.value),
                color="#2ca02c",
                alpha=0.65,
                edgecolor="white",
                density=True,
                label=f"Signal Modes (N={len(_mode_raw)})",
            )
            _ax1.axvline(_mode_med, color="#136013", linestyle="--", linewidth=2.0, label=f"Median Mode = {_mode_med:.2f}")
            _ax1.set_xlabel("ndbd Posterior Mode (counts)")
            _ax1.set_ylabel("Probability Density")
            _ax1.set_title("Fitted Signal Mode under H0")
            _ax1.legend(fontsize=8)

        plt.tight_layout()
        _dist_plot = mo.as_html(_fig)

    _dist_plot
    return


@app.cell
def _(apply_bias_cb, bias_calib_dict, data_results_df, mo, np, pd, plt, results_df):
    if results_df.empty or "co_shift" not in results_df.columns:
        _scan_ui = mo.md("*No sweep results to plot.*")
    else:
        _scan_rows = []
        for _cs, _grp in results_df.groupby("co_shift"):
            _p90 = _grp["ndbd_p90"].dropna().values
            if len(_p90) == 0:
                continue

            _pid = _grp["point_id"].iloc[0] if "point_id" in _grp.columns else ""
            _calib = None
            if apply_bias_cb.value and bias_calib_dict:
                if _pid in bias_calib_dict:
                    _calib = bias_calib_dict[_pid]
                else:
                    for _k, _v in bias_calib_dict.items():
                        if isinstance(_v, dict) and _v.get("co_shift") == _cs:
                            _calib = _v
                            break

            if _calib and "a" in _calib and "b" in _calib:
                _p90_use = np.clip((_p90 - float(_calib["b"])) / float(_calib["a"]), 0.0, None)
            else:
                _p90_use = _p90

            _scan_rows.append(
                {
                    "co_shift": float(_cs),
                    "median": float(np.median(_p90_use)),
                    "p16": float(np.percentile(_p90_use, 16)),
                    "p84": float(np.percentile(_p90_use, 84)),
                    "p025": float(np.percentile(_p90_use, 2.5)),
                    "p975": float(np.percentile(_p90_use, 97.5)),
                    "n_toys": len(_p90_use),
                }
            )

        _scan_df = pd.DataFrame(_scan_rows).sort_values("co_shift")

        _fig2, _ax2 = plt.subplots(figsize=(7.5, 4.5))

        # ±2σ band
        _ax2.fill_between(
            _scan_df["co_shift"],
            _scan_df["p025"],
            _scan_df["p975"],
            color="#a8c7fa",
            alpha=0.45,
            label="±2σ Containment (95%)",
        )
        # ±1σ band
        _ax2.fill_between(
            _scan_df["co_shift"],
            _scan_df["p16"],
            _scan_df["p84"],
            color="#2a75d3",
            alpha=0.35,
            label="±1σ Containment (68%)",
        )
        # Median
        _ax2.plot(
            _scan_df["co_shift"],
            _scan_df["median"],
            marker="o",
            color="#124182",
            linewidth=2.0,
            label="Median Expected 90% CL Limit",
        )

        # Overlay observed data limit if present
        if not data_results_df.empty and "co_shift" in data_results_df.columns and "ndbd_q90" in data_results_df.columns:
            _d_sorted = data_results_df.sort_values("co_shift")
            _ax2.plot(
                _d_sorted["co_shift"],
                _d_sorted["ndbd_q90"],
                marker="s",
                color="#d9381e",
                linewidth=1.8,
                linestyle="--",
                label="Observed Data Limit (ndbd_q90)",
            )

        _ax2.set_xlabel("Co-60 Energy Shift (keV)", fontsize=11)
        _ax2.set_ylabel("Signal Count (counts)", fontsize=11)
        _ax2.set_title("Expected 90% CL Sensitivity vs Co-60 Shift", fontsize=12)
        _ax2.legend(fontsize=9)
        _ax2.grid(True, alpha=0.3)
        plt.tight_layout()

        _scan_ui = mo.vstack(
            [
                mo.md("### Systematic Sensitivity Scan"),
                mo.as_html(_fig2),
                mo.md("#### Summary Statistics Table"),
                mo.ui.table(_scan_df),
            ]
        )

    _scan_ui
    return


if __name__ == "__main__":
    app.run()
