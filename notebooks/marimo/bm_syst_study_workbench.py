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
    from omegaconf import OmegaConf

    from qpym2.fit.bias import apply_bias_correction
    from qpym2.io.bkg_model import list_bm_variants, read_bm_variant_modes
    from qpym2.storage.backend import read_fit
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

    Interactive inspection, diagnostics, and marginalisation for Background Model (BM) systematics in QPyM2 v0.6:
    - **BM Threshold (`bm_threshold`)**: Evaluates the impact of changing the detector threshold from nominal 200 keV to 250 keV (`Threshold250` vs nominal `Threshold200`). This modifies the M2 fit ROI to $e_1, e_2 \ge 250\text{ keV}$ and overrides MC normalisations with the upstream JAGS modes.
    - **BM Binning (`bm_binning`)**: Evaluates upstream JAGS background model binning choices (`MinBinSize_5_20`, `MinBinSize_25_60`, `BinningIntraAlphaBands`) on the nominal geometry. The systematic shifts are aggregated via an **envelope** ($\Delta_{\text{left}} = \min_v \Delta_{\text{left}}^{(v)}$, $\Delta_{\text{right}} = \max_v \Delta_{\text{right}}^{(v)}$).
    - **Master Rate Marginalisation (`syst_rate`)**: Combines systematic shifts in quadrature into split-Gaussian widths ($\sigma_{\text{left}}, \sigma_{\text{right}}$), smears the bias-corrected nominal posterior, convolves with exposure and PSA efficiency, and extracts the 90% C.L. half-life limit ($T_{1/2}^{0\nu}$).
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
    results_raw_input = mo.ui.text(
        value="/data/cuore/Analysis/BMMC/inputs/ResultsRaw.h5",
        label="ResultsRaw.h5 Path (for upstream BM modes)",
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
            study_type_select,
            study_dir_input,
            results_raw_input,
            bias_calib_input,
        ]
    )
    return bias_calib_input, results_raw_input, study_dir_input, study_type_select


@app.cell
def _(
    Path,
    assemble_study,
    json,
    mo,
    pd,
    study_dir_input,
    study_type_select,
):
    # Prefer explicit study_dir_input unless user chose dropdown preset
    _chosen_dir = study_dir_input.value if study_dir_input.value else study_type_select.value
    _sdir = Path(_chosen_dir)

    study_obj = None
    status_df = pd.DataFrame()
    results_df = pd.DataFrame()
    syst_summary_dict = {}
    delta_table_df = pd.DataFrame()
    rate_summary_dict = {}

    if not _sdir.exists():
        _header_ui = mo.callout(mo.md(f"⚠️ Study directory `{_sdir}` does not exist."), kind="alert")
    else:
        try:
            study_obj = assemble_study(_sdir)
            status_df = study_obj.check_syst()
        except Exception as _e:
            _header_ui = mo.callout(mo.md(f"⚠️ Error assembling study or checking status: {_e}"), kind="alert")

        _res_path = _sdir / "results.csv"
        results_df = pd.read_csv(_res_path) if _res_path.is_file() else pd.DataFrame()

        _sum_path = _sdir / "syst_summary.json"
        if _sum_path.is_file():
            with open(_sum_path, "r", encoding="utf-8") as _f:
                syst_summary_dict = json.load(_f)

        _dt_path = _sdir / "delta_table.csv"
        delta_table_df = pd.read_csv(_dt_path) if _dt_path.is_file() else pd.DataFrame()

        _rate_sum_path = _sdir / "syst_rate_summary.json"
        if _rate_sum_path.is_file():
            with open(_rate_sum_path, "r", encoding="utf-8") as _f:
                rate_summary_dict = json.load(_f)

        _scfg = getattr(study_obj.cfg, "syst", None) if study_obj else None
        _category = getattr(_scfg, "category", getattr(study_obj.cfg.study, "type", "syst")) if study_obj else "unknown"

        _summary_md = f"""
        **Study Name:** `{study_obj.name if study_obj else 'N/A'}` | **Category / Type:** `{_category}`  
        **Study Directory:** `{_sdir}`  
        **Points Count:** `{len(study_obj.points) if study_obj else 0}` | **Collected Rows:** `{len(results_df)}`
        """

        _ui_elements = [mo.md(_summary_md)]
        if not status_df.empty:
            _ui_elements.extend([mo.md("#### Task Status"), mo.ui.table(status_df)])
        _header_ui = mo.vstack(_ui_elements)

    _header_ui
    return (
        delta_table_df,
        rate_summary_dict,
        results_df,
        status_df,
        study_obj,
        syst_summary_dict,
    )


@app.cell
def _(Path, list_bm_variants, mo, results_raw_input):
    _raw_path = Path(results_raw_input.value)
    available_variants = []
    if _raw_path.is_file():
        try:
            available_variants = list_bm_variants(_raw_path)
        except Exception:
            available_variants = ["Threshold250", "MinBinSize_5_20", "MinBinSize_25_60", "BinningIntraAlphaBands"]
    else:
        available_variants = ["Threshold250", "MinBinSize_5_20", "MinBinSize_25_60", "BinningIntraAlphaBands"]

    variant_selector = mo.ui.dropdown(
        options=available_variants,
        value="Threshold250" if "Threshold250" in available_variants else (available_variants[0] if available_variants else ""),
        label="Select BM Variant to Inspect Upstream Normalisations",
    )

    mo.vstack([
        mo.md("### Upstream Background Model Modes (`ResultsRaw.h5`)"),
        variant_selector,
    ])
    return available_variants, variant_selector


@app.cell
def _(
    Path,
    mo,
    pd,
    plt,
    read_bm_variant_modes,
    results_raw_input,
    variant_selector,
):
    _raw_path = Path(results_raw_input.value)
    _sel_variant = variant_selector.value

    if not _raw_path.is_file() or not _sel_variant:
        _bm_modes_ui = mo.md("_ResultsRaw.h5 not found or no variant selected._")
    else:
        try:
            modes_series = read_bm_variant_modes(_raw_path, _sel_variant)
            df_modes = pd.DataFrame({"Component": modes_series.index, "Mode": modes_series.values}).sort_values("Mode", ascending=False)

            # Bar plot of top 12 background components
            top_df = df_modes.head(12)
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.barh(top_df["Component"][::-1], top_df["Mode"][::-1], color="#3366cc", alpha=0.85)
            ax.set_xlabel("Fitted Component Mode (Counts)")
            ax.set_title(f"Upstream Mode Normalisations: {_sel_variant} (Top 12 Components)")
            ax.grid(axis="x", linestyle="--", alpha=0.5)
            plt.tight_layout()

            _bm_modes_ui = mo.vstack([
                mo.ui.table(df_modes.head(20)),
                mo.as_html(fig),
            ])
        except Exception as _e:
            _bm_modes_ui = mo.callout(mo.md(f"⚠️ Could not read modes for variant `{_sel_variant}`: {_e}"), kind="alert")

    _bm_modes_ui
    return ax, df_modes, fig, modes_series, top_df


@app.cell
def _(mo, syst_summary_dict):
    if not syst_summary_dict:
        _syst_summary_ui = mo.md("_No `syst_summary.json` found in the selected study directory._")
    else:
        rows = []
        for pid, sdata in syst_summary_dict.items():
            r = {
                "point_id": pid,
                "category": sdata.get("category", ""),
                "delta_mode": sdata.get("delta_mode", 0.0),
                "delta_left": sdata.get("delta_left", 0.0),
                "delta_right": sdata.get("delta_right", 0.0),
                "mode_nom": sdata.get("mode_nom", np.nan),
            }
            if "per_variant" in sdata:
                for vname, vvals in sdata["per_variant"].items():
                    r[f"{vname}_dl"] = vvals.get("delta_left", 0.0)
                    r[f"{vname}_dr"] = vvals.get("delta_right", 0.0)
            rows.append(r)

        df_sum = pd.DataFrame(rows)
        _syst_summary_ui = mo.vstack([
            mo.md("### Systematic Shifts Summary (`syst_summary.json`)"),
            mo.ui.table(df_sum),
        ])

    _syst_summary_ui
    return df_sum, pid, r, rows, sdata


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

    mo.vstack([
        exposure_slider,
        eff_mc_slider,
        psa_mean_slider,
        psa_sig_slider,
        n_smear_slider,
    ])
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
    # 1. Resolve systematic deltas table
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
        # Canonical baseline delta table
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
            SystDelta(name=sname, delta_right=dr, delta_left=dl)
            for sname, dl, dr in default_deltas
        ]

    # Combine in quadrature
    sigma_left, sigma_right = combine_deltas([(d.delta_left, d.delta_right) for d in deltas_list])

    # Mock or real nominal posterior centered at 14.6 counts
    rng = np.random.default_rng(42)
    nominal_bc = np.clip(rng.normal(loc=14.64, scale=4.8, size=25000), 0.0, None)

    # Perform split-Gaussian smearing
    smeared_post = smear_posterior(
        nominal_bc,
        sigma_right=float(abs(sigma_right)),
        sigma_left=float(abs(sigma_left)),
        n_smear=int(n_smear_slider.value),
        rng_seed=42,
    )

    # Convert to rate
    rate_post = to_rate(
        smeared_post,
        eff_mc=float(eff_mc_slider.value),
        exposure=float(exposure_slider.value),
        psa_mean=float(psa_mean_slider.value),
        psa_sig=float(psa_sig_slider.value),
        rng_seed=42,
        clip=True,
    )

    # Calculate statistics
    mode_nom = float(get_mode(nominal_bc))
    mode_smeared = float(get_mode(smeared_post))
    rate_mode = float(get_mode(rate_post))
    rate_90 = float(np.percentile(rate_post, 90))
    half_life_90 = float(np.log(2.0) / rate_90) if rate_90 > 0 else np.inf

    # Generate plots
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
            "hdi": (float(calc_intervals([rate_post])[0][1]), float(calc_intervals([rate_post])[0][2])),
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

    _marginalisation_ui = mo.vstack([
        mo.md(results_banner),
        mo.md("#### 1. Systematics Breakdown & Quadrature Sum"),
        mo.as_html(fig_dt),
        mo.md("#### 2. Posterior Smearing"),
        mo.as_html(fig_comb),
        mo.md("#### 3. Marginalised Decay Rate Posterior"),
        mo.as_html(fig_rate),
    ])

    _marginalisation_ui
    return (
        default_deltas,
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


if __name__ == "__main__":
    app.run()
