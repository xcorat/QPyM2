import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    from pathlib import Path
    from typing import Any, Dict, List, Optional, Tuple, Union

    import arviz as az
    import h5py
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import marimo as mo
    import numpy as np
    import pandas as pd
    from omegaconf import DictConfig

    from qpym2.fit.bias import apply_bias_correction
    from qpym2.plotting.syst_posterior import (
        plot_multi_posterior_overlay,
        plot_posterior_comparison,
    )
    from qpym2.storage.backend import read_fit
    from qpym2.storage.study import init_study_db
    from qpym2.studies._common import extract_logp_mode, load_bias_popt
    from qpym2.studies.assembly import assemble_study
    from qpym2.studies.syst_co_shift import (
        calculate_co_shift_weights,
        check_syst_co_shift,
        collect_syst_co_shift,
        fit_co_shift_parabola,
        plot_co_shift_parabola,
    )
    from qpym2.syst.core import compute_deltas
    from qpym2.syst.io import load_syst_result, save_syst_result
    from qpym2.utils.mathy import calc_intervals, get_mode, sample_weighted_mixture

    return (
        Any,
        Dict,
        DictConfig,
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
        calculate_co_shift_weights,
        check_syst_co_shift,
        collect_syst_co_shift,
        compute_deltas,
        extract_logp_mode,
        fit_co_shift_parabola,
        get_mode,
        h5py,
        init_study_db,
        json,
        load_bias_popt,
        load_syst_result,
        mo,
        np,
        pd,
        plot_co_shift_parabola,
        plot_multi_posterior_overlay,
        plot_posterior_comparison,
        plt,
        read_fit,
        sample_weighted_mixture,
        save_syst_result,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # QPyM2: $^{60}\text{Co}$ Energy Shift Marginalization Workbench

    Interactive inspection, diagnostics, and validation of the $^{60}\text{Co}$ Energy Calibration Shift systematic study:
    - **Physical Motivation**: Residual non-linearities and calibration discrepancies near the $0\nu\beta\beta$ ROI ($2527.5\text{ keV}$) are parameterized by an energy shift $\Delta_{\text{Co}}$ applied to the $^{60}\text{Co}$ multi-site background component.
    - **Log-Likelihood Parabola Fit**: The score $S_i = \text{mode}(\log \mathcal{L}^{(i)})$ is evaluated across candidate shifts $\Delta_{\text{Co}} \in [0.0, 0.8, 1.0, 1.2]\text{ keV}$ and modeled with a 2nd-degree polynomial:
      $$S(\Delta_{\text{Co}}) = a\Delta_{\text{Co}}^2 + b\Delta_{\text{Co}} + c \implies \Delta^* = -\frac{b}{2a}, \quad \sigma_{\Delta} = \frac{1}{\sqrt{-2a}}$$
    - **Softmax Posterior Weighting**: Candidate fits are weighted according to their relative likelihood / posterior evidence:
      $$w_i = \frac{\exp(S_i - \max(S))}{\sum_j \exp(S_j - \max(S))} \times \frac{\bar{N}}{N_i}$$
    - **Bias Correction & Mixture Marginalization**: Each shifted posterior is bias-corrected via $(N^{(i)} - b_i) / a_i$, and resampled according to weights $w_i$ to form the marginalized posterior $P_{\text{marg}}$.
    - **Systematic Deltas**: Evaluated against the nominal reference shift (default $\Delta_{\text{Co}} = 1.0\text{ keV}$):
      $$\Delta_{\text{mode}} = \text{mode}(P_{\text{marg}}) - \text{mode}(N_{\text{bc}}^{(\text{nom})})$$
      $$\Delta_{\text{right}} = \max\left(0, \text{HDI}_{\text{hi}}^{\text{marg}} - \text{HDI}_{\text{hi}}^{\text{nom}}\right), \quad \Delta_{\text{left}} = \min\left(0, \text{HDI}_{\text{lo}}^{\text{marg}} - \text{HDI}_{\text{lo}}^{\text{nom}}\right)$$
    """)
    return


@app.cell
def _(mo):
    study_dir_input = mo.ui.text(
        value="studies/example_co_shift",
        label="Study Directory",
        full_width=True,
    )
    bias_calib_input = mo.ui.text(
        value="studies/example_bias_study/bias_calibration.json",
        label="Bias Calibration JSON Path",
        full_width=True,
    )
    ref_shift_dropdown = mo.ui.dropdown(
        options=["1.0", "0.8", "1.2", "0.0"],
        value="1.0",
        label="Nominal Reference Shift (keV)",
    )
    estimator_dropdown = mo.ui.dropdown(
        options=["mode", "mean"],
        value="mode",
        label="Scoring & Shift Estimator",
    )
    lp_var_dropdown = mo.ui.dropdown(
        options=["pik", "lp"],
        value="pik",
        label="Likelihood Variable",
    )
    n_resamples_input = mo.ui.number(
        value=100000,
        start=10000,
        stop=500000,
        step=10000,
        label="Mixture Resample Draws",
    )

    mo.vstack(
        [
            mo.md("### Study Configuration & Controls"),
            study_dir_input,
            bias_calib_input,
            mo.hstack(
                [ref_shift_dropdown, estimator_dropdown, lp_var_dropdown, n_resamples_input],
                justify="start",
            ),
        ]
    )
    return (
        bias_calib_input,
        estimator_dropdown,
        lp_var_dropdown,
        n_resamples_input,
        ref_shift_dropdown,
        study_dir_input,
    )


@app.cell
def _(
    Path,
    assemble_study,
    json,
    mo,
    pd,
    study_dir_input,
):
    _sdir = Path(study_dir_input.value)
    study_obj = None
    status_df = pd.DataFrame()
    syst_summary_dict = {}

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

        _sum_path = _sdir / "syst_summary.json"
        if _sum_path.is_file():
            try:
                with open(_sum_path, "r", encoding="utf-8") as _f:
                    syst_summary_dict = json.load(_f)
            except Exception:
                syst_summary_dict = {}

        _scfg = getattr(study_obj.cfg, "syst", None) if study_obj else None
        _summary_md = f"""
        **Study Name:** `{study_obj.name if study_obj else 'N/A'}` | 🟣 **Category:** `co_shift`  
        **Study Directory:** `{_sdir}`  
        **Evaluation Points:** `{len(study_obj.points) if study_obj else 0}` | **Trace Status:** `{status_df['trace_exists'].sum() if not status_df.empty and 'trace_exists' in status_df.columns else 0}/{len(status_df)} Available`  
        **Saved Summary:** `{'Yes (syst_summary.json present)' if syst_summary_dict else 'Not yet saved'}`
        """

        _display_cols = [
            c
            for c in [
                "point_id",
                "co_shift",
                "fit_exists",
                "trace_exists",
                "fit_uuid",
                "status",
            ]
            if c in status_df.columns
        ]

        _ui_elements = [mo.md(_summary_md)]
        if not status_df.empty:
            _ui_elements.extend(
                [
                    mo.md("#### Co-Shift Evaluation Points & Cache Status"),
                    mo.ui.table(status_df[_display_cols] if not status_df.empty else status_df),
                ]
            )
        _header_ui = mo.vstack(_ui_elements)

    _header_ui
    return status_df, study_obj, syst_summary_dict


@app.cell
def _(
    Path,
    apply_bias_correction,
    az,
    bias_calib_input,
    calculate_co_shift_weights,
    compute_deltas,
    estimator_dropdown,
    extract_logp_mode,
    fit_co_shift_parabola,
    load_bias_popt,
    lp_var_dropdown,
    mo,
    n_resamples_input,
    np,
    pd,
    ref_shift_dropdown,
    sample_weighted_mixture,
    status_df,
    study_obj,
):
    shifts_list = []
    scores_list = []
    draw_counts_list = []
    raw_samples_dict = {}
    corrected_samples_dict = {}
    parabola_res = None
    weights_arr = None
    marginal_samples = None
    nominal_samples = None
    deltas = None
    analysis_error = None

    if study_obj is None or status_df.empty:
        analysis_error = "Study object or status DataFrame is unavailable."
    else:
        # Load bias calibration popt dict if available
        _calib_p = Path(bias_calib_input.value)
        _bias_popts = load_bias_popt(_calib_p) if _calib_p.is_file() else {}

        _ref_shift = float(ref_shift_dropdown.value)
        _estimator = estimator_dropdown.value
        _lp_var = lp_var_dropdown.value
        _n_resamples = int(n_resamples_input.value)

        for _, row in status_df.iterrows():
            if not row.get("trace_exists", False):
                continue
            _tr_p = Path(row["trace_path"])
            if not _tr_p.is_file():
                continue
            _shift = float(row["co_shift"])

            try:
                _idata = az.from_netcdf(str(_tr_p))
                _score = extract_logp_mode(_idata, var=_lp_var, estimator=_estimator)

                _post = _idata.posterior["ndbd"].values.flatten()
                _n_draws = len(_post)

                shifts_list.append(_shift)
                scores_list.append(_score)
                draw_counts_list.append(_n_draws)
                raw_samples_dict[_shift] = _post

                # Apply bias correction
                _popt = _bias_popts.get(_shift, [1.0, 0.0])
                _post_bc = apply_bias_correction(_post, slope=_popt[0], offset=_popt[1])
                corrected_samples_dict[_shift] = _post_bc
            except Exception as _ex:
                analysis_error = f"Error reading trace for shift {_shift}: {_ex}"
                break

        if not analysis_error and len(shifts_list) >= 3:
            try:
                # Sort by shift
                _sort_idx = np.argsort(shifts_list)
                shifts_arr = np.array(shifts_list)[_sort_idx]
                scores_arr = np.array(scores_list)[_sort_idx]
                draw_counts_arr = np.array(draw_counts_list)[_sort_idx]

                parabola_res = fit_co_shift_parabola(shifts_arr, scores_arr)
                weights_arr = calculate_co_shift_weights(scores_arr, draw_counts_arr)

                # Mixture sampling
                _sample_list = [corrected_samples_dict[s] for s in shifts_arr]
                marginal_samples = sample_weighted_mixture(
                    _sample_list, weights_arr, size=_n_resamples, random_state=42
                )

                # Nominal reference
                if _ref_shift in corrected_samples_dict:
                    nominal_samples = corrected_samples_dict[_ref_shift]
                else:
                    _closest_shift = shifts_arr[np.argmin(np.abs(shifts_arr - _ref_shift))]
                    nominal_samples = corrected_samples_dict[_closest_shift]

                deltas = compute_deltas(nominal_samples, marginal_samples, estimator=_estimator)
            except Exception as _ex2:
                analysis_error = f"Error during parabola / mixture analysis: {_ex2}"

    if analysis_error:
        _analysis_status_ui = mo.callout(mo.md(f"⚠️ {analysis_error}"), kind="alert")
    else:
        _analysis_status_ui = mo.callout(
            mo.md(
                f"✅ **Analysis Ready**: Loaded {len(shifts_list)} traces across shifts "
                f"`{sorted(shifts_list)} keV`. Reference shift: `{ref_shift_dropdown.value} keV`."
            ),
            kind="success",
        )

    _analysis_status_ui
    return (
        corrected_samples_dict,
        deltas,
        marginal_samples,
        nominal_samples,
        parabola_res,
        raw_samples_dict,
        scores_list,
        shifts_list,
        weights_arr,
    )


@app.cell
def _(
    deltas,
    estimator_dropdown,
    mo,
    parabola_res,
    plot_co_shift_parabola,
    plt,
    ref_shift_dropdown,
    shifts_list,
    scores_list,
    weights_arr,
):
    if parabola_res is None or weights_arr is None or deltas is None:
        _parabola_ui = mo.md("*Parabola and mixture metrics will appear once traces are loaded.*")
    else:
        _fig, _ax = plt.subplots(figsize=(8, 4.5))
        plot_co_shift_parabola(shifts_list, scores_list, parabola_res, ax=_ax)
        plt.tight_layout()

        _vx = parabola_res["vertex_shift"]
        _sig = parabola_res["sigma_shift"]
        _a, _b, _c = parabola_res["poly_coeffs"]

        _weights_table_rows = "".join(
            f"| `{s:.1f} keV` | `{sc:.2f}` | **`{w:.4f}`** ({w * 100:.1f}%) |\n"
            for s, sc, w in zip(sorted(shifts_list), scores_list, weights_arr)
        )

        _parabola_ui = mo.vstack(
            [
                mo.md("### 1. Log-Posterior Parabola & Softmax Weights"),
                mo.hstack(
                    [
                        mo.as_html(_fig),
                        mo.vstack(
                            [
                                mo.stat(
                                    label="Optimal Parabola Vertex Δ*",
                                    value=f"{_vx:.3f} keV",
                                    caption=f"1σ curvature width: ±{_sig:.3f} keV",
                                ),
                                mo.stat(
                                    label="Quadratic Curvature (a)",
                                    value=f"{_a:.3f}",
                                    caption="Expected a < 0 for local maximum",
                                ),
                                mo.md(
                                    f"| Shift | Log-Score | Softmax Weight |\n"
                                    f"| :--- | :--- | :--- |\n"
                                    f"{_weights_table_rows}"
                                ),
                            ]
                        ),
                    ],
                    justify="start",
                ),
            ]
        )

    _parabola_ui
    return


@app.cell
def _(
    calc_intervals,
    corrected_samples_dict,
    deltas,
    estimator_dropdown,
    get_mode,
    marginal_samples,
    mo,
    nominal_samples,
    np,
    plot_multi_posterior_overlay,
    plot_posterior_comparison,
    plt,
    ref_shift_dropdown,
    weights_arr,
):
    if marginal_samples is None or nominal_samples is None or deltas is None:
        _posterior_ui = mo.md("*Posterior distributions will appear once data is processed.*")
    else:
        # Multi-posterior overlay figure
        _sorted_shifts = sorted(corrected_samples_dict.keys())
        _samples_list = [corrected_samples_dict[s] for s in _sorted_shifts]
        _labels = [f"Δ = {s:.1f} keV (w = {w:.3f})" for s, w in zip(_sorted_shifts, weights_arr)]

        _fig_overlay, _ax_overlay = plt.subplots(figsize=(8, 4.5))
        plot_multi_posterior_overlay(
            _samples_list,
            labels=_labels,
            ax=_ax_overlay,
            title="Bias-Corrected Posteriors with Softmax Weights",
        )
        plt.tight_layout()

        # Comparison figure
        _fig_comp, _ax_comp = plt.subplots(figsize=(8, 4.5))
        plot_posterior_comparison(
            nominal_samples,
            marginal_samples,
            labels=(f"Nominal ({ref_shift_dropdown.value} keV)", "Co-Shift Marginalized"),
            deltas=deltas,
            ax=_ax_comp,
            title="Nominal vs Co-Shift Marginalized Posterior",
        )
        plt.tight_layout()

        _d_mode = deltas.delta_mode
        _d_left = deltas.delta_left
        _d_right = deltas.delta_right

        _nom_mode = get_mode(nominal_samples)
        _marg_mode = get_mode(marginal_samples)
        _nom_hdi = calc_intervals(nominal_samples, p_low=5.0, p_high=95.0)
        _marg_hdi = calc_intervals(marginal_samples, p_low=5.0, p_high=95.0)

        _posterior_ui = mo.vstack(
            [
                mo.md("### 2. Posterior Overlays & Marginalization Diagnostics"),
                mo.hstack([mo.as_html(_fig_overlay), mo.as_html(_fig_comp)], justify="start"),
                mo.md("#### Systematic Delta Shifts (counts)"),
                mo.hstack(
                    [
                        mo.stat(
                            label="Mode Shift (Δ_mode)",
                            value=f"{_d_mode:+.2f} counts",
                            caption=f"Nom: {_nom_mode:.1f} → Marg: {_marg_mode:.1f}",
                        ),
                        mo.stat(
                            label="Negative Shift (Δ_left)",
                            value=f"{_d_left:.2f} counts",
                            caption="Lower HDI edge difference (<= 0)",
                        ),
                        mo.stat(
                            label="Positive Shift (Δ_right)",
                            value=f"{_d_right:+.2f} counts",
                            caption="Upper HDI edge difference (>= 0)",
                        ),
                    ],
                    justify="start",
                ),
            ]
        )

    _posterior_ui
    return


@app.cell
def _(
    bias_calib_input,
    collect_syst_co_shift,
    estimator_dropdown,
    lp_var_dropdown,
    mo,
    n_resamples_input,
    ref_shift_dropdown,
    study_obj,
):
    save_run_btn = mo.ui.run_button(label="💾 Run collect_syst_co_shift & Save Artifacts to Disk")

    _save_feedback = None
    if save_run_btn.value:
        if study_obj is None:
            _save_feedback = mo.callout(
                mo.md("⚠️ Cannot save: Study object is not assembled."), kind="alert"
            )
        else:
            try:
                _res_summary = collect_syst_co_shift(
                    study_obj,
                    reference_shift=float(ref_shift_dropdown.value),
                    bias_calibration=bias_calib_input.value,
                    estimator=estimator_dropdown.value,
                    lp_var=lp_var_dropdown.value,
                    n_resamples=int(n_resamples_input.value),
                )
                _sdir = study_obj.dir
                _save_feedback = mo.callout(
                    mo.md(
                        f"🎉 **Successfully Collected & Saved Systematic Artifacts!**\n\n"
                        f"- Summary JSON: `{_sdir / 'syst_summary.json'}`\n"
                        f"- Marginal Posterior HDF5: `{_sdir / 'syst_co_shift_marginal.h5'}`\n"
                        f"- Parabola Plot: `{_sdir / 'plots/co_shift_logp_parabola.png'}`\n"
                        f"- Multi-Posterior Overlay: `{_sdir / 'plots/co_shift_posterior_overlay.png'}`\n"
                        f"- Comparison Plot: `{_sdir / 'plots/co_shift_marginal_comparison.png'}`\n\n"
                        f"**Result Deltas**: `Δ_mode = {_res_summary.get('delta_mode', 0.0):.2f}`, "
                        f"`Δ_left = {_res_summary.get('delta_left', 0.0):.2f}`, "
                        f"`Δ_right = {_res_summary.get('delta_right', 0.0):.2f}`"
                    ),
                    kind="success",
                )
            except Exception as _ex:
                _save_feedback = mo.callout(
                    mo.md(f"⚠️ Error executing `collect_syst_co_shift`: {_ex}"), kind="alert"
                )

    mo.vstack(
        [
            mo.md("### 3. Disk Persistence & Downstream Export"),
            mo.md(
                "Click below to execute `collect_syst_co_shift` with the configured parameters, "
                "writing `syst_summary.json` (for `marginalize.py`), Format v2 HDF5, and publication plots."
            ),
            save_run_btn,
            _save_feedback if _save_feedback is not None else mo.md(""),
        ]
    )
    return (save_run_btn,)


if __name__ == "__main__":
    app.run()
