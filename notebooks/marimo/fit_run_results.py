import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import shlex
    import subprocess
    import sys
    from pathlib import Path

    import arviz as az
    import marimo as mo
    import numpy as np
    import pandas as pd

    from qpym2.cfg import get_config, get_comp_cfg, get_h2uv_model
    from qpym2.cfg.config_builders import (
        build_input_config,
        get_artifact_store_dir,
        get_data_names,
        get_fit_dir,
        get_jags_path,
        get_plots_dir,
    )
    from qpym2.cfg.hydra_adapter import HydraConfigManager, resolve_config_dir
    from qpym2.fit.params import model_cfg as build_model_cfg
    from qpym2.io.fit import read_fit as read_fit_artifact, read_input
    from qpym2.plotting.projection_plots import plot_proj_h2
    from qpym2.storage.hasher import derive_id
    from qpym2.storage.schemas import FitParams, InputParams
    from qpym2.storage.store import ArtifactStore
    from qpym2.utils.hists import get_hist_settings, get_sum_hist
    from qpym2.utils.mathy import find_mode, get_mode

    return (
        ArtifactStore,
        FitParams,
        HydraConfigManager,
        InputParams,
        Path,
        az,
        build_input_config,
        build_model_cfg,
        derive_id,
        find_mode,
        get_artifact_store_dir,
        get_comp_cfg,
        get_config,
        get_data_names,
        get_fit_dir,
        get_h2uv_model,
        get_hist_settings,
        get_jags_path,
        get_mode,
        get_plots_dir,
        get_sum_hist,
        mo,
        np,
        os,
        pd,
        plot_proj_h2,
        read_fit_artifact,
        read_input,
        resolve_config_dir,
        shlex,
        subprocess,
        sys,
    )


@app.cell
def _(mo, resolve_config_dir):
    default_config_dir, config_dir_source = resolve_config_dir()
    config_dir_input = mo.ui.text(
        value=str(default_config_dir),
        label="Config directory",
        full_width=True,
    )
    mo.vstack(
        [
            mo.md(
                """
                # Run one QPyM2 fit and inspect its results

                Prepare or reuse the fit-input artifact, run one configured fit,
                then inspect posterior summaries, convergence, and data-versus-fit
                projections. The fit is executed by `qpym2.run.run_single_fit` and
                published to the configured artifact store.
                """
            ),
            config_dir_input,
            mo.md(f"Config directory source: `{config_dir_source}`."),
        ]
    )
    return (config_dir_input,)


@app.cell
def _(HydraConfigManager, config_dir_input, mo):
    HydraConfigManager.set_config_dir(config_dir_input.value)
    config_options = HydraConfigManager.list_options()

    environment_options = config_options.get("environments", ["local"])
    environment = mo.ui.dropdown(
        options=environment_options,
        value="local" if "local" in environment_options else environment_options[0],
        label="Environment",
    )
    analysis_options = config_options.get("analysis", ["base"])
    analysis = mo.ui.dropdown(
        options=analysis_options,
        value="tech_note" if "tech_note" in analysis_options else analysis_options[0],
        label="Analysis",
    )
    model_config_options = config_options.get("models", ["tech_note"])
    model_config = mo.ui.dropdown(
        options=model_config_options,
        value="tech_note" if "tech_note" in model_config_options else model_config_options[0],
        label="Model config set",
    )
    mo.hstack([environment, analysis, model_config])
    return analysis, environment, model_config


@app.cell
def _(analysis, config_dir_input, environment, get_config, mo, model_config):
    cfg = get_config(
        environment=environment.value,
        analysis=analysis.value,
        models=model_config.value,
        config_dir=config_dir_input.value,
    )
    model_key_options = list(cfg.models.keys())
    model_key = mo.ui.dropdown(
        options=model_key_options,
        value="tech_note_200" if "tech_note_200" in model_key_options else model_key_options[0],
        label="Model entry",
    )
    mo.hstack([mo.md("Model entries in this config:"), model_key])
    return cfg, model_key


@app.cell
def _(cfg, mo):
    co_cut_options = ["None", *list(cfg.cuts.co_cuts.keys())]
    ndbd_cutoff_options = ["None", *list(cfg.cuts.ndbd_cutoff.keys())]
    co_cut = mo.ui.dropdown(
        options=co_cut_options,
        value="gamma10_v150" if "gamma10_v150" in co_cut_options else co_cut_options[0],
        label="Co-60 cut",
    )
    ndbd_cutoff = mo.ui.dropdown(
        options=ndbd_cutoff_options,
        value="u1740" if "u1740" in ndbd_cutoff_options else ndbd_cutoff_options[0],
        label="0νββ cutoff",
    )
    co_shift = mo.ui.dropdown(
        options=[0.0, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5],
        value=0.8,
        label="Co-60 shift [keV]",
    )
    ndbd_shift = mo.ui.dropdown(
        options=[0.0, 0.3, 0.5, 0.8, 1.0],
        value=0.0,
        label="0νββ shift [keV]",
    )
    blinded = mo.ui.checkbox(value=False, label="Use blinded data")
    skip_roi = mo.ui.checkbox(value=False, label="Exclude an ROI interval")
    roi_min = mo.ui.number(value=2515.0, label="Excluded ROI min [keV]")
    roi_max = mo.ui.number(value=2540.0, label="Excluded ROI max [keV]")
    prior_width = mo.ui.dropdown(
        options=["5sig", "10sig", "syst", "2syst"],
        value="5sig",
        label="Prior width",
    )
    chains = mo.ui.number(value=8, start=1, step=1, label="Chains")
    draws = mo.ui.number(value=20000, start=1, step=1000, label="Draws per chain")
    single_core = mo.ui.checkbox(value=False, label="Run sampler single-core")
    force = mo.ui.checkbox(value=False, label="Force a new fit (ignore cached fit)")
    save_plots = mo.ui.checkbox(value=False, label="Save displayed figures")
    mo.vstack(
        [
            mo.md("## Fit point and sampler settings"),
            mo.hstack([co_cut, ndbd_cutoff]),
            mo.hstack([co_shift, ndbd_shift]),
            mo.hstack([blinded, skip_roi]),
            mo.hstack([roi_min, roi_max]),
            mo.hstack([prior_width, chains, draws]),
            mo.hstack([single_core, force, save_plots]),
        ]
    )
    return (
        blinded,
        chains,
        co_cut,
        co_shift,
        draws,
        force,
        ndbd_cutoff,
        ndbd_shift,
        prior_width,
        roi_max,
        roi_min,
        save_plots,
        single_core,
        skip_roi,
    )


@app.cell
def _(
    ArtifactStore,
    FitParams,
    InputParams,
    Path,
    analysis,
    blinded,
    build_input_config,
    build_model_cfg,
    cfg,
    chains,
    co_cut,
    co_shift,
    derive_id,
    draws,
    environment,
    force,
    get_artifact_store_dir,
    get_comp_cfg,
    get_data_names,
    get_fit_dir,
    get_h2uv_model,
    get_jags_path,
    get_plots_dir,
    mo,
    model_config,
    model_key,
    ndbd_cutoff,
    ndbd_shift,
    pd,
    prior_width,
    roi_max,
    roi_min,
    shlex,
    single_core,
    skip_roi,
    sys,
):
    selected_co_cut = None if co_cut.value == "None" else co_cut.value
    selected_ndbd_cutoff = None if ndbd_cutoff.value == "None" else ndbd_cutoff.value
    selected_skip_roi = (float(roi_min.value), float(roi_max.value)) if skip_roi.value else None
    if selected_skip_roi is not None and selected_skip_roi[0] >= selected_skip_roi[1]:
        raise ValueError("Excluded ROI minimum must be less than its maximum")

    h2uv_model = get_h2uv_model(cfg, model_key=model_key.value)
    component_cfg = get_comp_cfg(cfg, model_key=model_key.value)
    input_config = build_input_config(cfg, model_key=model_key.value)
    jags_path_value = get_jags_path(cfg)
    jags_path = Path(jags_path_value) if jags_path_value else None
    fit_root = Path(get_fit_dir(cfg))
    store_root = Path(get_artifact_store_dir(cfg))
    plots_dir = Path(get_plots_dir(cfg))
    data_name, data_blinded_name = get_data_names(cfg)
    mcfg = build_model_cfg(
        h2uv_model,
        co_cut=selected_co_cut,
        ndbd_cutoff=selected_ndbd_cutoff,
        co_shift=float(co_shift.value),
        ndbd_shift=float(ndbd_shift.value),
        skip_roi=selected_skip_roi or False,
        blinded=bool(blinded.value),
        smooth=getattr(h2uv_model, "smooth", True),
        co_cuts_dict=dict(cfg.cuts.co_cuts),
        ndbd_cuts_dict=dict(cfg.cuts.ndbd_cutoff),
    )

    input_params = InputParams(
        bin_size=float(h2uv_model.bin_size),
        esum_min=float(h2uv_model.esum_min),
        esum_max=float(h2uv_model.esum_max),
        threshold=float(h2uv_model.threshold),
        histogram_type="h2uv",
        co_cut=selected_co_cut,
        ndbd_cutoff=selected_ndbd_cutoff,
        skip_roi=selected_skip_roi,
        co_shift=float(co_shift.value),
        ndbd_shift=float(ndbd_shift.value),
        version=1,
    )
    input_id = derive_id("fit-input", input_params)
    fit_params = FitParams(
        components=tuple(c.get("name") for c in component_cfg if "name" in c),
        smooth=bool(h2uv_model.smooth),
        likelihood="poisson",
        priors={"mode": prior_width.value, "mkchain_samples": 100},
        co_cut=selected_co_cut,
        ndbd_cutoff=selected_ndbd_cutoff,
        co_shift=float(co_shift.value),
        ndbd_shift=float(ndbd_shift.value),
        skip_roi=selected_skip_roi,
        blinded=bool(blinded.value),
        version=1,
    )
    fit_id = derive_id("fit", fit_params, upstream={"input": input_id})
    store = ArtifactStore(store_root)
    input_artifact_path = store.get_path(input_id, "input.h5", artifact_type="fit-input")
    fit_artifact_dir = store.object_dir(fit_id, artifact_type="fit")

    input_command = [
        sys.executable,
        "-m",
        "qpym2.run.run_create_input",
        "--environment",
        environment.value,
        "--analysis",
        analysis.value,
        "--models",
        model_config.value,
        "--model-key",
        model_key.value,
        "--co-shifts",
        str(co_shift.value),
        "--ndbd-shift",
        str(ndbd_shift.value),
        "--legacy-path",
    ]
    if selected_co_cut is not None:
        input_command.extend(["--co-cut", selected_co_cut])
    if selected_ndbd_cutoff is not None:
        input_command.extend(["--ndbd-cutoff", selected_ndbd_cutoff])
    if selected_skip_roi is not None:
        input_command.extend(["--skip-roi", str(selected_skip_roi[0]), str(selected_skip_roi[1])])
    if blinded.value:
        input_command.append("--blinded")

    fit_command = [
        sys.executable,
        "-m",
        "qpym2.run.run_single_fit",
        "--environment",
        environment.value,
        "--analysis",
        analysis.value,
        "--models",
        model_config.value,
        "--model-key",
        model_key.value,
        "--co-shift",
        str(co_shift.value),
        "--ndbd-shift",
        str(ndbd_shift.value),
        "--prior-width",
        prior_width.value,
        "--chains",
        str(int(chains.value)),
        "--draws",
        str(int(draws.value)),
    ]
    if selected_co_cut is not None:
        fit_command.extend(["--co-cut", selected_co_cut])
    if selected_ndbd_cutoff is not None:
        fit_command.extend(["--ndbd-cutoff", selected_ndbd_cutoff])
    if selected_skip_roi is not None:
        fit_command.extend(["--skip-roi", str(selected_skip_roi[0]), str(selected_skip_roi[1])])
    if blinded.value:
        fit_command.extend(["--blinded", "true"])
    if single_core.value:
        fit_command.append("--single-core")
    if force.value:
        fit_command.append("--force")

    staged_data_name = (
        input_config.get("data_blinded_name", input_config["data_fname"])
        if blinded.value
        else input_config["data_fname"]
    )
    data_path = Path(input_config["outpath"]) / staged_data_name
    resolved_inputs = pd.DataFrame(
        [
            {"input": "staged data ROOT", "path": str(data_path), "reachable": data_path.is_file()},
            {
                "input": "JAGS background table",
                "path": str(jags_path) if jags_path else "Not configured",
                "reachable": bool(jags_path and jags_path.is_file()),
            },
            {
                "input": "fit-input artifact",
                "path": str(input_artifact_path),
                "reachable": input_artifact_path.is_file(),
            },
            {"input": "artifact store", "path": str(store_root), "reachable": store_root.is_dir()},
        ]
    )
    legacy_input_name = "input_noroi" if selected_skip_roi else "input"
    legacy_input_path = (
        fit_root / mcfg["dir_name"] / f"{legacy_input_name}{mcfg['shift_suffix']}.h5"
    )
    mo.vstack(
        [
            mo.md("## Resolved fit point and prerequisites"),
            mo.ui.table(resolved_inputs, selection=None),
            mo.md(
                f"Model bin size: **{h2uv_model.bin_size:g} keV**; ROI: "
                f"**{h2uv_model.esum_min:g}–{h2uv_model.esum_max:g} keV**.  \n"
                f"Fit input artifact: `{input_id}`  \n"
                f"Fit artifact: `{fit_id}`  \n"
                f"Legacy input path used by the fit runner: `{legacy_input_path}`  \n"
                f"Expected fit output: `{fit_artifact_dir}`"
            ),
            mo.md("Input preparation command:"),
            mo.md(f"```bash\n{shlex.join(input_command)}\n```"),
            mo.md("Fit command:"),
            mo.md(f"```bash\n{shlex.join(fit_command)}\n```"),
            mo.callout(
                mo.md(
                    "The run button first prepares or reuses this input and materializes "
                    "the legacy input file required by the single-fit runner."
                ),
                kind="info",
            ),
        ]
    )
    return (
        data_blinded_name,
        data_name,
        fit_artifact_dir,
        fit_command,
        fit_id,
        h2uv_model,
        input_artifact_path,
        input_command,
        plots_dir,
        store,
    )


app._unparsable_cell(
    r"""
    run_button = mo.ui.run_button(label="Prepare input and run this fit")
    last_run_click, set_last_run_click = mo.state(0)10.84.174.781010.84.174.78.84.174.78
    mo.vstack(
        [
            mo.md("## Run"),
            run_button,
            mo.md(
                "Fit inputs may require ROOT files and a JAGS table on the configured data paths."
            ),
        ]
    )
    """,
    name="_"
)


@app.cell
def _(
    config_dir_input,
    fit_command,
    input_command,
    last_run_click,
    mo,
    os,
    run_button,
    set_last_run_click,
    subprocess,
):
    if run_button.value and run_button.value > last_run_click():
        set_last_run_click(run_button.value)
        run_environment = os.environ.copy()
        run_environment["QPYM2_CONFIG_DIR"] = config_dir_input.value
        command_results = []
        for label, command in [
            ("Prepare fit input", input_command),
            ("Run fit", fit_command),
        ]:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                env=run_environment,
                check=False,
            )
            command_output = "\n".join(
                part for part in [result.stdout, result.stderr] if part
            ).strip()
            command_results.append((label, result.returncode, command_output))
            if result.returncode != 0:
                break

        failed = next((item for item in command_results if item[1] != 0), None)
        if failed is None:
            run_status = mo.callout(
                mo.md("Input preparation and fit runner completed successfully."),
                kind="success",
            )
        else:
            run_status = mo.callout(
                mo.md(f"{failed[0]} exited with status **{failed[1]}**."),
                kind="danger",
            )
        run_output = mo.vstack(
            [
                mo.md(f"### {label} — exit status {code}\n```text\n{output or '(no output)'}\n```")
                for label, code, output in command_results
            ]
        )
    else:
        run_status = mo.md("_No fit run requested in this session._")
        run_output = mo.md("")
    mo.vstack([run_status, run_output])
    return


@app.cell
def _(
    az,
    blinded,
    data_blinded_name,
    data_name,
    find_mode,
    fit_artifact_dir,
    fit_id,
    get_mode,
    input_artifact_path,
    mo,
    np,
    pd,
    read_fit_artifact,
    read_input,
    store,
):
    fit_trace_path = store.get_path(fit_id, "fit.nc", artifact_type="fit")
    fit_components_path = store.get_path(fit_id, "fit.h5", artifact_type="fit")
    if (
        store.exists(fit_id, artifact_type="fit")
        and fit_trace_path.is_file()
        and fit_components_path.is_file()
        and input_artifact_path.is_file()
    ):
        fit_components, trace, _ = read_fit_artifact(fullname=str(fit_artifact_dir / "fit"))
        _, signal_df, fit_data, edges = read_input(
            str(input_artifact_path),
            dataname=data_blinded_name if blinded.value else data_name,
        )
        fit_data = np.asarray(fit_data)
        variable_names = list(fit_components.index)
        posterior_modes = find_mode(trace, varnames=variable_names)
        posterior_summary = az.summary(trace, var_names=variable_names, round_to=3)
        result_table = pd.DataFrame(
            {
                "mode": posterior_modes,
                "integral": [
                    fit_components.loc[name, "integral"]
                    if "integral" in fit_components.columns
                    else None
                    for name in variable_names
                ],
                "mean": posterior_summary["mean"],
                "sd": posterior_summary["sd"],
                "r_hat": posterior_summary["r_hat"],
                "ess_bulk": posterior_summary["ess_bulk"],
            },
            index=variable_names,
        )
        divergence_count = (
            int(trace.sample_stats["diverging"].sum().values)
            if "diverging" in trace.sample_stats
            else 0
        )
        max_rhat = float(posterior_summary["r_hat"].max())
        min_ess = float(posterior_summary["ess_bulk"].min())
        if "ndbd" in trace.posterior:
            ndbd_samples = trace.posterior["ndbd"].values.flatten()
            ndbd_mode = float(get_mode(ndbd_samples))
            ndbd_p16, ndbd_p84, ndbd_p90 = np.percentile(ndbd_samples, [16, 84, 90])
            ndbd_summary = mo.md(
                f"**0νββ posterior:** mode `{ndbd_mode:.3f}` events; "
                f"central 68% interval `[{ndbd_p16:.3f}, {ndbd_p84:.3f}]`; "
                f"90th percentile `{ndbd_p90:.3f}` events."
            )
        else:
            ndbd_samples = None
            ndbd_summary = mo.md("No `ndbd` posterior variable is present in this fit.")

        convergence_kind = "success" if divergence_count == 0 and max_rhat < 1.01 else "warn"
        convergence_view = mo.callout(
            mo.md(
                f"Divergences: **{divergence_count}** · maximum r_hat: **{max_rhat:.4f}** · "
                f"minimum bulk ESS: **{min_ess:.0f}**"
            ),
            kind=convergence_kind,
        )
        results_view = mo.vstack(
            [
                mo.md(f"## Results for `{fit_id}`"),
                mo.md(f"Loaded fit files from `{fit_artifact_dir}`."),
                ndbd_summary,
                convergence_view,
                mo.ui.table(result_table, selection=None),
                mo.md(
                    f"Loaded **{len(fit_data):,} data events**, **{len(fit_components)} fitted components**, "
                    f"and **{len(edges[0]) - 1} × {len(edges[1]) - 1} histogram bins**."
                ),
            ]
        )
        posterior_norms = np.asarray(posterior_modes, dtype=float)
        results_available = True
    else:
        fit_components = None
        trace = None
        fit_data = None
        edges = None
        variable_names = []
        posterior_modes = None
        posterior_summary = None
        result_table = None
        divergence_count = None
        max_rhat = None
        min_ess = None
        ndbd_samples = None
        posterior_norms = None
        results_available = False
        results_view = mo.callout(
            mo.md(
                "The fit artifact is not available yet. Run the fit above, or check that "
                "the selected input and output paths are reachable."
            ),
            kind="info",
        )
    results_view
    return (
        fit_components,
        fit_data,
        posterior_norms,
        results_available,
        trace,
        variable_names,
    )


@app.cell
def _(
    az,
    fit_id,
    mo,
    np,
    plots_dir,
    results_available,
    save_plots,
    trace,
    variable_names,
):
    if results_available:
        trace_axes = az.plot_trace(
            trace,
            var_names=variable_names,
            compact=True,
            figsize=(14, max(3, 2.2 * len(variable_names))),
        )
        trace_fig = np.asarray(trace_axes, dtype=object).flat[0].figure
        trace_fig.suptitle("Posterior traces by component", y=1.01)
        trace_fig.tight_layout()
        if save_plots.value:
            plots_dir.mkdir(parents=True, exist_ok=True)
            trace_fig.savefig(plots_dir / f"{fit_id}_trace.png", dpi=160, bbox_inches="tight")
        trace_plot_view = trace_fig
    else:
        trace_plot_view = mo.md("_")
    trace_plot_view
    return


@app.cell
def _(
    az,
    fit_id,
    mo,
    np,
    plots_dir,
    results_available,
    save_plots,
    trace,
    variable_names,
):
    if results_available:
        pair_size = max(8, 2.2 * len(variable_names))
        pair_axes = az.plot_pair(
            trace,
            var_names=variable_names,
            kind="kde",
            marginals=True,
            figsize=(pair_size, pair_size),
        )
        pair_fig = np.asarray(pair_axes, dtype=object).flat[0].figure
        pair_fig.suptitle("Joint posterior distributions", y=1.01)
        pair_fig.tight_layout()
        if save_plots.value:
            plots_dir.mkdir(parents=True, exist_ok=True)
            pair_fig.savefig(
                plots_dir / f"{fit_id}_posterior_pair.png", dpi=160, bbox_inches="tight"
            )
        pair_plot_view = pair_fig
    else:
        pair_plot_view = mo.md("_")
    pair_plot_view
    return


@app.cell
def _(
    fit_components,
    fit_data,
    get_hist_settings,
    get_sum_hist,
    h2uv_model,
    mo,
    np,
    posterior_norms,
    results_available,
):
    if results_available:
        nbins, hist_range = get_hist_settings(h2uv_model)
        _data_hist, _x_edges, _y_edges = np.histogram2d(
            fit_data[:, 0],
            fit_data[:, 1],
            bins=nbins,
            range=hist_range,
        )
        _fit_hist = get_sum_hist(fit_components, posterior_norms)
        projection_data = (_data_hist, _fit_hist, _x_edges, _y_edges)
    else:
        projection_data = None
    projection_status = mo.md("_") if projection_data is None else mo.md("")
    projection_status
    return (projection_data,)


@app.cell
def _(
    fit_id,
    h2uv_model,
    mo,
    plot_proj_h2,
    plots_dir,
    projection_data,
    results_available,
    save_plots,
):
    if results_available:
        _data_hist, _fit_hist, _x_edges, _y_edges = projection_data
        fig_esum, _, _ = plot_proj_h2(
            _data_hist,
            _fit_hist,
            _x_edges,
            _y_edges,
            h2uv_model,
            axis="esum",
            title="E-sum projection with pulls",
            log_scale=True,
        )
        if save_plots.value:
            plots_dir.mkdir(parents=True, exist_ok=True)
            fig_esum.savefig(
                plots_dir / f"{fit_id}_esum_projection.png", dpi=160, bbox_inches="tight"
            )
        esum_plot_view = fig_esum
    else:
        esum_plot_view = mo.md("_")
    esum_plot_view
    return


@app.cell
def _(
    fit_id,
    h2uv_model,
    mo,
    plot_proj_h2,
    plots_dir,
    projection_data,
    results_available,
    save_plots,
):
    if results_available:
        _data_hist, _fit_hist, _x_edges, _y_edges = projection_data
        bin_size = float(getattr(h2uv_model, "bin_size", getattr(h2uv_model, "binsize", 1.0)))
        fig_ediff, _, _ = plot_proj_h2(
            _data_hist,
            _fit_hist,
            _x_edges,
            _y_edges,
            h2uv_model,
            axis="ediff",
            rebin=max(1, int(round(10.0 / bin_size))),
            title="E-difference projection with pulls",
            log_scale=False,
        )
        if save_plots.value:
            plots_dir.mkdir(parents=True, exist_ok=True)
            fig_ediff.savefig(
                plots_dir / f"{fit_id}_ediff_projection.png",
                dpi=160,
                bbox_inches="tight",
            )
        ediff_plot_view = fig_ediff
    else:
        ediff_plot_view = mo.md("_")
    ediff_plot_view
    return


@app.cell
def _():
    from pydantic_ai import Agent

    # Try switching from 'openai:gpt-6-luna' to 'openai-chat:gpt-6-luna'
    agent = Agent("openai-chat:gpt-6-luna")

    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
