import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    import os
    import sys
    from pathlib import Path

    import arviz as az
    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    import pandas as pd

    from qpym2.cfg.config_builders_new import (
        get_comp_cfg,
        get_fit_data_dir,
        get_h2uv_model,
        get_jags_path,
        get_model_storage_coords,
    )
    from qpym2.cfg.hydra_adapter_new import (
        HydraConfigManager,
        get_config,
        resolve_config_dir,
    )
    from qpym2.run.run_single_storage_fit import run_single_storage_fit
    from qpym2.storage.backend import fit_exists, read_fit, read_observed_data
    from qpym2.storage.fits import (
        build_fit_hash,
        resolve_fit_template_uuids,
        setup_storage_fit,
    )
    from qpym2.storage.inputs import resolve_storage_coords
    from qpym2.storage.paths import (
        get_output_dir,
        get_parameter_dir,
    )
    from qpym2.studies import assemble_study

    return (
        HydraConfigManager,
        Path,
        assemble_study,
        az,
        build_fit_hash,
        fit_exists,
        get_comp_cfg,
        get_config,
        get_fit_data_dir,
        get_h2uv_model,
        get_jags_path,
        get_model_storage_coords,
        get_output_dir,
        get_parameter_dir,
        json,
        mo,
        np,
        os,
        pd,
        plt,
        read_fit,
        read_observed_data,
        resolve_config_dir,
        resolve_fit_template_uuids,
        resolve_storage_coords,
        run_single_storage_fit,
        setup_storage_fit,
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
                # QPyM2 v0.6: Bayesian Fit Execution & Inspection Workbench

                Inspect and execute PyMC Bayesian fits within the modern QPyM2 storage architecture.
                Features deterministic parameter-set UUID linking, `<output_UUID>` caching,
                posterior convergence inspection, 2D/1D data vs fit comparisons, and study tracking.
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

    env_options = config_options.get("environments", ["local"])
    environment = mo.ui.dropdown(
        options=env_options,
        value="local" if "local" in env_options else env_options[0],
        label="Environment",
    )
    analysis_options = config_options.get("analysis", ["base"])
    analysis = mo.ui.dropdown(
        options=analysis_options,
        value="base" if "base" in analysis_options else analysis_options[0],
        label="Analysis profile",
    )
    model_sets = config_options.get("models", ["tech_note"])
    model_set = mo.ui.dropdown(
        options=model_sets,
        value="tech_note" if "tech_note" in model_sets else model_sets[0],
        label="Model config set",
    )

    mo.hstack([environment, analysis, model_set])
    return analysis, environment, model_set


@app.cell
def _(analysis, config_dir_input, environment, get_config, mo, model_set):
    cfg = get_config(
        environment=environment.value,
        analysis=analysis.value,
        models=model_set.value,
        config_dir=config_dir_input.value,
    )
    model_keys = list(cfg.models.keys())
    model_key = mo.ui.dropdown(
        options=model_keys,
        value="tech_note_200" if "tech_note_200" in model_keys else model_keys[0],
        label="Model entry",
    )
    mo.hstack([mo.md("**Model Entry:**"), model_key])
    return cfg, model_key


@app.cell
def _(cfg, mo):
    co_cuts = ["None", *list(getattr(cfg.cuts, "co_cuts", {}).keys())]
    ndbd_cuts = ["None", *list(getattr(cfg.cuts, "ndbd_cutoff", {}).keys())]

    co_cut = mo.ui.dropdown(
        options=co_cuts,
        value="gamma10_v150" if "gamma10_v150" in co_cuts else co_cuts[0],
        label="Co-60 cut",
    )
    ndbd_cutoff = mo.ui.dropdown(
        options=ndbd_cuts,
        value="u1740" if "u1740" in ndbd_cuts else ndbd_cuts[0],
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
    blinded = mo.ui.checkbox(value=False, label="Blinded data")
    prior_width = mo.ui.dropdown(
        options=["5sig", "10sig", "syst", "2syst"],
        value="5sig",
        label="Prior width",
    )
    chains = mo.ui.number(value=8, start=1, step=1, label="Chains")
    draws = mo.ui.number(value=20000, start=100, step=500, label="Draws per chain")
    single_core = mo.ui.checkbox(value=False, label="Single core")
    save_trace = mo.ui.checkbox(value=True, label="Save .nc trace")
    force = mo.ui.checkbox(value=False, label="Force re-run (ignore cache)")

    mo.vstack(
        [
            mo.md("### Parameter Point & MCMC Options"),
            mo.hstack([co_cut, ndbd_cutoff, co_shift, ndbd_shift]),
            mo.hstack([blinded, prior_width, chains, draws]),
            mo.hstack([single_core, save_trace, force]),
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
        save_trace,
        single_core,
    )


@app.cell
def _(
    blinded,
    cfg,
    chains,
    co_cut,
    co_shift,
    draws,
    mo,
    model_key,
    ndbd_cutoff,
    ndbd_shift,
    pd,
    prior_width,
    setup_storage_fit,
):
    # Preview parameter hashes and target paths without running sampler
    sel_co_cut = None if co_cut.value == "None" else co_cut.value
    sel_ndbd_cut = None if ndbd_cutoff.value == "None" else ndbd_cutoff.value

    fit_setup = setup_storage_fit(
        cfg=cfg,
        model_key=model_key.value,
        co_cut=sel_co_cut,
        ndbd_cutoff=sel_ndbd_cut,
        co_shift=float(co_shift.value),
        ndbd_shift=float(ndbd_shift.value),
        blinded=blinded.value,
        prior_width=prior_width.value,
        chains=int(chains.value),
        draws=int(draws.value),
        dry_run=True,
    )

    output_uuid = fit_setup["output_uuid"]
    fit_prefix = fit_setup["fit_prefix"]
    h5_exists = fit_setup["is_cached"]
    nc_exists = fit_setup["nc_path"].exists()
    param_dir = fit_setup["param_dir"]
    uuid_base = fit_setup["uuid_base"]
    uuid_co = fit_setup["uuid_co"]
    uuid_ndbd = fit_setup["uuid_ndbd"]
    fit_params = fit_setup["fit_params"]

    run_btn = mo.ui.run_button(label="Execute / Load Fit")

    info_df = pd.DataFrame(
        [
            {"Property": "Parameter Directory", "Value": str(param_dir)},
            {"Property": "Base UUID (invariants)", "Value": uuid_base},
            {"Property": "Co-60 UUID", "Value": uuid_co},
            {"Property": "0νββ UUID", "Value": uuid_ndbd},
            {"Property": "Output UUID", "Value": output_uuid},
            {
                "Property": "Fit HDF5 Target",
                "Value": f"{fit_prefix.name}.h5 ({'EXISTS' if h5_exists else 'MISSING'})",
            },
            {
                "Property": "Trace NetCDF Target",
                "Value": f"{fit_prefix.name}.nc ({'EXISTS' if nc_exists else 'MISSING'})",
            },
        ]
    )

    mo.vstack(
        [
            mo.md("### Output Hash Preview & Cache Status"),
            mo.ui.table(info_df),
            run_btn,
        ]
    )
    return (
        fit_params,
        fit_prefix,
        h5_exists,
        nc_exists,
        output_uuid,
        param_dir,
        run_btn,
        sel_co_cut,
        sel_ndbd_cut,
        uuid_base,
        uuid_co,
        uuid_ndbd,
    )


@app.cell
def _(
    blinded,
    cfg,
    chains,
    co_shift,
    draws,
    fit_prefix,
    force,
    h5_exists,
    mo,
    model_key,
    ndbd_shift,
    prior_width,
    read_fit,
    run_btn,
    run_single_storage_fit,
    save_trace,
    sel_co_cut,
    sel_ndbd_cut,
    single_core,
):
    if not run_btn.value and not h5_exists:
        mo.output.replace(mo.md("*Click 'Execute / Load Fit' to run or inspect fit results.*"))
        table, trace, data, metadata = None, None, None, None
    else:
        if h5_exists and not force.value:
            mo.output.replace(mo.md(f"**Loaded cached fit from:** `{fit_prefix}.h5`"))
            table, trace, data, metadata = read_fit(fit_prefix)
        else:
            mo.output.replace(mo.md(f"**Running Bayesian fit for {fit_prefix.name}...**"))
            res = run_single_storage_fit(
                cfg=cfg,
                model_key=model_key.value,
                co_cut=sel_co_cut,
                ndbd_cutoff=sel_ndbd_cut,
                co_shift=float(co_shift.value),
                ndbd_shift=float(ndbd_shift.value),
                blinded=blinded.value,
                prior_width=prior_width.value,
                chains=int(chains.value),
                draws=int(draws.value),
                single_core=single_core.value,
                save_trace=save_trace.value,
                force=force.value,
                verbose=False,
            )
            table, trace, data, metadata = (
                res.get("table"),
                res.get("trace"),
                res.get("data"),
                res.get("metadata"),
            )

    return data, metadata, table, trace


@app.cell
def _(mo, np, pd, table, trace):
    if table is None:
        mo.output.replace(mo.md(""))
        res_df = None
    else:
        _records = []
        for _var in table.index:
            _rec = {"Component": _var}
            if "integral" in table.columns:
                _rec["Prior Integral"] = f"{float(table.loc[_var, 'integral']):.1f}"
            if trace is not None and hasattr(trace, "posterior") and _var in trace.posterior:
                _vals = trace.posterior[_var].values
                _rec["Posterior Mean"] = f"{float(_vals.mean()):.2f}"
                _rec["Posterior Std"] = f"{float(_vals.std()):.2f}"
                _rec["90% Interval"] = (
                    f"[{float(np.percentile(_vals, 5)):.1f}, {float(np.percentile(_vals, 95)):.1f}]"
                )
            _records.append(_rec)

        res_df = pd.DataFrame(_records)
        mo.output.replace(
            mo.vstack(
                [
                    mo.md("### Fitted Component Summary"),
                    mo.ui.table(res_df),
                ]
            )
        )
    return (res_df,)


@app.cell
def _(mo, table, trace):
    if trace is None or not hasattr(trace, "posterior"):
        mo.output.replace(mo.md("*Trace (.nc) not available for detailed posterior plots.*"))
        plot_var = None
    else:
        _var_options = [v for v in table.index if v in trace.posterior]
        plot_var = mo.ui.dropdown(
            options=_var_options,
            value=_var_options[0] if _var_options else None,
            label="Plot variable",
        )
        mo.output.replace(mo.hstack([mo.md("**Posterior Diagnostics:**"), plot_var]))
    return (plot_var,)


@app.cell
def _(np, plot_var, plt, trace):
    if plot_var is None or plot_var.value is None or trace is None:
        diag_fig = None
    else:
        _vname = plot_var.value
        _vals = trace.posterior[_vname].values.flatten()
        diag_fig, _axes = plt.subplots(1, 2, figsize=(11, 4))

        # Posterior distribution
        _axes[0].hist(_vals, bins=40, color="steelblue", edgecolor="black", alpha=0.7, density=True)
        _q05, _q50, _q95 = np.percentile(_vals, [5, 50, 95])
        _axes[0].axvline(_q50, color="darkred", linestyle="--", label=f"Median: {_q50:.1f}")
        _axes[0].axvline(
            _q05, color="gray", linestyle=":", label=f"90% CI: [{_q05:.1f}, {_q95:.1f}]"
        )
        _axes[0].axvline(_q95, color="gray", linestyle=":")
        _axes[0].set_xlabel(f"{_vname} events")
        _axes[0].set_ylabel("Posterior Density")
        _axes[0].set_title(f"Posterior Distribution: {_vname}")
        _axes[0].legend(frameon=True)

        # MCMC chain trace
        for _c in range(trace.posterior.sizes.get("chain", 1)):
            _chain_vals = trace.posterior[_vname].sel(chain=_c).values
            _axes[1].plot(_chain_vals, alpha=0.6, label=f"Chain {_c}")
        _axes[1].set_xlabel("Draw")
        _axes[1].set_ylabel("Sample Value")
        _axes[1].set_title(f"Trace Plot: {_vname}")

        plt.tight_layout()
    diag_fig
    return (diag_fig,)


@app.cell
def _(assemble_study, mo, pd, plt):
    study_input = mo.ui.text(
        value="studies/example_co_shift", label="Study directory", full_width=True
    )
    check_study_btn = mo.ui.run_button(label="Inspect Study")

    mo.vstack(
        [
            mo.md("### Study Inspection & Result Aggregation"),
            study_input,
            check_study_btn,
        ]
    )
    return (check_study_btn, study_input)


@app.cell
def _(Path, assemble_study, check_study_btn, mo, pd, plt, study_input):
    if not check_study_btn.value:
        mo.output.replace(mo.md("*Click 'Inspect Study' to assemble and view study points.*"))
    else:
        _sdir = Path(study_input.value)
        if not _sdir.exists():
            mo.output.replace(mo.md(f"⚠️ Study directory `{_sdir}` does not exist."))
        else:
            study = assemble_study(_sdir)
            fits_df = study.check_fits()

            _res_csv = _sdir / "results.csv"
            res_plot = None
            if _res_csv.exists():
                _rdf = pd.read_csv(_res_csv)
                if "co_shift" in _rdf.columns and "ndbd_mean" in _rdf.columns:
                    study_fig, _ax = plt.subplots(figsize=(6, 3.5))
                    _ax.errorbar(
                        _rdf["co_shift"],
                        _rdf["ndbd_mean"],
                        yerr=_rdf["ndbd_std"] if "ndbd_std" in _rdf.columns else None,
                        fmt="o-",
                        capsize=3,
                        color="navy",
                    )
                    _ax.set_xlabel("Co-60 Shift [keV]")
                    _ax.set_ylabel("Fitted 0νββ Events")
                    _ax.set_title(f"{study.name}: 0νββ Yield vs Co-60 Shift")
                    plt.tight_layout()
                    res_plot = study_fig

            mo.output.replace(
                mo.vstack(
                    [
                        mo.md(
                            f"**Assembled Study:** `{study.name}` ({len(study)} evaluation points)"
                        ),
                        mo.ui.table(
                            fits_df[
                                [
                                    "task_id",
                                    "point_id",
                                    "co_shift",
                                    "output_uuid",
                                    "is_fitted",
                                    "status",
                                ]
                            ]
                        ),
                        res_plot
                        if res_plot is not None
                        else mo.md(
                            "*No results.csv found. Run study to generate aggregated plot.*"
                        ),
                    ]
                )
            )
    return


if __name__ == "__main__":
    app.run()
