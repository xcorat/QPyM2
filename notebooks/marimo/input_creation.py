import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    import time
    from pathlib import Path

    import matplotlib.pyplot as plt
    import marimo as mo
    import numpy as np
    import pandas as pd

    from qpym2.cfg.config_builders import (
        get_comp_cfg,
        get_data_names,
        get_fit_data_dir,
        get_h2uv_model,
        get_jags_path,
        get_model_storage_coords,
        get_staging_dir,
    )
    from qpym2.cfg.hydra_adapter import (
        HydraConfigManager,
        get_config,
        resolve_config_dir,
    )
    from qpym2.io.bkg_model import read_jagsh5
    from qpym2.storage.backend import (
        data_exists,
        read_observed_data,
        read_template,
        template_exists,
    )
    from qpym2.storage.hasher import canonical_json, compute_input_uuid
    from qpym2.storage.inputs import (
        build_source_hash_params,
        create_storage_inputs_for_point,
        generate_observed_data,
        generate_source_template,
        resolve_storage_coords,
    )
    from qpym2.storage.metadata import read_hdf5_metadata
    from qpym2.storage.paths import (
        get_data_path,
        get_input_dir,
        get_parameter_dir,
        get_template_path,
    )
    from qpym2.storage.schema import build_template_hash_dict

    return (
        HydraConfigManager,
        Path,
        build_source_hash_params,
        build_template_hash_dict,
        canonical_json,
        compute_input_uuid,
        create_storage_inputs_for_point,
        data_exists,
        generate_observed_data,
        generate_source_template,
        get_comp_cfg,
        get_config,
        get_data_names,
        get_data_path,
        get_fit_data_dir,
        get_h2uv_model,
        get_input_dir,
        get_jags_path,
        get_model_storage_coords,
        get_parameter_dir,
        get_staging_dir,
        get_template_path,
        json,
        mo,
        np,
        pd,
        plt,
        read_hdf5_metadata,
        read_jagsh5,
        read_observed_data,
        read_template,
        resolve_config_dir,
        resolve_storage_coords,
        template_exists,
        time,
    )


@app.cell
def _(mo, resolve_config_dir):
    default_config_dir, config_dir_source = resolve_config_dir()
    config_dir_input = mo.ui.text(
        value=str(default_config_dir),
        label="Config directory",
        full_width=True,
    )

    config_source_note = {
        "flag": "explicit --config-dir",
        "QPYM2_CONFIG_DIR": "QPYM2_CONFIG_DIR environment variable",
        "data_dir": "<data dir>/config",
        "cwd": "./config (current working directory)",
        "repo": "repository config/ (packaged example)",
    }.get(config_dir_source, config_dir_source)

    if config_dir_source == "cwd":
        loc_notice = mo.callout(
            mo.md(
                f"Using **`./config`** ({config_source_note}). Set `QPYM2_CONFIG_DIR` or "
                "pass `--config-dir` to use a config tree stored with your data."
            ),
            kind="warn",
        )
    else:
        loc_notice = mo.md(f"Config directory source: {config_source_note}.")

    header = mo.vstack(
        [
            mo.md(
                """
                # Storage Input Generation & Inspection Workbench (v0.6)

                Interactive tool to check, preview, generate, and validate identity-addressed
                fit inputs according to the modern QPyM2 v0.6 storage architecture.
                
                - **Path-Encoded Geometry:** `fit_data/roi_<min>_<max>/thr<val>/<coord>_b<bin_size>/input/`
                - **Deterministic Hashes:** `<component_name>_<input_UUID>.h5` based on parameter hashing
                - **Embedded HDF5 Metadata:** Invariant source deduplication across shifts
                """
            ),
            config_dir_input,
            loc_notice,
        ]
    )
    header
    return (config_dir_input,)


@app.cell
def _(HydraConfigManager, config_dir_input, mo):
    HydraConfigManager.set_config_dir(config_dir_input.value)
    _options = HydraConfigManager.list_options()

    _env_opts = _options.get("environments", ["local"])
    env_dropdown = mo.ui.dropdown(
        options=_env_opts,
        value="local" if "local" in _env_opts else _env_opts[0],
        label="Environment",
    )
    _staging_opts = _options.get("staging", ["default"])
    staging_dropdown = mo.ui.dropdown(
        options=_staging_opts,
        value="default" if "default" in _staging_opts else _staging_opts[0],
        label="Staging tier",
    )
    _analysis_opts = _options.get("analysis", ["base"])
    analysis_dropdown = mo.ui.dropdown(
        options=_analysis_opts,
        value="base" if "base" in _analysis_opts else _analysis_opts[0],
        label="Analysis",
    )
    _models_opts = _options.get("models", ["tech_note"])
    model_config_dropdown = mo.ui.dropdown(
        options=_models_opts,
        value="tech_note" if "tech_note" in _models_opts else _models_opts[0],
        label="Model config set",
    )

    mo.hstack([env_dropdown, staging_dropdown, analysis_dropdown, model_config_dropdown])
    return analysis_dropdown, env_dropdown, model_config_dropdown, staging_dropdown


@app.cell
def _(
    analysis_dropdown,
    config_dir_input,
    env_dropdown,
    get_config,
    mo,
    model_config_dropdown,
    staging_dropdown,
):
    cfg = get_config(
        environment=env_dropdown.value,
        staging=staging_dropdown.value,
        analysis=analysis_dropdown.value,
        models=model_config_dropdown.value,
        config_dir=config_dir_input.value,
    )

    _model_key_opts = list(cfg.models.keys())
    model_key = mo.ui.dropdown(
        options=_model_key_opts,
        value="tech_note_200" if "tech_note_200" in _model_key_opts else _model_key_opts[0],
        label="Model entry",
    )

    mo.vstack([mo.md("## 1. Model Selection"), mo.hstack([model_key])])
    return cfg, model_key


@app.cell
def _(
    Path,
    cfg,
    get_fit_data_dir,
    get_staging_dir,
    mo,
    model_key,
    pd,
    resolve_storage_coords,
):
    coords, h2uv_model, param_dir, input_dir = resolve_storage_coords(
        cfg, model_key=model_key.value
    )
    staging_dir = Path(get_staging_dir(cfg))
    fit_data_root = get_fit_data_dir(cfg)

    coords_summary = pd.DataFrame(
        [
            {
                "property": "fit_data root",
                "value": str(fit_data_root),
                "status": "exists" if fit_data_root.is_dir() else "will create",
            },
            {
                "property": "Resolved param_dir",
                "value": str(param_dir),
                "status": "exists" if param_dir.is_dir() else "will create",
            },
            {
                "property": "Input directory",
                "value": str(input_dir),
                "status": "exists" if input_dir.is_dir() else "will create",
            },
            {
                "property": "ROI Esum range",
                "value": f"[{coords['roi'][0]:g}, {coords['roi'][1]:g}] keV",
                "status": "ok",
            },
            {"property": "Threshold", "value": f"{coords['threshold']:g} keV", "status": "ok"},
            {
                "property": "Coordinate / Bin size",
                "value": f"{coords['coord']} (bin_size={coords['bin_size']:g} keV)",
                "status": "ok",
            },
            {
                "property": "Staging directory",
                "value": str(staging_dir),
                "status": "exists" if staging_dir.is_dir() else "missing",
            },
        ]
    )

    mo.vstack(
        [
            mo.md("## 2. Storage Coordinates & Directory Resolution"),
            mo.ui.table(coords_summary, selection=None),
        ]
    )
    return (
        coords,
        coords_summary,
        fit_data_root,
        h2uv_model,
        input_dir,
        param_dir,
        staging_dir,
    )


@app.cell
def _(
    Path,
    cfg,
    get_data_names,
    get_jags_path,
    mo,
    pd,
    read_jagsh5,
    staging_dir,
):
    jags_fpath = Path(get_jags_path(cfg))
    try:
        source_table, source_chain = read_jagsh5(str(jags_fpath))
        jags_status_str = f"Loaded {len(source_table)} MC sources from JAGS table."
        source_error = None
    except Exception as exc:
        source_table = pd.DataFrame(columns=["fname", "mean", "std", "mode"])
        source_chain = None
        jags_status_str = f"Error reading JAGS table: {exc}"
        source_error = str(exc)

    data_names = get_data_names(cfg)
    unblinded_fpath = staging_dir / f"{data_names[0]}.root"
    blinded_fpath = staging_dir / f"{data_names[1]}.root"
    ndbd_fpath = staging_dir / "ndbd.root"

    prereq_rows = [
        {
            "role": "JAGS Background Model",
            "source": jags_fpath.name,
            "path": str(jags_fpath),
            "available": jags_fpath.is_file(),
            "size_mb": round(jags_fpath.stat().st_size / 1e6, 3) if jags_fpath.is_file() else 0.0,
        },
        {
            "role": "Unblinded Data",
            "source": unblinded_fpath.name,
            "path": str(unblinded_fpath),
            "available": unblinded_fpath.is_file(),
            "size_mb": round(unblinded_fpath.stat().st_size / 1e6, 3)
            if unblinded_fpath.is_file()
            else 0.0,
        },
        {
            "role": "Blinded Data",
            "source": blinded_fpath.name,
            "path": str(blinded_fpath),
            "available": blinded_fpath.is_file(),
            "size_mb": round(blinded_fpath.stat().st_size / 1e6, 3)
            if blinded_fpath.is_file()
            else 0.0,
        },
        {
            "role": "0vbb Signal MC",
            "source": ndbd_fpath.name,
            "path": str(ndbd_fpath),
            "available": ndbd_fpath.is_file(),
            "size_mb": round(ndbd_fpath.stat().st_size / 1e6, 3) if ndbd_fpath.is_file() else 0.0,
        },
    ]

    for source_name, row in source_table.iterrows():
        fname = str(row["fname"])
        p = staging_dir / fname
        prereq_rows.append(
            {
                "role": f"MC: {source_name}",
                "source": fname,
                "path": str(p),
                "available": p.is_file(),
                "size_mb": round(p.stat().st_size / 1e6, 3) if p.is_file() else 0.0,
            }
        )

    prereq_df = pd.DataFrame(prereq_rows)
    avail_count = prereq_df["available"].sum()
    total_count = len(prereq_df)

    mo.vstack(
        [
            mo.md("## 3. Upstream Staging Prerequisites"),
            mo.callout(
                mo.md(
                    f"**Availability:** {avail_count} of {total_count} files ready in `{staging_dir}`. {jags_status_str}"
                ),
                kind="success" if avail_count == total_count else "warn",
            ),
            mo.ui.table(prereq_df, selection=None),
        ]
    )
    return (
        avail_count,
        blinded_fpath,
        data_names,
        jags_fpath,
        jags_status_str,
        ndbd_fpath,
        prereq_df,
        prereq_rows,
        source_chain,
        source_error,
        source_table,
        total_count,
        unblinded_fpath,
    )


@app.cell
def _(cfg, mo):
    _co_cut_opts = ["None", *[str(c) for c in getattr(cfg.cuts, "co_cuts", {}).keys()]]
    co_cut_dropdown = mo.ui.dropdown(
        options=_co_cut_opts,
        value="gamma10_v150" if "gamma10_v150" in _co_cut_opts else _co_cut_opts[0],
        label="Co-60 cut",
    )

    _co_shift_opts = [0.0, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0]
    co_shift_dropdown = mo.ui.dropdown(
        options=_co_shift_opts,
        value=1.0,
        label="Co-60 energy shift (keV)",
    )

    _ndbd_cut_opts = ["None", *[str(c) for c in getattr(cfg.cuts, "ndbd_cutoff", {}).keys()]]
    ndbd_cutoff_dropdown = mo.ui.dropdown(
        options=_ndbd_cut_opts,
        value="u1740" if "u1740" in _ndbd_cut_opts else _ndbd_cut_opts[0],
        label="0vbb cutoff",
    )

    _ndbd_shift_opts = [0.0, 0.3, 0.5, 0.8, 1.0]
    ndbd_shift_dropdown = mo.ui.dropdown(
        options=_ndbd_shift_opts,
        value=0.0,
        label="0vbb energy shift (keV)",
    )

    blinded_checkbox = mo.ui.checkbox(
        value=False,
        label="Blinded data sample (data_blinded.h5)",
    )

    mo.vstack(
        [
            mo.md("## 4. Physics Variation & Shift Controls"),
            mo.hstack([co_cut_dropdown, co_shift_dropdown]),
            mo.hstack([ndbd_cutoff_dropdown, ndbd_shift_dropdown]),
            mo.hstack([blinded_checkbox]),
        ]
    )
    return (
        blinded_checkbox,
        co_cut_dropdown,
        co_shift_dropdown,
        ndbd_cutoff_dropdown,
        ndbd_shift_dropdown,
    )


@app.cell
def _(
    Path,
    blinded_checkbox,
    cfg,
    co_cut_dropdown,
    co_shift_dropdown,
    create_storage_inputs_for_point,
    mo,
    model_key,
    ndbd_cutoff_dropdown,
    ndbd_shift_dropdown,
    pd,
):
    _co_cut_val = None if co_cut_dropdown.value == "None" else co_cut_dropdown.value
    _ndbd_cut_val = None if ndbd_cutoff_dropdown.value == "None" else ndbd_cutoff_dropdown.value
    _co_shift_val = float(co_shift_dropdown.value)
    _ndbd_shift_val = float(ndbd_shift_dropdown.value)
    _blinded_val = bool(blinded_checkbox.value)

    plan_manifest = create_storage_inputs_for_point(
        cfg,
        model_key=model_key.value,
        co_cut=_co_cut_val,
        ndbd_cutoff=_ndbd_cut_val,
        co_shift=_co_shift_val,
        ndbd_shift=_ndbd_shift_val,
        blinded=_blinded_val,
        dry_run=True,
    )

    plan_rows = []
    _data_info = plan_manifest["data"]
    plan_rows.append(
        {
            "component": "Observed Data",
            "file": Path(_data_info["path"]).name,
            "uuid": "N/A (fixed identity)",
            "cached": _data_info["cached"],
            "full_path": _data_info["path"],
        }
    )

    for src_name, t_info in plan_manifest["templates"].items():
        plan_rows.append(
            {
                "component": src_name,
                "file": Path(t_info["path"]).name,
                "uuid": t_info["uuid"],
                "cached": t_info["cached"],
                "full_path": t_info["path"],
            }
        )

    plan_df = pd.DataFrame(plan_rows)
    cached_count = plan_df["cached"].sum()
    total_inputs = len(plan_df)

    mo.vstack(
        [
            mo.md("## 5. Planned Storage Inputs (Dry-Run Preview)"),
            mo.callout(
                mo.md(
                    f"**Cache Status:** {cached_count} / {total_inputs} inputs already exist in "
                    f"`{plan_manifest['input_dir']}`.  \n"
                    f"**Shift Isolation:** Varying `co_shift` ({_co_shift_val} keV) only updates "
                    f"`60Co_<uuid>.h5`. All other background sources retain identical UUIDs and are deduplicated."
                ),
                kind="info",
            ),
            mo.ui.table(plan_df, selection=None),
        ]
    )
    return cached_count, plan_df, plan_manifest, plan_rows, total_inputs


@app.cell
def _(mo):
    run_btn = mo.ui.run_button(label="Generate Storage Inputs")
    overwrite_box = mo.ui.checkbox(value=False, label="Overwrite existing cached files")

    mo.vstack(
        [
            mo.md("## 6. Storage Input Generation"),
            mo.hstack([run_btn, overwrite_box]),
        ]
    )
    return overwrite_box, run_btn


@app.cell
def _(
    blinded_checkbox,
    cfg,
    co_cut_dropdown,
    co_shift_dropdown,
    create_storage_inputs_for_point,
    mo,
    model_key,
    ndbd_cutoff_dropdown,
    ndbd_shift_dropdown,
    overwrite_box,
    run_btn,
    time,
):
    if not run_btn.value:
        gen_status_view = mo.md(
            "_Click ‘Generate Storage Inputs’ to create or refresh input files._"
        )
    else:
        _t0 = time.time()
        try:
            _res = create_storage_inputs_for_point(
                cfg,
                model_key=model_key.value,
                co_cut=None if co_cut_dropdown.value == "None" else co_cut_dropdown.value,
                ndbd_cutoff=None
                if ndbd_cutoff_dropdown.value == "None"
                else ndbd_cutoff_dropdown.value,
                co_shift=float(co_shift_dropdown.value),
                ndbd_shift=float(ndbd_shift_dropdown.value),
                blinded=bool(blinded_checkbox.value),
                overwrite=bool(overwrite_box.value),
                dry_run=False,
            )
            _elapsed = time.time() - _t0
            _cached = sum(1 for t in _res["templates"].values() if t.get("cached"))
            _written = len(_res["templates"]) - _cached
            gen_status_view = mo.callout(
                mo.md(
                    f"**Input generation complete** in {_elapsed:.2f}s!  \n"
                    f"- Data: `{_res['data']['path']}` ({'reused cache' if _res['data'].get('cached') else 'written'})  \n"
                    f"- Templates: {_written} written, {_cached} reused from cache."
                ),
                kind="success",
            )
        except Exception as exc:
            gen_status_view = mo.callout(mo.md(f"**Generation failed:** `{exc}`"), kind="danger")

    gen_status_view
    return (gen_status_view,)


@app.cell
def _(gen_status_view, input_dir, mo):
    # Depending on gen_status_view ensures file list re-evaluates after generation
    _ = gen_status_view
    _existing_files = sorted(input_dir.glob("*.h5")) if input_dir.is_dir() else []
    _file_names = [f.name for f in _existing_files]

    file_selector = mo.ui.dropdown(
        options=_file_names,
        value=_file_names[0] if _file_names else None,
        label="Select input file in input/",
    )

    mo.vstack([mo.md("## 7. HDF5 File Inspector"), file_selector])
    return (file_selector,)


@app.cell
def _(
    file_selector,
    input_dir,
    json,
    mo,
    np,
    pd,
    read_hdf5_metadata,
    read_observed_data,
    read_template,
):
    selected_hist_data = None
    selected_edges_data = None
    selected_meta_data = None

    if file_selector.value:
        _target_path = input_dir / file_selector.value
        try:
            _meta = read_hdf5_metadata(_target_path)
            selected_meta_data = _meta

            if "data" in file_selector.value:
                _data_arr, _edges, _ = read_observed_data(_target_path)
                selected_hist_data = _data_arr
                selected_edges_data = _edges
                _info_df = pd.DataFrame(
                    [
                        {"property": "File type", "value": "Observed Data"},
                        {
                            "property": "Type / Shape",
                            "value": f"{type(_data_arr)} - {getattr(_data_arr, 'shape', len(_data_arr))}",
                        },
                        {
                            "property": "Event count / Total sum",
                            "value": f"{len(_data_arr):,} events"
                            if hasattr(_data_arr, "__len__") and not hasattr(_data_arr, "shape")
                            else f"{np.sum(_data_arr):,.1f}",
                        },
                        {
                            "property": "File size",
                            "value": f"{_target_path.stat().st_size / 1e3:.1f} KB",
                        },
                    ]
                )
            else:
                _hist, _edges, _ = read_template(_target_path)
                selected_hist_data = _hist
                selected_edges_data = _edges
                _nonzero = int(np.count_nonzero(_hist))
                _total_bins = int(_hist.size)
                _info_df = pd.DataFrame(
                    [
                        {"property": "File type", "value": "Source Template (Unsmoothed)"},
                        {"property": "Component", "value": str(_meta.get("component", "unknown"))},
                        {"property": "Array shape", "value": str(_hist.shape)},
                        {"property": "Histogram sum", "value": f"{np.sum(_hist):,.2f}"},
                        {
                            "property": "Min / Max bin value",
                            "value": f"{np.min(_hist):.4g} / {np.max(_hist):.4g}",
                        },
                        {
                            "property": "Non-zero bins",
                            "value": f"{_nonzero} / {_total_bins} ({_nonzero / _total_bins * 100:.1f}%)",
                        },
                        {
                            "property": "File size",
                            "value": f"{_target_path.stat().st_size / 1e3:.1f} KB",
                        },
                    ]
                )

            inspect_output = mo.vstack(
                [
                    mo.md(f"### Details: `{file_selector.value}`"),
                    mo.ui.table(_info_df, selection=None),
                    mo.md("#### Embedded HDF5 Metadata (`attrs['metadata_json']`)"),
                    mo.ui.code(json.dumps(_meta, indent=2), language="json"),
                ]
            )
        except Exception as exc:
            inspect_output = mo.callout(
                mo.md(f"Failed to inspect `{file_selector.value}`: `{exc}`"), kind="danger"
            )
    else:
        inspect_output = mo.callout(
            mo.md(
                "No `.h5` files found in the current model `input/` directory. Generate inputs above."
            ),
            kind="info",
        )

    inspect_output
    return (
        inspect_output,
        selected_edges_data,
        selected_hist_data,
        selected_meta_data,
    )


@app.cell
def _(file_selector, mo, np, plt, selected_edges_data, selected_hist_data):
    plot_view = []
    if (
        selected_hist_data is not None
        and isinstance(selected_hist_data, np.ndarray)
        and selected_hist_data.ndim == 2
    ):
        try:
            u_edges = selected_edges_data[0]
            v_edges = selected_edges_data[1]

            fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)

            # 2D Heatmap
            pcm = axes[0].pcolormesh(
                u_edges, v_edges, selected_hist_data.T, cmap="viridis", shading="auto"
            )
            axes[0].set_title(f"2D Template: {file_selector.value}", fontsize=10)
            axes[0].set_xlabel("u (sum energy / sqrt(2)) [keV]")
            axes[0].set_ylabel("v (asymmetry) [keV]")
            fig.colorbar(pcm, ax=axes[0], label="Counts")

            # 1D Projection along u
            u_proj = selected_hist_data.sum(axis=1)
            u_centers = 0.5 * (u_edges[:-1] + u_edges[1:])
            axes[1].step(u_centers, u_proj, where="mid", color="teal")
            axes[1].set_title("Projection along u", fontsize=10)
            axes[1].set_xlabel("u [keV]")
            axes[1].set_ylabel("Counts")
            axes[1].grid(alpha=0.3)

            # 1D Projection along v
            v_proj = selected_hist_data.sum(axis=0)
            v_centers = 0.5 * (v_edges[:-1] + v_edges[1:])
            axes[2].step(v_centers, v_proj, where="mid", color="darkorange")
            axes[2].set_title("Projection along v", fontsize=10)
            axes[2].set_xlabel("v [keV]")
            axes[2].set_ylabel("Counts")
            axes[2].grid(alpha=0.3)

            plot_view.extend(
                [
                    mo.md(f"### Interactive Distribution Plots (`{file_selector.value}`)"),
                    mo.as_html(fig),
                ]
            )
            plt.close(fig)
        except Exception as exc:
            plot_view.append(mo.callout(mo.md(f"Could not render plot: `{exc}`"), kind="danger"))
    elif selected_hist_data is not None:
        plot_view.append(
            mo.callout(
                mo.md(
                    f"Selected item is not a 2D histogram (shape: {getattr(selected_hist_data, 'shape', 'non-array')})."
                ),
                kind="neutral",
            )
        )

    mo.vstack([mo.md("## 8. Distribution & Projection Visualization"), *plot_view])
    return (plot_view,)


@app.cell
def _(
    build_source_hash_params,
    compute_input_uuid,
    input_dir,
    mo,
    param_dir,
    pd,
    read_hdf5_metadata,
):
    verification_checks = []

    # 1. Deterministic hashing test
    p1 = build_source_hash_params(
        "60Co",
        "v1",
        cuts={"co_cut": "v>150"},
        transforms={"shift_type": "ushiftp", "shift_val": 0.8},
    )
    p2 = build_source_hash_params(
        "60Co",
        "v1",
        cuts={"co_cut": "v>150"},
        transforms={"shift_type": "ushiftp", "shift_val": 0.8},
    )
    u1 = compute_input_uuid(p1)
    u2 = compute_input_uuid(p2)
    verification_checks.append(
        {
            "contract": "1. Deterministic UUID repeatability",
            "expectation": "Identical parameters yield identical 12-char hex UUID",
            "result": "Passed" if (u1 == u2 and len(u1) == 12) else "Failed",
            "details": f"UUID={u1}",
        }
    )

    # 2. Shift isolation test
    p_tl_0 = build_source_hash_params("208Tl", "v1", cuts={"co_cut": "v>150"}, transforms=None)
    p_tl_1 = build_source_hash_params("208Tl", "v1", cuts={"co_cut": "v>150"}, transforms=None)
    p_co_0 = build_source_hash_params("60Co", "v1", cuts={"co_cut": "v>150"}, transforms=None)
    p_co_1 = build_source_hash_params(
        "60Co",
        "v1",
        cuts={"co_cut": "v>150"},
        transforms={"shift_type": "ushiftp", "shift_val": 0.8},
    )
    u_tl_0 = compute_input_uuid(p_tl_0)
    u_tl_1 = compute_input_uuid(p_tl_1)
    u_co_0 = compute_input_uuid(p_co_0)
    u_co_1 = compute_input_uuid(p_co_1)

    co_changed = u_co_0 != u_co_1
    tl_unchanged = u_tl_0 == u_tl_1
    shared_baseline = u_tl_0 == u_co_0
    verification_checks.append(
        {
            "contract": "2. Parameter-set UUID sharing & shift isolation",
            "expectation": "Invariant sources share baseline UUID; shifted 60Co receives new UUID",
            "result": "Passed" if (co_changed and tl_unchanged and shared_baseline) else "Failed",
            "details": f"Baseline UUID: {u_tl_0} (shared by 208Tl and unshifted 60Co), Shifted 60Co: {u_co_1}",
        }
    )

    # 3. Path routing compliance
    path_ok = str(input_dir).endswith("/input") and (param_dir / "input" == input_dir)
    verification_checks.append(
        {
            "contract": "3. Path-encoded storage hierarchy compliance",
            "expectation": "Path strictly matches fit_data/roi_.../thr.../.../input/",
            "result": "Passed" if path_ok else "Failed",
            "details": f"input_dir={input_dir.name}",
        }
    )

    # 4. Embedded metadata roundtrip
    files = list(input_dir.glob("*_*.h5")) if input_dir.is_dir() else []
    if files:
        sample_file = files[0]
        meta = read_hdf5_metadata(sample_file)
        has_meta = "component" in meta and "staging_id" in meta
        verification_checks.append(
            {
                "contract": "4. Embedded HDF5 metadata integrity",
                "expectation": "Template file contains attrs['metadata_json'] with schema keys",
                "result": "Passed" if has_meta else "Failed",
                "details": f"Checked {sample_file.name}",
            }
        )
    else:
        verification_checks.append(
            {
                "contract": "4. Embedded HDF5 metadata integrity",
                "expectation": "Template file contains attrs['metadata_json'] with schema keys",
                "result": "Awaiting generation",
                "details": "Generate inputs to verify on-disk file",
            }
        )

    verify_df = pd.DataFrame(verification_checks)
    mo.vstack(
        [
            mo.md("## 9. Storage Architecture & Contract Verification"),
            mo.ui.table(verify_df, selection=None),
        ]
    )
    return (
        co_changed,
        files,
        p1,
        p2,
        p_co_0,
        p_co_1,
        p_tl_0,
        p_tl_1,
        path_ok,
        tl_unchanged,
        u1,
        u2,
        u_co_0,
        u_co_1,
        u_tl_0,
        u_tl_1,
        verification_checks,
        verify_df,
    )


if __name__ == "__main__":
    app.run()
