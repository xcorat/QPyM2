import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import pandas as pd

    from qpym2.cfg import get_config
    from qpym2.cfg.config_builders import (
        build_input_config,
        get_artifact_store_dir,
        get_comp_cfg,
        get_data_names,
        get_jags_path,
        get_h2uv_model,
        get_staging_dir,
    )
    from qpym2.cfg.hydra_adapter import HydraConfigManager, resolve_config_dir
    from qpym2.fit.create_m2uv import build_input
    from qpym2.fit.params import model_cfg
    from qpym2.io.bkg_model import read_jagsh5
    from qpym2.io.fit_data.input_set import InputSet
    from qpym2.io.fit_data.objects import Observation, PriorInput, Template
    from qpym2.storage.fit_data_store import FitDataStore
    from qpym2.storage.hasher import hash_file
    from qpym2.storage.identity_profiles import encode_id

    return (
        FitDataStore,
        HydraConfigManager,
        InputSet,
        Observation,
        Path,
        PriorInput,
        Template,
        build_input,
        build_input_config,
        encode_id,
        get_artifact_store_dir,
        get_comp_cfg,
        get_config,
        get_data_names,
        get_jags_path,
        get_h2uv_model,
        get_staging_dir,
        hash_file,
        model_cfg,
        mo,
        np,
        pd,
        read_jagsh5,
        resolve_config_dir,
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
                # Fit Data Layer: Input Storage & Validation

                Interactive workbench to exercise and validate the new Fit Data Layer
                (`FitDataStore`, `Observation`, `Template`, `InputSet`, and sidecar-backed identities)
                using the real CUORE M2 configuration and physics parameters.
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
    _options = HydraConfigManager.list_options()

    _env_opts = _options.get("environments", ["local"])
    env_dropdown = mo.ui.dropdown(
        options=_env_opts,
        value="local" if "local" in _env_opts else _env_opts[0],
        label="Environment",
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

    mo.hstack([env_dropdown, analysis_dropdown, model_config_dropdown])
    return analysis_dropdown, env_dropdown, model_config_dropdown


@app.cell
def _(
    Path,
    analysis_dropdown,
    config_dir_input,
    env_dropdown,
    get_artifact_store_dir,
    get_config,
    get_staging_dir,
    mo,
    model_config_dropdown,
):
    cfg = get_config(
        environment=env_dropdown.value,
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

    _default_store_dir = Path(get_artifact_store_dir(cfg))
    store_dir_input = mo.ui.text(
        value=str(_default_store_dir),
        label="FitDataStore root directory",
        full_width=True,
    )

    staging_path = Path(get_staging_dir(cfg))

    mo.vstack([mo.md("## Model Selection & Store Root"), mo.hstack([model_key, store_dir_input])])
    return cfg, model_key, staging_path, store_dir_input


@app.cell
def _(
    FitDataStore,
    Path,
    cfg,
    get_data_names,
    get_comp_cfg,
    get_jags_path,
    get_h2uv_model,
    build_input_config,
    mo,
    model_key,
    pd,
    read_jagsh5,
    staging_path,
    store_dir_input,
):
    fit_store = FitDataStore(store_dir_input.value)
    store_path = Path(store_dir_input.value)

    _folders = ["observations", "templates", "inputs", "fits", ".tmp"]
    _folder_status = []
    for _folder in _folders:
        _p = store_path / _folder
        _count = len(list(_p.rglob("*.json"))) if _p.is_dir() else 0
        _folder_status.append(
            {
                "directory": _folder,
                "exists": _p.is_dir(),
                "JSON sidecars": _count,
            }
        )

    store_status_df = pd.DataFrame(_folder_status)

    h2uv_model = get_h2uv_model(cfg, model_key.value)
    comp_configs = get_comp_cfg(cfg, model_key.value)
    comp_names = [c["name"] for c in comp_configs]

    _co_cut_opts = ["None", *[str(c) for c in cfg.cuts.co_cuts.keys()]]
    co_cut_input = mo.ui.dropdown(
        options=_co_cut_opts,
        value="gamma10_v150" if "gamma10_v150" in _co_cut_opts else _co_cut_opts[0],
        label="Co-60 cut",
    )
    _co_shift_opts = [0.0, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0]
    co_shift_input = mo.ui.dropdown(
        options=_co_shift_opts,
        value=1.0,
        label="Co-60 energy shift (keV)",
    )

    _ndbd_cut_opts = ["None", *[str(c) for c in cfg.cuts.ndbd_cutoff.keys()]]
    ndbd_cutoff_input = mo.ui.dropdown(
        options=_ndbd_cut_opts,
        value="u1740" if "u1740" in _ndbd_cut_opts else _ndbd_cut_opts[0],
        label="0νββ cutoff",
    )
    ndbd_shift_input = mo.ui.dropdown(
        options=[0.0, 0.3, 0.5, 0.8, 1.0],
        value=0.0,
        label="0νββ energy shift (keV)",
    )

    skip_roi_input = mo.ui.dropdown(
        options=["None", "2515-2535 (blinding window)"],
        value="None",
        label="Skip ROI (2515-2535 keV)",
    )

    dataset_input = mo.ui.text(value="1ton", label="Dataset name")
    blinded_input = mo.ui.checkbox(value=False, label="Blinded data")

    input_cfg = build_input_config(cfg, model_key.value)
    jags_path = Path(get_jags_path(cfg))
    data_names = get_data_names(cfg)
    try:
        source_table, source_chain = read_jagsh5(str(jags_path))
        source_error = None
    except Exception as exc:
        source_table = pd.DataFrame(columns=["fname", "mean", "std", "mode"])
        source_chain = None
        source_error = str(exc)

    data_paths = {
        "unblinded": staging_path / f"{data_names[0]}.root",
        "blinded": staging_path / f"{data_names[1]}.root",
    }
    source_status = [
        {"input": "JAGS background model", "path": str(jags_path), "available": jags_path.is_file()},
        *[
            {"input": f"{label} staged data", "path": str(path), "available": path.is_file()}
            for label, path in data_paths.items()
        ],
        {
            "input": "NDBD staged MC",
            "path": str(staging_path / "ndbd.root"),
            "available": (staging_path / "ndbd.root").is_file(),
        },
    ]
    source_status.extend(
        {
            "input": f"Background MC: {source}",
            "path": str(staging_path / str(row["fname"])),
            "available": (staging_path / str(row["fname"])).is_file(),
        }
        for source, row in source_table.iterrows()
    )
    source_status_df = pd.DataFrame(source_status)

    mo.vstack(
        [
            mo.md("### Store Layout Status"),
            mo.ui.table(store_status_df, selection=None),
            mo.md("### Required Real Inputs"),
            mo.ui.table(source_status_df, selection=None),
            mo.md(
                f"JAGS source table: `{source_error or f'{len(source_table)} MC sources loaded'}`"
            ),
            mo.md("## Resolved Physics Parameters"),
            mo.md(
                f"- **Model:** `{model_key.value}` (`{h2uv_model.name}`)  \n"
                f"- **Energy ROI:** `esum_min = {h2uv_model.esum_min:g} keV`, `esum_max = {h2uv_model.esum_max:g} keV`  \n"
                f"- **Threshold:** `{h2uv_model.threshold:g} keV`, **Bin Size:** `{h2uv_model.bin_size:g} keV`  \n"
                f"- **Active fit components:** {', '.join(f'`{c}`' for c in comp_names)}  \n"
                f"- **Source-level MC templates:** {len(source_table)} background sources + `ndbd`  \n"
                f"- **Staging Path:** `{staging_path}`"
            ),
            mo.md("### Cut & Blinding Controls (ROI/binning come from the selected model)"),
            mo.hstack([dataset_input, blinded_input]),
            mo.hstack([skip_roi_input]),
            mo.hstack([co_cut_input, co_shift_input]),
            mo.hstack([ndbd_cutoff_input, ndbd_shift_input]),
        ]
    )
    return (
        blinded_input,
        co_cut_input,
        co_shift_input,
        comp_names,
        data_names,
        data_paths,
        dataset_input,
        fit_store,
        h2uv_model,
        input_cfg,
        jags_path,
        ndbd_cutoff_input,
        ndbd_shift_input,
        source_error,
        source_status_df,
        source_chain,
        source_table,
        skip_roi_input,
    )


@app.cell
def _(
    blinded_input,
    co_cut_input,
    co_shift_input,
    dataset_input,
    h2uv_model,
    mo,
    ndbd_cutoff_input,
    ndbd_shift_input,
    skip_roi_input,
):
    _sel_co_cut = None if co_cut_input.value == "None" else co_cut_input.value
    _sel_ndbd_cutoff = None if ndbd_cutoff_input.value == "None" else ndbd_cutoff_input.value
    _sel_skip_roi = [2515.0, 2535.0] if "2515" in skip_roi_input.value else None
    binning_spec = {
        "bin_size": float(h2uv_model.bin_size),
        "esum_min": float(h2uv_model.esum_min),
        "esum_max": float(h2uv_model.esum_max),
        "threshold": float(h2uv_model.threshold),
    }
    if _sel_skip_roi is not None:
        binning_spec["skip_roi"] = _sel_skip_roi

    mo.vstack(
        [
            mo.md("## Real Input Selection"),
            mo.md(
                f"**Dataset label:** `{dataset_input.value}`  \n"
                f"**Observed sample:** `{'blinded' if blinded_input.value else 'unblinded'}` staging ROOT  \n"
                f"**Co-60 cut/shift:** `{_sel_co_cut or 'none'}` / `{co_shift_input.value:g} keV`  \n"
                f"**NDBD cutoff/shift:** `{_sel_ndbd_cutoff or 'none'}` / `{ndbd_shift_input.value:g} keV`  \n"
                f"**Skip ROI:** `{_sel_skip_roi or 'none'}`  \n"
                f"**Config-derived binning:** `{binning_spec}`"
            ),
            mo.md(
                "No synthetic fallback is used. Artifact IDs are derived only after the staged "
                "data, JAGS table, and all staged MC inputs have been read successfully."
            ),
        ]
    )
    return binning_spec

@app.cell
def _(mo):
    publish_button = mo.ui.run_button(label="Build from staging ROOT")
    last_publish_click, set_last_publish_click = mo.state(0)
    mo.vstack(
        [
            mo.md("## Storage Action"),
            mo.md(
                "Build real observed data and source-level MC templates using the configured "
                "staging ROOT files, then publish validated objects and JSON sidecars."
            ),
            publish_button,
        ]
    )
    return last_publish_click, publish_button, set_last_publish_click


@app.cell
def _(
    InputSet,
    Observation,
    PriorInput,
    Template,
    binning_spec,
    build_input,
    cfg,
    co_cut_input,
    co_shift_input,
    data_paths,
    dataset_input,
    encode_id,
    blinded_input,
    fit_store,
    hash_file,
    h2uv_model,
    input_cfg,
    jags_path,
    last_publish_click,
    model_cfg,
    mo,
    ndbd_cutoff_input,
    ndbd_shift_input,
    np,
    pd,
    publish_button,
    set_last_publish_click,
    skip_roi_input,
    source_error,
    source_chain,
    source_table,
    staging_path,
):
    publish_status = mo.md("_Choose ‘Build from staging ROOT’ to create real fit inputs._")
    published_files = []
    published_input_set_id = None
    published_input_set_params = None
    published_template_ids = {}

    if publish_button.value and publish_button.value > last_publish_click():
        set_last_publish_click(publish_button.value)
        try:
            if source_error:
                raise RuntimeError(f"Could not read the configured JAGS table: {source_error}")
            if source_chain is None:
                raise RuntimeError("The configured JAGS source chain is unavailable")

            _selected_blinding = "blinded" if bool(blinded_input.value) else "unblinded"
            _data_path = data_paths[_selected_blinding]
            _ndbd_path = staging_path / "ndbd.root"
            _source_paths = {
                str(source): staging_path / str(row["fname"])
                for source, row in source_table.iterrows()
            }
            _required_paths = [jags_path, _data_path, _ndbd_path, *_source_paths.values()]
            _missing = [path for path in _required_paths if not path.is_file()]
            if _missing:
                raise FileNotFoundError(
                    "Required staged inputs are missing: "
                    + ", ".join(str(path) for path in _missing)
                )

            _sel_co_cut = None if co_cut_input.value == "None" else co_cut_input.value
            _sel_ndbd_cutoff = (
                None if ndbd_cutoff_input.value == "None" else ndbd_cutoff_input.value
            )
            _sel_skip_roi = (
                (2515.0, 2535.0) if "2515" in skip_roi_input.value else False
            )
            import copy

            _hist_model = copy.deepcopy(h2uv_model)
            _hist_model.skip_roi = None
            _selected_input_cfg = dict(input_cfg)
            _selected_input_cfg["data_fname"] = _data_path.name
            _fit_cfg = model_cfg(
                _hist_model,
                co_cut=_sel_co_cut,
                ndbd_cutoff=_sel_ndbd_cutoff,
                co_shift=float(co_shift_input.value),
                ndbd_shift=float(ndbd_shift_input.value),
                skip_roi=_sel_skip_roi,
                blinded=False,
                smooth=True,
                common_cuts=[],
                co_cuts_dict=dict(cfg.cuts.co_cuts),
                ndbd_cuts_dict=dict(cfg.cuts.ndbd_cutoff),
            )

            _mctable, _signal_df, _data, _b_data, _edges = build_input(
                _fit_cfg, _selected_input_cfg, str(jags_path), hist_cache={}
            )
            if _data is None or len(_data) == 0:
                raise ValueError("The configured staged data selection produced no events")
            if _mctable.empty or _signal_df.empty:
                raise ValueError("The configured staging inputs produced empty MC templates")

            _data_source_id = hash_file(_data_path)
            _jags_source_id = hash_file(jags_path)
            _prior_input = PriorInput(
                data=source_chain,
                source="JAGS background model",
                kind="mkchain",
                source_id=_jags_source_id,
                metadata={"source_file": str(jags_path)},
            )
            _prior_input_id = _prior_input.derive_id()
            _obs_params = {
                "dataset": dataset_input.value,
                "selection": f"co_cut={_sel_co_cut or 'none'}",
                "blinded": bool(blinded_input.value),
                "coords": ["u", "v"],
                "esum_min": float(h2uv_model.esum_min),
                "esum_max": float(h2uv_model.esum_max),
                "threshold": float(h2uv_model.threshold),
                "source_id": _data_source_id,
            }
            if _sel_skip_roi:
                _obs_params["skip_roi"] = list(_sel_skip_roi)
            _obs = Observation(
                data=np.asarray(_data),
                dataset=dataset_input.value,
                selection=_obs_params["selection"],
                blinded=bool(blinded_input.value),
                skip_roi=_sel_skip_roi or None,
                esum_min=float(h2uv_model.esum_min),
                esum_max=float(h2uv_model.esum_max),
                threshold=float(h2uv_model.threshold),
                source_id=_data_source_id,
                coords=("u", "v"),
                metadata={"source_file": str(_data_path)},
            )
            _obs_id = _obs.derive_id()

            _templates = []
            _template_ids = {}
            _bkg_rows = []
            for _source, _row in _mctable.iterrows():
                _fname = str(_row["fname"])
                _source_path = _source_paths[str(_source)]
                _source_id = hash_file(_source_path)
                _is_co60 = "co60" in _fname.lower()
                _template = Template(
                    hist=np.asarray(_row["mchist"]),
                    source=str(_source),
                    role="background",
                    edges=_edges,
                    bin_size=float(h2uv_model.bin_size),
                    esum_min=float(h2uv_model.esum_min),
                    esum_max=float(h2uv_model.esum_max),
                    threshold=float(h2uv_model.threshold),
                    co_cut=_sel_co_cut,
                    skip_roi=_sel_skip_roi or None,
                    co_shift=float(co_shift_input.value) if _is_co60 else 0.0,
                    source_id=_source_id,
                    metadata={"source_file": _fname},
                )
                _tid = _template.derive_id()
                _templates.append((_template, _tid))
                _template_ids[str(_source)] = _tid
                _bkg_rows.append(
                    {"source": str(_source), "mean": _row["mean"], "mode": _row["mode"]}
                )

            _signal_source_id = hash_file(_ndbd_path)
            _signal_row = _signal_df.iloc[0]
            _signal_template = Template(
                hist=np.asarray(_signal_row["mchist"]),
                source="ndbd",
                role="signal",
                edges=_edges,
                bin_size=float(h2uv_model.bin_size),
                esum_min=float(h2uv_model.esum_min),
                esum_max=float(h2uv_model.esum_max),
                threshold=float(h2uv_model.threshold),
                co_cut=_sel_co_cut,
                ndbd_cutoff=_sel_ndbd_cutoff,
                skip_roi=_sel_skip_roi or None,
                ndbd_shift=float(ndbd_shift_input.value),
                source_id=_signal_source_id,
                metadata={"source_file": _ndbd_path.name},
            )
            _signal_template_id = _signal_template.derive_id()

            _input_set = InputSet(
                observation_id=_obs_id,
                template_ids=_template_ids,
                signal_template_id=_signal_template_id,
                binning=binning_spec,
                edges=_edges,
                bkg_model_params=pd.DataFrame(_bkg_rows).set_index("source"),
                signal_params=_signal_df[["mean", "mode"]].copy(),
                prior_input_id=_prior_input_id,
                metadata={"jags_file": str(jags_path)},
            )
            _input_set_params = _input_set.get_identity_params()
            _input_set_id = encode_id("input_set", _input_set_params)

            def _publish_or_reuse(_kind, _object, _object_id):
                if fit_store.exists(_kind, _object_id):
                    if not fit_store.verify(_kind, _object_id):
                        raise ValueError(f"Existing {_kind} artifact failed validation: {_object_id}")
                    return fit_store.resolve_path(_kind, _object_id), True
                return _object.publish(fit_store, object_id=_object_id), False

            _objects = [("observation", _obs, _obs_id)]
            _objects.extend(("template", obj, oid) for obj, oid in _templates)
            _objects.append(("template", _signal_template, _signal_template_id))
            _objects.append(("prior_input", _prior_input, _prior_input_id))
            _objects.append(("input_set", _input_set, _input_set_id))
            for _kind, _object, _object_id in _objects:
                _path, _reused = _publish_or_reuse(_kind, _object, _object_id)
                _sidecar = fit_store.sidecar_path(_kind, _object_id)
                published_files.append(
                    {
                        "kind": _kind,
                        "id": _object_id,
                        "payload": str(_path),
                        "sidecar": str(_sidecar),
                        "state": "reused" if _reused else "published",
                    }
                )

            published_input_set_id = _input_set_id
            published_input_set_params = _input_set_params
            published_template_ids = _template_ids | {"ndbd": _signal_template_id}
            publish_status = mo.callout(
                mo.md(
                    f"Built inputs from real staging files. Observation events: {len(_data):,}; "
                    f"source MC templates: {len(_templates)} background + 1 signal."
                ),
                kind="success",
            )
        except Exception as exc:
            publish_status = mo.callout(
                mo.md(f"Input build failed; no synthetic fallback was used.\n\n`{exc}`"),
                kind="danger",
            )

    _published_table = pd.DataFrame(published_files)
    mo.vstack(
        [
            publish_status,
            mo.ui.table(_published_table, selection=None) if published_files else mo.md(""),
        ]
    )
    return (
        published_input_set_id,
        published_input_set_params,
        published_template_ids,
    )

@app.cell
def _(InputSet, fit_store, mo, np, pd, published_input_set_id):
    if published_input_set_id and fit_store.exists("input_set", published_input_set_id):
        _loaded_set = InputSet.load(fit_store, published_input_set_id)
        _manifest = fit_store.read_manifest("input_set", published_input_set_id)
        _mctable, _signal_df, _resolved_data, _edges = _loaded_set.resolve(fit_store)

        _source_rows = []
        for _source, _row in _mctable.iterrows():
            _source_rows.append(
                {
                    "MC source": _source,
                    "prior mean": _row.get("mean"),
                    "prior mode": _row.get("mode"),
                    "histogram sum": float(np.asarray(_row["mchist"]).sum()),
                }
            )
        _source_table = pd.DataFrame(_source_rows)
        _signal_sum = float(np.asarray(_signal_df.iloc[0]["mchist"]).sum())
        _edge_str = ", ".join(f"{len(axis) - 1} bins" for axis in _edges) if _edges else "None"
        _partition = _manifest.get("partition", {})

        inspection_view = mo.vstack(
            [
                mo.md(f"## Real InputSet (`{published_input_set_id}`)"),
                mo.md(
                    f"**Observed events:** {len(_resolved_data):,}  \
"
                    f"**Background MC sources:** {len(_mctable)}  \
"
                    f"**Signal histogram sum:** {_signal_sum:,.1f}  \
"
                    f"**Histogram geometry:** {_edge_str}  \
"
                    f"**Partition:** `{_partition}`  \
"
                    f"**Sidecar kind:** `{_manifest.get('kind')}`"
                ),
                mo.ui.table(_source_table, selection=None),
            ]
        )
    else:
        inspection_view = mo.callout(
            mo.md("No InputSet was published. Check the real staging inputs and build status above."),
            kind="info",
        )

    inspection_view
    return


@app.cell
def _(
    co_shift_input,
    fit_store,
    mo,
    ndbd_shift_input,
    pd,
    published_input_set_id,
    published_input_set_params,
    published_template_ids,
):
    _validation_items = []
    if published_input_set_id and published_input_set_params:
        try:
            with fit_store.stage_object(
                "input_set",
                published_input_set_id,
                overwrite=False,
                identity=published_input_set_params,
            ) as _tmp:
                _tmp.write_text("this duplicate publication must not be reached")
            _immutability = "Failed"
        except FileExistsError:
            _immutability = "Passed"
        _validation_items.append(
            {
                "contract": "Immutable publication",
                "result": _immutability,
            }
        )

        _co_shift_ok = True
        _ndbd_shift_ok = True
        for _source, _template_id in published_template_ids.items():
            _manifest = fit_store.read_manifest("template", _template_id)
            _identity = _manifest["identity"]
            _source_file = str(_manifest.get("metadata", {}).get("source_file", ""))
            if _source == "ndbd":
                _ndbd_shift_ok &= float(_identity["ndbd_shift"]) == float(ndbd_shift_input.value)
            elif "co60" in _source_file.lower():
                _co_shift_ok &= float(_identity["co_shift"]) == float(co_shift_input.value)
            else:
                _co_shift_ok &= float(_identity["co_shift"]) == 0.0

        _validation_items.extend(
            [
                {"contract": "Only Co-60 source receives Co shift", "result": "Passed" if _co_shift_ok else "Failed"},
                {"contract": "NDBD source receives configured NDBD shift", "result": "Passed" if _ndbd_shift_ok else "Failed"},
            ]
        )
    else:
        _validation_items.append(
            {"contract": "Real inputs available", "result": "Build required before validation"}
        )

    _val_df = pd.DataFrame(_validation_items)
    mo.vstack([mo.md("## Contract & Architecture Validation"), mo.ui.table(_val_df, selection=None)])
    return


if __name__ == "__main__":
    app.run()
