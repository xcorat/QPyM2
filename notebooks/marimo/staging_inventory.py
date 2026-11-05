import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import hashlib
    import json
    from pathlib import Path

    import ROOT
    import marimo as mo
    import pandas as pd

    from qpym2.cfg.config_builders import (
        build_data_config,
        build_mc_config,
        build_ndbd_config,
        build_staging_config,
        get_artifact_store_dir,
        get_jags_path,
    )
    from qpym2.cfg.hydra_adapter import (
        HydraConfigManager,
        get_config,
        resolve_config_dir,
    )

    return (
        HydraConfigManager,
        Path,
        ROOT,
        build_data_config,
        build_mc_config,
        build_ndbd_config,
        build_staging_config,
        get_artifact_store_dir,
        get_config,
        get_jags_path,
        hashlib,
        json,
        mo,
        pd,
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
                # Config & Staging Inventory (v0.6)

                Inspect the staging tier and input datasets described by the QPyM2 configuration.
                Paths resolve against `infrastructure.data_root` (or `QPYM2_DATA_DIR`).
                Staging preselection is decoupled from downstream model ROI and fit options.
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
    options = HydraConfigManager.list_options()

    env_opts = options.get("environments", ["local"])
    environment = mo.ui.dropdown(
        options=env_opts,
        value="local" if "local" in env_opts else env_opts[0],
        label="Environment",
    )

    staging_opts = options.get("staging", ["default"])
    staging = mo.ui.dropdown(
        options=staging_opts,
        value="default" if "default" in staging_opts else staging_opts[0],
        label="Staging Tier",
    )

    analysis_opts = options.get("analysis", ["base"])
    analysis = mo.ui.dropdown(
        options=analysis_opts,
        value="base" if "base" in analysis_opts else analysis_opts[0],
        label="Analysis",
    )

    mo.hstack([environment, staging, analysis])
    return analysis, environment, staging


@app.cell
def _(analysis, config_dir_input, environment, get_config, mo, staging):
    cfg = get_config(
        environment=environment.value,
        staging=staging.value,
        analysis=analysis.value,
        config_dir=config_dir_input.value,
    )
    data_root = cfg.infrastructure.data_root
    mo.callout(
        mo.md(
            f"Loaded configuration for **`{environment.value}`** / **`{staging.value}`** / "
            f"**`{analysis.value}`** (data root: `{data_root}`)."
        ),
        kind="neutral",
    )
    return cfg, data_root


@app.cell
def _(Path, build_staging_config, cfg, data_root, mo, pd):
    staging_config = build_staging_config(cfg)
    staging_dir = Path(staging_config["outpath"])

    cuts_dict = cfg.staging.get("cuts", {})

    staging_summary = pd.DataFrame(
        [
            {"key": "data_root", "value": str(data_root)},
            {"key": "staging_dir", "value": str(staging_dir)},
            {"key": "output tree", "value": staging_config["treename"]},
            {"key": "output columns", "value": ", ".join(staging_config["out_cols"])},
            {"key": "broad_preselection", "value": str(cuts_dict.get("broad_preselection", "N/A"))},
            {"key": "combined staging cuts", "value": staging_config["cuts"]},
        ]
    )

    mo.vstack(
        [
            mo.md("## 1. Resolved staging configuration"),
            mo.md(
                "_Staging is model-independent: no ROI, threshold, or energy shifts are applied here._"
            ),
            mo.ui.table(staging_summary, selection=None),
        ]
    )
    return staging_config, staging_dir


@app.cell
def _(
    Path,
    build_data_config,
    build_mc_config,
    build_ndbd_config,
    cfg,
    get_jags_path,
    mo,
    pd,
):
    source_rows = []

    for kind, builder in [
        ("Background MC", build_mc_config),
        ("0vbb Signal MC", build_ndbd_config),
        ("Experimental Data", build_data_config),
    ]:
        source_config = builder(cfg)
        raw_inpath = source_config.get("inpath")
        inpath = Path(raw_inpath) if raw_inpath else Path(".")
        pattern = source_config.get("fname") or "*"
        matches = len(list(inpath.glob(pattern))) if inpath.is_dir() else 0
        source_rows.append(
            {
                "stage": kind,
                "input directory": str(inpath),
                "exists": inpath.is_dir(),
                "pattern": pattern,
                "matches": matches,
            }
        )

    jags_resolved = get_jags_path(cfg)
    jags_path = Path(jags_resolved) if jags_resolved else Path("Out.h5")
    source_rows.append(
        {
            "stage": "JAGS source table",
            "input directory": str(jags_path.parent),
            "exists": jags_path.is_file(),
            "pattern": jags_path.name,
            "matches": int(jags_path.is_file()),
        }
    )

    sources_df = pd.DataFrame(source_rows)
    mo.vstack(
        [
            mo.md("## 2. Upstream raw inputs"),
            mo.ui.table(sources_df, selection=None),
            mo.md("_Paths are resolved against `infrastructure.data_root`._"),
        ]
    )
    return


@app.cell
def _(mo, staging_dir):
    staging_dir_input = mo.ui.text(
        value=str(staging_dir),
        label="Staging directory to inspect",
        full_width=True,
    )
    scan_entries = mo.ui.checkbox(
        value=False,
        label="Read event counts from ROOT headers (all files)",
    )
    mo.vstack([staging_dir_input, scan_entries])
    return scan_entries, staging_dir_input


@app.cell
def _(Path, ROOT, mo, pd, scan_entries, staging_dir_input):
    scan_dir = Path(staging_dir_input.value).expanduser()
    staging_files = (
        sorted(p for p in scan_dir.glob("*.root") if p.is_file()) if scan_dir.is_dir() else []
    )

    def _classify(fname: str) -> str:
        if fname in ("data.root", "Reduced_1ton.root"):
            return "Experimental Data"
        if fname == "data_blinded.root":
            return "Experimental Data (Blinded)"
        if fname.startswith("ndbd"):
            return "0vbb Signal MC"
        if fname.startswith("mu-"):
            return "Muon MC"
        return "Background MC"

    if staging_files:
        rows = []
        entries_dict = {}
        if scan_entries.value:
            for p in staging_files:
                try:
                    rf = ROOT.TFile.Open(str(p))
                    if rf and not rf.IsZombie():
                        t = rf.Get("uvtree")
                        entries_dict[p.name] = t.GetEntries() if t else -1
                        rf.Close()
                    else:
                        entries_dict[p.name] = -1
                except Exception:
                    entries_dict[p.name] = -1

        for path in staging_files:
            row = {
                "file": path.name,
                "category": _classify(path.name),
                "size_mb": round(path.stat().st_size / 1e6, 3),
            }
            if scan_entries.value:
                row["entries"] = entries_dict.get(path.name, -1)
            row["mtime"] = pd.Timestamp(path.stat().st_mtime, unit="s")
            rows.append(row)

        staging_inventory = pd.DataFrame(rows)
        total_mb = round(staging_inventory["size_mb"].sum(), 3)

        cat_agg = {"files": ("file", "count"), "total_mb": ("size_mb", "sum")}
        if scan_entries.value and "entries" in staging_inventory.columns:
            cat_agg["total_entries"] = ("entries", "sum")
            total_events = staging_inventory["entries"].sum()
            callout_text = f"Found **{len(staging_files)} files** ({total_mb} MB, **{total_events:,} events**) in `{scan_dir}`."
        else:
            callout_text = f"Found **{len(staging_files)} files** ({total_mb} MB) in `{scan_dir}`."

        cat_summary = staging_inventory.groupby("category").agg(**cat_agg).reset_index()
        cat_summary["total_mb"] = cat_summary["total_mb"].round(2)

        inventory_view = mo.ui.table(staging_inventory, selection=None)
        cat_view = mo.ui.table(cat_summary, selection=None)
        status = mo.callout(mo.md(callout_text), kind="success")
        content = mo.vstack(
            [
                mo.md("## 3. Staging inventory"),
                status,
                mo.md("### Category breakdown"),
                cat_view,
                mo.md("### All staged files"),
                inventory_view,
            ]
        )
    elif scan_dir.is_dir():
        content = mo.vstack(
            [
                mo.md("## 3. Staging inventory"),
                mo.callout(
                    mo.md(f"Directory exists but contains no ROOT files: `{scan_dir}`."),
                    kind="warn",
                ),
            ]
        )
    else:
        content = mo.vstack(
            [
                mo.md("## 3. Staging inventory"),
                mo.callout(
                    mo.md(
                        f"Staging directory not reachable: `{scan_dir}`. "
                        "Edit the path above or select an accessible environment."
                    ),
                    kind="warn",
                ),
            ]
        )

    content
    return scan_dir, staging_files


@app.cell
def _(mo, staging_files):
    read_jags_sources = mo.ui.switch(
        value=False,
        label="Check expected MC outputs from JAGS source table",
    )
    file_options = [f.name for f in staging_files] if staging_files else []
    file_to_inspect = mo.ui.dropdown(
        options=file_options,
        value=file_options[0] if file_options else None,
        label="Select ROOT file to inspect",
    )
    inspect_root = mo.ui.switch(
        value=False,
        label="Open selected ROOT file",
    )
    mo.hstack([read_jags_sources, file_to_inspect, inspect_root])
    return file_to_inspect, inspect_root, read_jags_sources


@app.cell
def _(
    Path,
    build_data_config,
    build_ndbd_config,
    cfg,
    get_jags_path,
    mo,
    pd,
    read_jags_sources,
    scan_dir,
):
    expected_names = {
        build_data_config(cfg)["out_fname"],
        build_ndbd_config(cfg)["out_fname"],
    }
    source_note = "config-only expectations"

    if read_jags_sources.value:
        jags_p = get_jags_path(cfg)
        if jags_p and Path(jags_p).is_file():
            try:
                from qpym2.io.bkg_model import read_jagsh5

                source_table, _ = read_jagsh5(jags_p)
                expected_names.update(source_table["fname"].tolist())
                source_note = "config + JAGS source table"
            except Exception as exc:
                source_note = f"JAGS source table read error ({exc})"
        else:
            source_note = f"JAGS source table not reachable ({jags_p})"

    expected_df = pd.DataFrame(
        [
            {
                "expected file": name,
                "present": (Path(scan_dir) / name).is_file(),
                "path": str(Path(scan_dir) / name),
            }
            for name in sorted(expected_names)
        ]
    )
    mo.vstack(
        [
            mo.md(f"## 4. Completeness check ({source_note})"),
            mo.ui.table(expected_df, selection=None),
        ]
    )
    return


@app.cell
def _(ROOT, file_to_inspect, inspect_root, mo, scan_dir, staging_config):
    if not inspect_root.value:
        root_view = mo.md("_ROOT inspection disabled._")
    elif not file_to_inspect.value:
        root_view = mo.callout(mo.md("No staging file selected."), kind="warn")
    else:
        target_path = scan_dir / file_to_inspect.value
        try:
            root_file = ROOT.TFile.Open(str(target_path))
            treename = staging_config["treename"]
            tree = root_file.Get(treename) if root_file and not root_file.IsZombie() else None
            if tree:
                branches = [branch.GetName() for branch in tree.GetListOfBranches()]
                root_view = mo.vstack(
                    [
                        mo.md(
                            f"**File:** `{file_to_inspect.value}` -> **Tree:** `{treename}` "
                            f"(**{tree.GetEntries():,}** entries)"
                        ),
                        mo.md("Branches: " + ", ".join(f"`{name}`" for name in branches)),
                    ]
                )
            else:
                root_view = mo.callout(
                    mo.md(f"Tree `{treename}` not found in `{file_to_inspect.value}`."),
                    kind="warn",
                )
            if root_file:
                root_file.Close()
        except Exception as exc:
            root_view = mo.callout(mo.md(f"ROOT read failed: `{exc}`"), kind="danger")

    mo.vstack([mo.md("## 5. ROOT tree check"), root_view])
    return


if __name__ == "__main__":
    app.run()
