"""
Unit tests for qpym2.storage (path routing, schema, hashing, metadata, backend, study).
"""

import numpy as np
import pandas as pd

from qpym2.storage.schema import (
    normalize_cut_string,
    normalize_float,
    build_template_hash_dict,
    build_fit_hash_dict,
)
from qpym2.storage.hasher import (
    canonical_json,
    hash_dict,
    compute_input_uuid,
)
from qpym2.storage.paths import (
    format_roi_dir,
    format_thr_dir,
    format_coord_bin_dir,
    get_parameter_dir,
    get_input_dir,
    get_output_dir,
    get_template_path,
    get_data_path,
    get_fit_path,
    get_study_dir,
    get_study_paths,
)
from qpym2.storage.backend import (
    write_template,
    read_template,
    write_observed_data,
    read_observed_data,
    write_fit,
    read_fit,
    template_exists,
    data_exists,
    fit_exists,
)
from qpym2.storage.study import (
    init_study_db,
    record_study_point,
    write_study_manifest,
    read_study_manifest,
    collect_study_results,
)


def test_path_routing(tmp_path):
    assert format_roi_dir(2400.0, 2570.0) == "roi_2400_2570"
    assert format_thr_dir(200) == "thr200"
    assert format_thr_dir(200.0) == "thr200"
    assert format_thr_dir(150.5) == "thr150.5"
    assert format_coord_bin_dir("h2uv", 2.5) == "h2uv_b2.5"
    assert format_coord_bin_dir("h2uv", 1.0) == "h2uv_b1"

    param_dir = get_parameter_dir(tmp_path, (2400, 2570), 200, "h2uv", 2.5)
    expected = tmp_path / "roi_2400_2570" / "thr200" / "h2uv_b2.5"
    assert param_dir == expected

    inp_dir = get_input_dir(param_dir)
    out_dir = get_output_dir(param_dir)
    assert inp_dir == expected / "input"
    assert out_dir == expected / "output"

    tpl_path = get_template_path(inp_dir, "co60", "a3f81c9d04e2")
    assert tpl_path == inp_dir / "co60_a3f81c9d04e2.h5"

    data_path = get_data_path(inp_dir, blinded=False)
    assert data_path == inp_dir / "data.h5"

    bdata_path = get_data_path(inp_dir, blinded=True)
    assert bdata_path == inp_dir / "data_blinded.h5"

    fit_h5 = get_fit_path(out_dir, "b7e291c04f81", ext=".h5")
    assert fit_h5 == out_dir / "fit_b7e291c04f81.h5"

    fit_nc = get_fit_path(out_dir, "b7e291c04f81", ext=".nc")
    assert fit_nc == out_dir / "fit_b7e291c04f81.nc"

    study_dir = get_study_dir(tmp_path, "sensitivity_v1")
    assert study_dir == tmp_path / "studies" / "sensitivity_v1"

    spaths = get_study_paths(study_dir, "sensitivity_v1")
    assert spaths["manifest"] == study_dir / "sensitivity_v1.manifest"
    assert spaths["sqlite"] == study_dir / "sensitivity_v1.sqlite"
    assert spaths["results_csv"] == study_dir / "results.csv"
    assert spaths["plots_dir"] == study_dir / "plots"


def test_schema_and_normalization():
    assert normalize_cut_string("   ") is None
    assert normalize_cut_string("true") is None
    assert normalize_cut_string("  v>150 && true && e1 > 1000  ") == "v>150 && e1 > 1000"

    assert normalize_float(0.800000000001) == 0.8
    assert normalize_float(None) is None

    tpl_dict = build_template_hash_dict(
        component="co60",
        staging_id="v5.1",
        cuts={"common": "  v>150  "},
        transforms={"shift_val": 0.80000001},
        grouping={"filter_regex": "60Co"},
        norm="mode",
    )
    assert tpl_dict["component"] == "co60"
    assert tpl_dict["cuts"]["common"] == "v>150"
    assert tpl_dict["transforms"]["shift_val"] == 0.8

    fit_dict = build_fit_hash_dict(
        model_name="m2",
        data_source="data",
        templates={"u238": "uuid_u", "co60": "uuid_co"},
        smoothing={"enabled": True},
        priors={"sigma_scale": 1.0},
        roi_settings={"skip_roi": [2515.0, 2540.0]},
        mcmc={"chains": 8, "draws": 20000},
    )
    # Check templates are sorted by key
    assert list(fit_dict["templates"].keys()) == ["co60", "u238"]
    assert fit_dict["roi_settings"]["skip_roi"] == [2515.0, 2540.0]


def test_deterministic_hashing():
    d1 = {"b": 2, "a": 1, "nested": {"y": 20, "x": 10}}
    d2 = {"a": 1, "nested": {"x": 10, "y": 20}, "b": 2}
    # Key ordering in dict must produce identical JSON and hash
    assert canonical_json(d1) == canonical_json(d2)
    assert hash_dict(d1) == hash_dict(d2)
    assert len(hash_dict(d1, length=12)) == 12

    # Template hashing deduplication check:
    # co60 with shift 0.0 vs shift 0.8
    co_shift0 = build_template_hash_dict(
        component="co60", staging_id="v5.1",
        cuts={"common": "v>150"}, transforms={"shift_val": 0.0}
    )
    co_shift8 = build_template_hash_dict(
        component="co60", staging_id="v5.1",
        cuts={"common": "v>150"}, transforms={"shift_val": 0.8}
    )
    # u238 has no shift transform
    u_shift0 = build_template_hash_dict(
        component="u238", staging_id="v5.1",
        cuts={"common": "v>150"}, transforms=None
    )
    u_shift8 = build_template_hash_dict(
        component="u238", staging_id="v5.1",
        cuts={"common": "v>150"}, transforms=None
    )

    uuid_co0 = compute_input_uuid(co_shift0)
    uuid_co8 = compute_input_uuid(co_shift8)
    uuid_u0 = compute_input_uuid(u_shift0)
    uuid_u8 = compute_input_uuid(u_shift8)

    # Co60 UUID must change when shift changes
    assert uuid_co0 != uuid_co8
    # U238 UUID must be identical across runs
    assert uuid_u0 == uuid_u8


def test_template_roundtrip(tmp_path):
    tpl_file = tmp_path / "co60_123456789abc.h5"
    assert not template_exists(tpl_file)

    hist_arr = np.arange(12, dtype=np.float64).reshape((3, 4))
    edges = [np.array([0, 1, 2, 3]), np.array([10, 20, 30, 40, 50])]
    meta = {"component": "co60", "shift": 0.8, "provenance": "test"}

    write_template(tpl_file, hist_arr, edges=edges, metadata=meta)
    assert template_exists(tpl_file)

    r_hist, r_edges, r_meta = read_template(tpl_file)
    np.testing.assert_array_equal(r_hist, hist_arr)
    assert len(r_edges) == 2
    np.testing.assert_array_equal(r_edges[0], edges[0])
    np.testing.assert_array_equal(r_edges[1], edges[1])
    assert r_meta["component"] == "co60"
    assert r_meta["shift"] == 0.8


def test_observed_data_roundtrip(tmp_path):
    data_file = tmp_path / "data.h5"
    assert not data_exists(data_file)

    data_arr = np.array([[100.5, 200.5], [101.2, 202.1]], dtype=np.float64)
    meta = {"dataset": "nominal", "events": 2}

    write_observed_data(data_file, data_arr, metadata=meta)
    assert data_exists(data_file)

    r_data, r_edges, r_meta = read_observed_data(data_file)
    np.testing.assert_array_equal(r_data, data_arr)
    assert r_meta["dataset"] == "nominal"
    assert r_meta["events"] == 2


def test_fit_roundtrip(tmp_path):
    fit_prefix = tmp_path / "fit_a1b2c3d4e5f6"
    assert not fit_exists(fit_prefix)

    table = pd.DataFrame({
        "component": ["co60", "u238", "ndbd"],
        "mean": [12.5, 45.2, 0.5],
        "std": [3.1, 6.7, 1.2],
    })
    data_evs = np.array([1.0, 2.0, 3.0])
    meta = {"output_uuid": "a1b2c3d4e5f6", "model": "m2"}

    # Write without trace
    h5_p, nc_p = write_fit(fit_prefix, table, trace=None, data=data_evs, metadata=meta, save_trace=False)
    assert h5_p.exists()
    assert nc_p is None
    assert fit_exists(fit_prefix, check_trace=False)
    assert not fit_exists(fit_prefix, check_trace=True)

    r_table, r_trace, r_data, r_meta = read_fit(fit_prefix)
    pd.testing.assert_frame_equal(r_table, table)
    np.testing.assert_array_equal(r_data, data_evs)
    assert r_trace is None
    assert r_meta["output_uuid"] == "a1b2c3d4e5f6"


def test_study_management(tmp_path):
    study_dir = tmp_path / "studies" / "test_study"
    paths = get_study_paths(study_dir, "test_study")

    # 1. Manifest
    manifest_data = {"study_name": "test_study", "points": 2}
    write_study_manifest(paths["manifest"], manifest_data)
    loaded_manifest = read_study_manifest(paths["manifest"])
    assert loaded_manifest["study_name"] == "test_study"

    # 2. SQLite Index
    init_study_db(paths["sqlite"])
    record_study_point(paths["sqlite"], {
        "point_id": "pt_001",
        "study_name": "test_study",
        "task_id": 1,
        "output_uuid": "pt_001",
        "status": "pending",
        "summary_json": {"prior_mean": 10},
    })

    # 3. Task Results CSV collection
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)
    df1 = pd.DataFrame([{"task_id": 1, "output_uuid": "pt_001", "mean_nev": 10.2}])
    df2 = pd.DataFrame([{"task_id": 2, "output_uuid": "pt_002", "mean_nev": 20.5}])
    df1.to_csv(tasks_dir / "task_1.csv", index=False)
    df2.to_csv(tasks_dir / "task_2.csv", index=False)

    merged = collect_study_results(study_dir, "test_study")
    assert len(merged) == 2
    assert set(merged["output_uuid"]) == {"pt_001", "pt_002"}
    assert paths["results_csv"].exists()
