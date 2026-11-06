"""
Pure unit tests for Lineshape Systematic (LSS: sscl and dq) study execution and collection.

Does not require external NERSC datasets, live ROOT data, or PyROOT.
"""

from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import pytest
from omegaconf import OmegaConf

from qpym2.studies.point import StudyPoint
from qpym2.studies.study import Study
from qpym2.syst.io import load_syst_result


class MockPosteriorVar:
    def __init__(self, values):
        self.values = np.asarray(values)


class MockPosterior:
    def __init__(self, vars_dict):
        self._vars = vars_dict
        self.data_vars = list(vars_dict.keys())

    def __getitem__(self, key):
        return self._vars[key]


class MockTrace:
    def __init__(self, ndbd_samples):
        self.posterior = MockPosterior(
            {
                "ndbd": MockPosteriorVar(ndbd_samples),
                "co60": MockPosteriorVar(ndbd_samples * 1.5),
            }
        )


def test_check_syst_lss_status(tmp_path):
    """Verify StudyPoint.check_syst_lss detects pending, in-progress, and completed states for sscl and dq."""
    study_dir = tmp_path / "lss_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    point = StudyPoint(point_id="p0_cs1p0", task_id=0, co_shift=1.0)
    datasets = list(range(3601, 3616))  # 15 datasets

    # 1. sscl: Pending
    st_sscl = point.check_syst_lss(study_dir, category="lss_sscl", datasets=datasets)
    assert st_sscl["is_completed"] is False
    assert st_sscl["status"] == "PENDING"
    assert st_sscl["n_samples_target"] == 15
    assert st_sscl["n_samples_completed"] == 0

    # 2. sscl: In-progress (fewer than 15 rows)
    task_csv = tasks_dir / "syst_lss_sscl_p0_cs1p0.csv"
    pd.DataFrame({"ds": [3601, 3602], "ndbd_mode": [12.0, 12.5]}).to_csv(task_csv, index=False)
    st_sscl = point.check_syst_lss(study_dir, category="lss_sscl", datasets=datasets)
    assert st_sscl["is_completed"] is False
    assert st_sscl["status"] == "IN_PROGRESS"
    assert st_sscl["n_samples_completed"] == 2

    # 3. sscl: Completed (15 rows)
    pd.DataFrame({"ds": datasets, "ndbd_mode": [12.0] * 15}).to_csv(task_csv, index=False)
    st_sscl = point.check_syst_lss(study_dir, category="lss_sscl", datasets=datasets)
    assert st_sscl["is_completed"] is True
    assert st_sscl["status"] == "COMPLETED"
    assert st_sscl["n_samples_completed"] == 15

    # 4. dq (pm1): Pending
    st_dq = point.check_syst_lss(study_dir, category="lss_dq", datasets=datasets, dq_mode="pm1")
    assert st_dq["is_completed"] is False
    assert st_dq["status"] == "PENDING"
    assert st_dq["n_samples_target"] == 30  # 15 * 2

    # 5. dq (pm1): Completed (30 rows)
    task_dq_csv = tasks_dir / "syst_lss_dq_p0_cs1p0.csv"
    dq_rows = []
    for d in datasets:
        dq_rows.append({"ds": d, "variant": "dqres_1r", "ndbd_mode": 12.1})
        dq_rows.append({"ds": d, "variant": "dqres_1l", "ndbd_mode": 11.9})
    pd.DataFrame(dq_rows).to_csv(task_dq_csv, index=False)
    st_dq = point.check_syst_lss(study_dir, category="lss_dq", datasets=datasets, dq_mode="pm1")
    assert st_dq["is_completed"] is True
    assert st_dq["status"] == "COMPLETED"
    assert st_dq["n_samples_completed"] == 30


def test_study_point_run_lss_sscl_with_mock(tmp_path, monkeypatch):
    """Verify StudyPoint.run_syst_lss executes sscl fits across datasets and saves CSV/H5 tasks."""
    study_dir = tmp_path / "lss_sscl_study"
    study_dir.mkdir()

    point = StudyPoint(point_id="p0", task_id=0, co_shift=1.0)
    datasets = [3601, 3602]

    # Mock storage setup
    mock_comps = pd.DataFrame(
        {"integral": [100.0, 5.0], "prior": [None, None]},
        index=["co60", "ndbd"],
    )
    mock_setup = {
        "data": np.ones((4, 4)) * 5.0,
        "fit_comps": mock_comps,
        "edges": (np.linspace(2400, 2570, 5), np.linspace(0, 100, 5)),
    }
    monkeypatch.setattr("qpym2.storage.fits.setup_storage_fit", lambda **kwargs: mock_setup)

    class MockHM:
        umin = 1600.0
        umax = 1900.0
        vmin = 0.0
        vmax = 200.0
        binsize = 5.0
        threshold = 100.0
        smooth = None
        skip_roi = None

    monkeypatch.setattr(
        "qpym2.storage.fits.resolve_storage_coords",
        lambda *args, **kwargs: (None, MockHM(), None, None),
    )
    monkeypatch.setattr(
        "qpym2.storage.fits.resolve_fit_template_uuids",
        lambda *args, **kwargs: ({}, "base_uuid", "co_uuid", "ndbd_uuid"),
    )
    monkeypatch.setattr(
        "qpym2.storage.fits.assemble_fit_components",
        lambda **kwargs: (pd.DataFrame(index=["co60"]), pd.DataFrame(index=["ndbd"]), None, None),
    )
    monkeypatch.setattr(
        "qpym2.cfg.config_builders.get_comp_cfg",
        lambda *args, **kwargs: [{"name": "co60", "filter_regex": "co60"}],
    )
    monkeypatch.setattr(
        "qpym2.cfg.config_builders.get_jags_path",
        lambda *args, **kwargs: "dummy_jags.h5",
    )

    # Mock lss_params
    mock_lss_params = {
        "sigma0_2d": np.ones((2, 100)),
        "delta_p1_1r": np.array([0.05, 0.06]),
    }

    # Mock evdata
    mock_evdata = {
        "co60": {
            "u0": np.array([1700.0, 1720.0]),
            "v0": np.array([50.0, 60.0]),
            "ch1": np.array([1, 2]),
            "ch2": np.array([3, 4]),
            "ds": np.array([3601, 3602]),
        },
        "ndbd": {
            "u0": np.array([1780.0, 1785.0]),
            "v0": np.array([10.0, 15.0]),
            "ch1": np.array([1, 2]),
            "ch2": np.array([3, 4]),
            "ds": np.array([3601, 3602]),
        },
    }

    # Mock fit function
    def mock_fit(mt, sig, data, hm, comp_cfg, priors):
        samples = np.abs(np.random.normal(loc=13.0, scale=1.5, size=200))
        return MockTrace(samples)

    dummy_cfg = OmegaConf.create({})
    csv_path, h5_path = point.run_syst_lss(
        cfg=dummy_cfg,
        study_dir=study_dir,
        category="lss_sscl",
        datasets=datasets,
        lss_params=mock_lss_params,
        evdata=mock_evdata,
        fit_fn=mock_fit,
        verbose=False,
    )

    assert Path(csv_path).is_file()
    assert Path(h5_path).is_file()

    df = pd.read_csv(csv_path)
    assert len(df) == 2
    assert list(df["ds"].values) == [3601, 3602]
    assert list(df["variant"].values) == ["sscale_1r", "sscale_1r"]
    assert "ndbd_mode" in df.columns
    assert "ndbd_p16" in df.columns
    assert "ndbd_p90" in df.columns

    with h5py.File(h5_path, "r") as h5f:
        ndbd_arr = h5f["ndbd_samples"][:]
        assert ndbd_arr.shape == (2, 200)
        assert h5f.attrs["point_id"] == "p0"
        assert h5f.attrs["category"] == "lss_sscl"


def test_study_point_run_lss_dq_with_mock(tmp_path, monkeypatch):
    """Verify StudyPoint.run_syst_lss executes dq fits for both 1r and 1l variants in pm1 mode."""
    study_dir = tmp_path / "lss_dq_study"
    study_dir.mkdir()

    point = StudyPoint(point_id="p0", task_id=0, co_shift=1.0)
    datasets = [3601, 3602]

    # Mock storage setup
    mock_comps = pd.DataFrame(
        {"integral": [100.0, 5.0], "prior": [None, None]},
        index=["co60", "ndbd"],
    )
    mock_setup = {
        "data": np.ones((4, 4)) * 5.0,
        "fit_comps": mock_comps,
        "edges": (np.linspace(2400, 2570, 5), np.linspace(0, 100, 5)),
    }
    monkeypatch.setattr("qpym2.storage.fits.setup_storage_fit", lambda **kwargs: mock_setup)

    class MockHM:
        umin = 1600.0
        umax = 1900.0
        vmin = 0.0
        vmax = 200.0
        binsize = 5.0
        threshold = 100.0
        smooth = None
        skip_roi = None

    monkeypatch.setattr(
        "qpym2.storage.fits.resolve_storage_coords",
        lambda *args, **kwargs: (None, MockHM(), None, None),
    )
    monkeypatch.setattr(
        "qpym2.storage.fits.resolve_fit_template_uuids",
        lambda *args, **kwargs: ({}, "base_uuid", "co_uuid", "ndbd_uuid"),
    )
    monkeypatch.setattr(
        "qpym2.storage.fits.assemble_fit_components",
        lambda **kwargs: (pd.DataFrame(index=["co60"]), pd.DataFrame(index=["ndbd"]), None, None),
    )
    monkeypatch.setattr(
        "qpym2.cfg.config_builders.get_comp_cfg",
        lambda *args, **kwargs: [{"name": "co60", "filter_regex": "co60"}],
    )
    monkeypatch.setattr(
        "qpym2.cfg.config_builders.get_jags_path",
        lambda *args, **kwargs: "dummy_jags.h5",
    )

    # Mock lss_params
    mock_lss_params = {
        "q_pars_arr": np.zeros((2, 3)),
        "q_pars_sample_1r": np.ones((2, 3)) * 0.1,
        "q_pars_sample_1l": np.ones((2, 3)) * -0.1,
    }

    # Mock evdata
    mock_evdata = {
        "co60": {
            "u0": np.array([1700.0, 1720.0]),
            "v0": np.array([50.0, 60.0]),
            "ch1": np.array([1, 2]),
            "ch2": np.array([3, 4]),
            "ds": np.array([3601, 3602]),
        },
        "ndbd": {
            "u0": np.array([1780.0, 1785.0]),
            "v0": np.array([10.0, 15.0]),
            "ch1": np.array([1, 2]),
            "ch2": np.array([3, 4]),
            "ds": np.array([3601, 3602]),
        },
    }

    def mock_fit(mt, sig, data, hm, comp_cfg, priors):
        samples = np.abs(np.random.normal(loc=14.0, scale=1.5, size=200))
        return MockTrace(samples)

    dummy_cfg = OmegaConf.create({})
    csv_path, h5_path = point.run_syst_lss(
        cfg=dummy_cfg,
        study_dir=study_dir,
        category="lss_dq",
        datasets=datasets,
        dq_mode="pm1",
        lss_params=mock_lss_params,
        evdata=mock_evdata,
        fit_fn=mock_fit,
        verbose=False,
    )

    assert Path(csv_path).is_file()
    assert Path(h5_path).is_file()

    df = pd.read_csv(csv_path)
    assert len(df) == 4  # 2 datasets * 2 variants
    assert set(df["variant"].unique()) == {"dqres_1r", "dqres_1l"}
    assert set(df["ds"].unique()) == {3601, 3602}


def test_study_collect_syst_results_lss_sscl(tmp_path):
    """Verify Study.collect_syst_results performs signed quadrature summation and writes Format v2 HDF5 for lss_sscl."""
    study_dir = tmp_path / "lss_sscl_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    p0 = StudyPoint(point_id="p0", task_id=0, co_shift=1.0)
    datasets = list(range(3601, 3616))  # 15 datasets

    # Synthetic task CSV with known shifts:
    # Nominal mode = 10.0
    # ds 3601 has mode = 13.0 -> delta = +3.0
    # ds 3602 has mode = 6.0  -> delta = -4.0
    # ds 3603..3615 have mode = 10.0 -> delta = 0.0
    # Expected: delta_right = sqrt(3^2) = 3.0, delta_left = -sqrt(4^2) = -4.0
    modes = [10.0] * 15
    modes[0] = 13.0
    modes[1] = 6.0

    df_task = pd.DataFrame(
        {
            "ds": datasets,
            "variant": ["sscale_1r"] * 15,
            "category": ["lss_sscl"] * 15,
            "point_id": ["p0"] * 15,
            "task_id": [0] * 15,
            "co_shift": [1.0] * 15,
            "ndbd_shift": [0.0] * 15,
            "co_cut": ["gamma10_v150"] * 15,
            "ndbd_cutoff": ["u1740"] * 15,
            "ndbd_mode": modes,
            "ndbd_p16": [m - 2.0 for m in modes],
            "ndbd_p90": [m + 3.0 for m in modes],
            "ndbd_mode_bc": [np.nan] * 15,
            "ndbd_p16_bc": [np.nan] * 15,
            "ndbd_p90_bc": [np.nan] * 15,
        }
    )
    df_task.to_csv(tasks_dir / "syst_lss_sscl_p0.csv", index=False)

    study_cfg = OmegaConf.create(
        {
            "study": {"name": "test_lss_sscl"},
            "syst": {"category": "lss_sscl"},
        }
    )
    study = Study(
        name="test_lss_sscl",
        cfg=study_cfg,
        points=[p0],
        study_dir=study_dir,
    )

    nom_post = np.random.normal(loc=10.0, scale=1.0, size=1000)
    from qpym2.utils.mathy import calc_intervals, signed_quad_sum

    nom_mode = calc_intervals([nom_post], method="hdi")[0][0]
    sp, sm = signed_quad_sum(np.array(modes) - nom_mode)
    exp_right = float(sp)
    exp_left = -float(sm)

    res = study.collect_syst_results(nominal_post=nom_post, verbose=False)

    assert "summary" in res
    assert "p0" in res["summary"]
    summary_p0 = res["summary"]["p0"]

    assert summary_p0["delta_right"] == pytest.approx(exp_right, abs=1e-5)
    assert summary_p0["delta_left"] == pytest.approx(exp_left, abs=1e-5)
    assert summary_p0["category"] == "lss_sscl"

    # Verify Format v2 HDF5 result
    h5_path = study_dir / "syst_lss_sscl_p0.h5"
    assert h5_path.is_file()

    loaded = load_syst_result(h5_path)
    assert loaded["name"] == "lss_sscl"
    assert loaded["per_post"] is not None
    assert len(loaded["per_post"]) == 15
    assert "delta_mode" in loaded["per_post"].columns


def test_study_collect_syst_results_lss_dq(tmp_path):
    """Verify Study.collect_syst_results performs envelope selection and quadrature summation for lss_dq."""
    study_dir = tmp_path / "lss_dq_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    p0 = StudyPoint(point_id="p0", task_id=0, co_shift=1.0)

    # Nominal mode = 10.0
    # ds 3601: 1r=12.0 (+2.0), 1l=9.0 (-1.0) -> dpos=2.0, dneg=-1.0
    # ds 3602: 1r=7.0 (-3.0), 1l=8.0 (-2.0) -> dpos=0.0, dneg=-3.0
    # Expected delta_right = sqrt(2^2 + 0) = 2.0
    # Expected delta_left = -sqrt((-1)^2 + (-3)^2) = -sqrt(10) ≈ -3.162
    dq_rows = [
        {"ds": 3601, "variant": "dqres_1r", "ndbd_mode": 12.0, "ndbd_p16": 10.0, "ndbd_p90": 15.0},
        {"ds": 3601, "variant": "dqres_1l", "ndbd_mode": 9.0, "ndbd_p16": 7.0, "ndbd_p90": 12.0},
        {"ds": 3602, "variant": "dqres_1r", "ndbd_mode": 7.0, "ndbd_p16": 5.0, "ndbd_p90": 10.0},
        {"ds": 3602, "variant": "dqres_1l", "ndbd_mode": 8.0, "ndbd_p16": 6.0, "ndbd_p90": 11.0},
    ]
    for r in dq_rows:
        r.update(
            {
                "category": "lss_dq",
                "point_id": "p0",
                "task_id": 0,
                "co_shift": 1.0,
                "ndbd_shift": 0.0,
                "co_cut": "gamma10_v150",
                "ndbd_cutoff": "u1740",
                "ndbd_mode_bc": np.nan,
                "ndbd_p16_bc": np.nan,
                "ndbd_p90_bc": np.nan,
            }
        )
    pd.DataFrame(dq_rows).to_csv(tasks_dir / "syst_lss_dq_p0.csv", index=False)

    study_cfg = OmegaConf.create(
        {
            "study": {"name": "test_lss_dq"},
            "syst": {"category": "lss_dq", "dq_mode": "pm1"},
        }
    )
    study = Study(
        name="test_lss_dq",
        cfg=study_cfg,
        points=[p0],
        study_dir=study_dir,
    )

    nom_post = np.random.normal(loc=10.0, scale=0.8, size=1000)
    from qpym2.utils.mathy import calc_intervals

    nom_mode = calc_intervals([nom_post], method="hdi")[0][0]
    d1r, d1l = 12.0 - nom_mode, 9.0 - nom_mode
    d2r, d2l = 7.0 - nom_mode, 8.0 - nom_mode
    exp_right = float(np.sqrt(max(d1r, d1l, 0.0) ** 2 + max(d2r, d2l, 0.0) ** 2))
    exp_left = -float(np.sqrt(min(d1r, d1l, 0.0) ** 2 + min(d2r, d2l, 0.0) ** 2))

    res = study.collect_syst_results(nominal_post=nom_post, verbose=False)

    summary_p0 = res["summary"]["p0"]
    assert summary_p0["delta_right"] == pytest.approx(exp_right, abs=1e-5)
    assert summary_p0["delta_left"] == pytest.approx(exp_left, abs=1e-5)
    assert summary_p0["category"] == "lss_dq"

    h5_path = study_dir / "syst_lss_dq_p0.h5"
    assert h5_path.is_file()
    loaded = load_syst_result(h5_path)
    assert loaded["name"] == "lss_dq"
    assert loaded["per_post"] is not None
    assert len(loaded["per_post"]) == 4


def test_lss_data_path_resolution():
    """Verify resolution order and validation for LSS data paths."""
    from qpym2.cfg.syst_config import get_lss_data_path, build_lss_file_paths

    # 1. Null/empty raises ValueError
    cfg_empty = OmegaConf.create({"infrastructure": {"data_files": {"lss_data_path": None}}})
    with pytest.raises(ValueError, match="LSS data path is not configured"):
        get_lss_data_path(cfg_empty)

    with pytest.raises(ValueError, match="lss_data_path must be a valid directory path"):
        build_lss_file_paths(None, [3601])

    with pytest.raises(ValueError, match="lss_data_path must be a valid directory path"):
        build_lss_file_paths("", [3601])

    # 2. Resolved from infrastructure
    cfg_infra = OmegaConf.create(
        {"infrastructure": {"data_files": {"lss_data_path": "/path/from/infra"}}}
    )
    assert get_lss_data_path(cfg_infra) == "/path/from/infra"

    # 3. Resolved from syst override
    cfg_syst = OmegaConf.create(
        {
            "infrastructure": {"data_files": {"lss_data_path": None}},
            "syst": {"lss_data_path": "/path/from/syst"},
        }
    )
    assert get_lss_data_path(cfg_syst) == "/path/from/syst"
