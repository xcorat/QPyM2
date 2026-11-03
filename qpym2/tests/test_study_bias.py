"""
Pure unit tests for bias toy study execution and collection in the storage-study layer.

Does not require external NERSC datasets, live ROOT data, or PyROOT.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from omegaconf import OmegaConf

from qpym2.studies.point import StudyPoint
from qpym2.studies.study import Study


class MockPosteriorVar:
    def __init__(self, values):
        self.values = np.asarray(values)


class MockLogLikelihood:
    def __init__(self, pik_values):
        self.pik = MockPosteriorVar(pik_values)


class MockTrace:
    def __init__(self, ndbd_samples, logp_samples=None):
        self.posterior = {
            "ndbd": MockPosteriorVar(ndbd_samples),
            "co60": MockPosteriorVar(ndbd_samples * 2),
        }
        default_logp = np.random.normal(loc=-100.0, scale=1.0, size=100)
        self.log_likelihood = MockLogLikelihood(logp_samples if logp_samples is not None else default_logp)


def test_check_bias_toys_status(tmp_path):
    """Verify StudyPoint.check_bias_toys detects pending, in-progress, and completed states."""
    study_dir = tmp_path / "bias_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    point = StudyPoint(point_id="p0_shift_0.0", task_id=0, co_shift=0.0)

    # 1. Pending (file does not exist)
    status = point.check_bias_toys(study_dir, n_toys=10)
    assert status["is_completed"] is False
    assert status["status"] == "PENDING"
    assert status["n_toys_completed"] == 0

    # 2. In-progress (fewer rows than target)
    task_csv = tasks_dir / "toys_bias_p0_shift_0.0.csv"
    pd.DataFrame({"toy_id": [0, 1, 2], "ndbd_inj_": [10.0, 20.0, 30.0]}).to_csv(task_csv, index=False)
    status = point.check_bias_toys(study_dir, n_toys=10)
    assert status["is_completed"] is False
    assert status["status"] == "IN_PROGRESS"
    assert status["n_toys_completed"] == 3

    # 3. Completed (rows >= n_toys)
    pd.DataFrame({"toy_id": range(10), "ndbd_inj_": range(10)}).to_csv(task_csv, index=False)
    status = point.check_bias_toys(study_dir, n_toys=10)
    assert status["is_completed"] is True
    assert status["status"] == "COMPLETED"
    assert status["n_toys_completed"] == 10


def test_study_point_run_bias_toys_with_mock(tmp_path, monkeypatch):
    """Verify StudyPoint.run_bias_toys generates pseudo-data and writes task CSV."""
    study_dir = tmp_path / "bias_study"
    study_dir.mkdir()

    point = StudyPoint(point_id="nominal", task_id=0, co_shift=0.0)

    # Mock setup_storage_fit
    mock_comps = pd.DataFrame(
        {
            "integral": [100.0, 5.0],
            "pdf_hist": [np.ones((4, 4)) / 16.0, np.ones((4, 4)) / 16.0],
        },
        index=["co60", "ndbd"],
    )
    mock_edges = (np.linspace(2400, 2570, 5), np.linspace(0, 100, 5))
    mock_setup = {
        "fit_comps": mock_comps,
        "edges": mock_edges,
    }

    monkeypatch.setattr(
        "qpym2.storage.fits.setup_storage_fit",
        lambda **kwargs: mock_setup,
    )

    # Fast mock fit function: returns predictable mode and upp34 around injection
    def mock_fit(data, comps, xe, ye):
        # Return samples centered around roughly 25
        samples = np.random.normal(loc=25.0, scale=3.0, size=500)
        return MockTrace(samples)

    dummy_cfg = OmegaConf.create({})
    csv_path = point.run_bias_toys(
        cfg=dummy_cfg,
        study_dir=study_dir,
        n_toys=4,
        bias_min=10.0,
        bias_max=50.0,
        seed=123,
        fit_fn=mock_fit,
        verbose=False,
    )

    assert Path(csv_path).is_file()
    df = pd.read_csv(csv_path)
    assert len(df) == 4
    assert "ndbd_inj_" in df.columns
    assert "ndbd_mode" in df.columns
    assert "ndbd_upp34" in df.columns
    assert "pull" in df.columns
    assert "logp_mode_" in df.columns
    assert "co60_mode" in df.columns
    assert df.loc[0, "type_"] == "bias"
    assert df.loc[0, "point_id"] == "nominal"


def test_study_collect_bias_results(tmp_path):
    """Verify Study.collect_bias_results merges tasks and extracts (a, b) calibration."""
    study_dir = tmp_path / "bias_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    p0 = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)
    p1 = StudyPoint(point_id="p1", task_id=1, co_shift=0.8)

    # Create synthetic bias data for p0: y = 1.05 * x - 3.0
    rng = np.random.default_rng(42)
    inj_0 = rng.uniform(10, 150, size=30)
    mode_0 = 1.05 * inj_0 - 3.0 + rng.normal(0, 0.5, size=30)
    upp34_0 = mode_0 + 2.5  # constant upper error = 2.5

    df0 = pd.DataFrame({
        "toy_id": range(30),
        "task_id": 0,
        "point_id": "p0",
        "co_shift": 0.0,
        "type_": "bias",
        "ndbd_inj_": inj_0,
        "ndbd_mode": mode_0,
        "ndbd_p16": mode_0 - 2.5,
        "ndbd_p90": mode_0 + 4.0,
        "ndbd_upp34": upp34_0,
        "logp_mode_": -50.0,
    })
    df0.to_csv(tasks_dir / "toys_bias_p0.csv", index=False)

    # Create synthetic bias data for p1: y = 0.98 * x - 1.5
    inj_1 = rng.uniform(10, 150, size=30)
    mode_1 = 0.98 * inj_1 - 1.5 + rng.normal(0, 0.5, size=30)
    upp34_1 = mode_1 + 2.0

    df1 = pd.DataFrame({
        "toy_id": range(30),
        "task_id": 1,
        "point_id": "p1",
        "co_shift": 0.8,
        "type_": "bias",
        "ndbd_inj_": inj_1,
        "ndbd_mode": mode_1,
        "ndbd_p16": mode_1 - 2.0,
        "ndbd_p90": mode_1 + 3.5,
        "ndbd_upp34": upp34_1,
        "logp_mode_": -52.0,
    })
    df1.to_csv(tasks_dir / "toys_bias_p1.csv", index=False)

    study = Study(
        name="test_bias_study",
        cfg=OmegaConf.create({"bias": {"xmin": 12.0}}),
        points=[p0, p1],
        study_dir=study_dir,
    )

    out = study.collect_bias_results(xmin=12.0, verbose=False)
    results_df = out["results_df"]
    calibrations = out["calibration"]

    # Check results.csv
    assert len(results_df) == 60
    assert (study_dir / "results.csv").is_file()

    # Check calibration.json
    assert (study_dir / "bias_calibration.json").is_file()
    assert "p0" in calibrations
    assert "p1" in calibrations

    # Verify fitted slope roughly matches true slopes (1.05 and 0.98)
    assert abs(calibrations["p0"]["a"] - 1.05) < 0.05
    assert abs(calibrations["p0"]["b"] - (-3.0)) < 2.0
    assert abs(calibrations["p1"]["a"] - 0.98) < 0.05
    assert abs(calibrations["p1"]["b"] - (-1.5)) < 2.0

    # Verify SQLite index has point entries
    assert (study_dir / "test_bias_study.sqlite").is_file()
