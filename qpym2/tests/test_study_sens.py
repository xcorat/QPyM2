"""
Pure unit tests for sensitivity toy study execution and collection in the storage-study layer.

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


def test_check_sens_toys_status(tmp_path):
    """Verify StudyPoint.check_sens_toys detects pending, in-progress, and completed states."""
    study_dir = tmp_path / "sens_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    point = StudyPoint(point_id="p0_shift_0.0", task_id=0, co_shift=0.0)

    # 1. Pending (file does not exist)
    status = point.check_sens_toys(study_dir, n_toys=10)
    assert status["is_completed"] is False
    assert status["status"] == "PENDING"
    assert status["n_toys_completed"] == 0

    # 2. In-progress (fewer rows than target)
    task_csv = tasks_dir / "toys_sens_p0_shift_0.0.csv"
    pd.DataFrame({"toy_id": [0, 1, 2], "ndbd_inj_": [0.0, 0.0, 0.0]}).to_csv(task_csv, index=False)
    status = point.check_sens_toys(study_dir, n_toys=10)
    assert status["is_completed"] is False
    assert status["status"] == "IN_PROGRESS"
    assert status["n_toys_completed"] == 3

    # 3. Completed (rows >= n_toys)
    pd.DataFrame({"toy_id": range(10), "ndbd_inj_": [0.0] * 10}).to_csv(task_csv, index=False)
    status = point.check_sens_toys(study_dir, n_toys=10)
    assert status["is_completed"] is True
    assert status["status"] == "COMPLETED"
    assert status["n_toys_completed"] == 10


def test_study_point_run_sens_toys_with_mock(tmp_path, monkeypatch):
    """Verify StudyPoint.run_sens_toys generates pseudo-data under H0 and writes task CSV."""
    study_dir = tmp_path / "sens_study"
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

    # Fast mock fit function: returns predictable upper limit and mode
    def mock_fit(data, comps, xe, ye):
        # Half-normal-like distribution around 0 for background-only pseudo-data
        samples = np.abs(np.random.normal(loc=0.0, scale=2.5, size=500))
        return MockTrace(samples)

    dummy_cfg = OmegaConf.create({})
    csv_path = point.run_sens_toys(
        cfg=dummy_cfg,
        study_dir=study_dir,
        n_toys=4,
        seed=123,
        fit_fn=mock_fit,
        verbose=False,
    )

    assert Path(csv_path).is_file()
    df = pd.read_csv(csv_path)
    assert len(df) == 4
    assert "ndbd_inj_" in df.columns
    assert (df["ndbd_inj_"] == 0.0).all()
    assert "ndbd_mode" in df.columns
    assert "ndbd_p90" in df.columns
    assert "ndbd_upp34" in df.columns
    assert "pull" in df.columns
    assert "logp_mode_" in df.columns
    assert "co60_mode" in df.columns
    assert df.loc[0, "type_"] == "sens"
    assert df.loc[0, "point_id"] == "nominal"


def test_study_collect_sens_results(tmp_path):
    """Verify Study.collect_sens_results merges tasks, computes quantiles, applies bias correction, and writes JSON/DB."""
    study_dir = tmp_path / "sens_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    p0 = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)
    p1 = StudyPoint(point_id="p1", task_id=1, co_shift=0.8)

    # Synthetic sensitivity toy results for p0: ndbd_p90 centered around 10.0
    rng = np.random.default_rng(42)
    p90_0 = 10.0 + rng.normal(0, 1.0, size=25)
    mode_0 = np.abs(rng.normal(0, 0.5, size=25))
    df0 = pd.DataFrame({
        "toy_id": range(25),
        "task_id": 0,
        "point_id": "p0",
        "co_shift": 0.0,
        "type_": "sens",
        "ndbd_inj_": 0.0,
        "ndbd_mode": mode_0,
        "ndbd_p16": 0.0,
        "ndbd_p90": p90_0,
        "ndbd_upp34": 3.0,
        "logp_mode_": -50.0,
    })
    df0.to_csv(tasks_dir / "toys_sens_p0.csv", index=False)

    # Synthetic sensitivity toy results for p1: ndbd_p90 centered around 12.0
    p90_1 = 12.0 + rng.normal(0, 1.0, size=25)
    mode_1 = np.abs(rng.normal(0, 0.5, size=25))
    df1 = pd.DataFrame({
        "toy_id": range(25),
        "task_id": 1,
        "point_id": "p1",
        "co_shift": 0.8,
        "type_": "sens",
        "ndbd_inj_": 0.0,
        "ndbd_mode": mode_1,
        "ndbd_p16": 0.0,
        "ndbd_p90": p90_1,
        "ndbd_upp34": 3.5,
        "logp_mode_": -52.0,
    })
    df1.to_csv(tasks_dir / "toys_sens_p1.csv", index=False)

    # Mock bias calibration: y = 1.0 * x + 1.0 -> bc = (y - 1.0) / 1.0
    bias_calib = {
        "p0": {"a": 1.0, "b": 1.0, "co_shift": 0.0},
        "p1": {"a": 1.0, "b": 2.0, "co_shift": 0.8},
    }

    study = Study(
        name="test_sens_study",
        cfg=OmegaConf.create({"sens": {"bias_calibration": None}}),
        points=[p0, p1],
        study_dir=study_dir,
    )

    out = study.collect_sens_results(bias_calibration=bias_calib, verbose=False)
    results_df = out["results_df"]
    summaries = out["summary"]

    assert len(results_df) == 50
    assert (study_dir / "results.csv").is_file()
    assert (study_dir / "sens_summary.json").is_file()
    assert (study_dir / "test_sens_study.sqlite").is_file()

    assert "p0" in summaries
    assert "p1" in summaries

    p0_summary = summaries["p0"]
    assert p0_summary["n_toys"] == 25
    assert abs(p0_summary["sens_p90_median"] - float(np.median(p90_0))) < 1e-4
    assert p0_summary["sens_p90_p16"] <= p0_summary["sens_p90_median"] <= p0_summary["sens_p90_p84"]
    assert p0_summary["sens_p90_p025"] <= p0_summary["sens_p90_p975"]

    # Verify bias correction was applied: bc_median = (median - 1.0)
    assert abs(p0_summary["sens_bc_p90_median"] - (p0_summary["sens_p90_median"] - 1.0)) < 1e-4
    p1_summary = summaries["p1"]
    assert abs(p1_summary["sens_bc_p90_median"] - (p1_summary["sens_p90_median"] - 2.0)) < 1e-4
