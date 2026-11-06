"""
Pure unit tests for MC-stats (bmmc) systematic study execution and collection.

Does not require external NERSC datasets, live ROOT data, or PyROOT.
"""

from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from omegaconf import OmegaConf

from qpym2.io.bkg_model import get_mcsyst_samples
from qpym2.studies.point import StudyPoint
from qpym2.studies.study import Study
from qpym2.syst.io import load_syst_result


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
            "co60": MockPosteriorVar(ndbd_samples * 1.5),
        }
        default_logp = np.random.normal(loc=-100.0, scale=1.0, size=100)
        self.log_likelihood = MockLogLikelihood(logp_samples if logp_samples is not None else default_logp)


def test_get_mcsyst_samples_direct(tmp_path):
    """Verify get_mcsyst_samples in stat_only and covariance perturbation modes."""
    # Create synthetic chain DataFrame
    chain_df = pd.DataFrame({
        "co60": [10.0, 11.0, 12.0, 13.0, 14.0],
        "u238": [50.0, 52.0, 49.0, 51.0, 53.0],
    })

    # 1. Stat-only mode (pure draws from chain)
    samples_stat = get_mcsyst_samples(chain=chain_df, corrmat_fpath=None, n=3, seed=42, stat_only=True)
    assert len(samples_stat) == 3
    assert list(samples_stat.columns) == ["co60", "u238"]
    # All values must originate from the chain rows
    for _, row in samples_stat.iterrows():
        matches = ((chain_df["co60"] == row["co60"]) & (chain_df["u238"] == row["u238"])).any()
        assert matches

    # 2. Seed determinism
    samples_stat_2 = get_mcsyst_samples(chain=chain_df, corrmat_fpath=None, n=3, seed=42, stat_only=True)
    pd.testing.assert_frame_equal(samples_stat, samples_stat_2)

    # 3. Covariance perturbation mode
    corrmat_csv = tmp_path / "corrmat.csv"
    cov_df = pd.DataFrame(
        [[1.0, 0.2], [0.2, 2.0]],
        index=["co60", "u238"],
        columns=["co60", "u238"],
    )
    cov_df.to_csv(corrmat_csv)

    samples_cov = get_mcsyst_samples(chain=chain_df, corrmat_fpath=corrmat_csv, n=4, seed=123, stat_only=False)
    assert len(samples_cov) == 4
    assert list(samples_cov.columns) == ["co60", "u238"]
    # Perturbed samples should not be exact chain rows
    assert not np.allclose(samples_cov.values, chain_df.iloc[:4].values)


def test_check_syst_mc_status(tmp_path):
    """Verify StudyPoint.check_syst_mc detects pending, in-progress, and completed states."""
    study_dir = tmp_path / "bmmc_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    point = StudyPoint(point_id="p0_shift_0.0", task_id=0, co_shift=0.0)

    # 1. Pending (neither CSV nor H5 exists)
    status = point.check_syst_mc(study_dir, n_samples=10)
    assert status["is_completed"] is False
    assert status["status"] == "PENDING"
    assert status["n_samples_completed"] == 0

    # 2. In-progress (CSV exists with fewer rows than target or H5 missing)
    task_csv = tasks_dir / "syst_mc_p0_shift_0.0.csv"
    pd.DataFrame({"sample_idx": [0, 1, 2], "ndbd_mode": [5.0, 5.2, 4.8]}).to_csv(task_csv, index=False)
    status = point.check_syst_mc(study_dir, n_samples=10)
    assert status["is_completed"] is False
    assert status["status"] == "IN_PROGRESS"
    assert status["n_samples_completed"] == 3

    # 3. Completed (rows >= n_samples and H5 exists)
    task_h5 = tasks_dir / "syst_mc_p0_shift_0.0_samples.h5"
    pd.DataFrame({"sample_idx": range(10), "ndbd_mode": [5.0] * 10}).to_csv(task_csv, index=False)
    with h5py.File(task_h5, "w") as h5f:
        h5f.create_dataset("ndbd_samples", data=np.ones((10, 50)))
        h5f.create_dataset("max_logp", data=np.ones(10))

    status = point.check_syst_mc(study_dir, n_samples=10)
    assert status["is_completed"] is True
    assert status["status"] == "COMPLETED"
    assert status["n_samples_completed"] == 10


def test_study_point_run_syst_mc_with_mock(tmp_path, monkeypatch):
    """Verify StudyPoint.run_syst_mc perturbs templates, executes fits, and writes CSV/H5 tasks."""
    study_dir = tmp_path / "bmmc_study"
    study_dir.mkdir()

    point = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)

    # Mock storage resolution functions
    mock_comps = pd.DataFrame(
        {
            "integral": [100.0, 5.0],
            "pdf_hist": [np.ones((4, 4)) / 16.0, np.ones((4, 4)) / 16.0],
            "prior": [None, None],
        },
        index=["co60", "ndbd"],
    )
    mock_edges = (np.linspace(2400, 2570, 5), np.linspace(0, 100, 5))
    mock_setup = {
        "data": np.ones((4, 4)) * 5.0,
        "fit_comps": mock_comps,
        "edges": mock_edges,
    }

    monkeypatch.setattr("qpym2.storage.fits.setup_storage_fit", lambda **kwargs: mock_setup)

    class MockHM:
        smooth = None

    monkeypatch.setattr(
        "qpym2.storage.fits.resolve_storage_coords",
        lambda *args, **kwargs: (None, MockHM(), None, None),
    )
    monkeypatch.setattr(
        "qpym2.storage.fits.resolve_fit_template_uuids",
        lambda *args, **kwargs: ({}, "base_uuid", "co_uuid", "ndbd_uuid"),
    )

    mock_mctable = pd.DataFrame(
        {
            "mchist": [np.ones((4, 4)) / 16.0],
            "norm": [100.0],
            "mean": [100.0],
            "std": [10.0],
            "mode": [100.0],
        },
        index=["co60_source"],
    )
    mock_signal = pd.DataFrame(
        {
            "pdf_hist": [np.ones((4, 4)) / 16.0],
            "integral": [5.0],
            "prior": [None],
        },
        index=["ndbd"],
    )
    monkeypatch.setattr(
        "qpym2.storage.fits.assemble_fit_components",
        lambda **kwargs: (mock_mctable, mock_signal, None, None),
    )
    monkeypatch.setattr(
        "qpym2.cfg.config_builders.get_comp_cfg",
        lambda *args, **kwargs: [{"name": "co60", "filter_regex": "co60"}],
    )
    monkeypatch.setattr(
        "qpym2.cfg.config_builders.get_jags_path",
        lambda *args, **kwargs: "dummy_jags.h5",
    )
    monkeypatch.setattr(
        "qpym2.cfg.config_builders.get_corrmat_path",
        lambda *args, **kwargs: None,
    )

    # Mock get_mcsyst_samples
    def mock_get_mcsyst(chain, corrmat_fpath, n, seed, stat_only):
        return pd.DataFrame({"co60_source": [95.0 + i for i in range(n)]})

    monkeypatch.setattr("qpym2.io.bkg_model.get_mcsyst_samples", mock_get_mcsyst)

    # Fast mock fit function
    def mock_fit(data, comps, xe, ye):
        samples = np.abs(np.random.normal(loc=12.0, scale=2.0, size=200))
        return MockTrace(samples, logp_samples=np.random.normal(loc=-45.0, scale=0.5, size=50))

    dummy_cfg = OmegaConf.create({})
    csv_path, h5_path = point.run_syst_mc(
        cfg=dummy_cfg,
        study_dir=study_dir,
        n_samples=3,
        seed=42,
        stat_only=True,
        fit_fn=mock_fit,
        verbose=False,
    )

    assert Path(csv_path).is_file()
    assert Path(h5_path).is_file()

    df = pd.read_csv(csv_path)
    assert len(df) == 3
    assert "sample_idx" in df.columns
    assert "point_id" in df.columns
    assert "ndbd_mode" in df.columns
    assert "max_logp" in df.columns
    assert "ndbd_p90" in df.columns
    assert df.loc[0, "point_id"] == "p0"

    with h5py.File(h5_path, "r") as h5f:
        ndbd_arr = h5f["ndbd_samples"][:]
        max_logp_arr = h5f["max_logp"][:]
        assert ndbd_arr.shape == (3, 200)
        assert len(max_logp_arr) == 3
        assert h5f.attrs["point_id"] == "p0"
        assert h5f.attrs["n_samples"] == 3


def test_study_collect_syst_results(tmp_path):
    """Verify Study.collect_syst_results combines mixture posteriors, computes deltas, and updates DB/HDF5."""
    study_dir = tmp_path / "bmmc_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    p0 = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)
    p1 = StudyPoint(point_id="p1", task_id=1, co_shift=0.8)

    rng = np.random.default_rng(42)

    # Mock task CSV and HDF5 for p0 (samples centered around 10.0)
    n_samples = 4
    draws_per_fit = 300
    p0_samples = rng.normal(loc=10.0, scale=2.0, size=(n_samples, draws_per_fit))
    p0_logp = np.array([-40.0, -41.0, -39.5, -42.0])

    df0 = pd.DataFrame({
        "sample_idx": range(n_samples),
        "task_id": 0,
        "point_id": "p0",
        "co_shift": 0.0,
        "ndbd_shift": 0.0,
        "co_cut": "None",
        "ndbd_cutoff": "None",
        "stat_only": False,
        "ndbd_mode": [10.0] * n_samples,
        "ndbd_p16": [8.0] * n_samples,
        "ndbd_p90": [12.0] * n_samples,
        "max_logp": p0_logp,
        "ndbd_mode_bc": [np.nan] * n_samples,
        "ndbd_p16_bc": [np.nan] * n_samples,
        "ndbd_p90_bc": [np.nan] * n_samples,
    })
    df0.to_csv(tasks_dir / "syst_mc_p0.csv", index=False)
    with h5py.File(tasks_dir / "syst_mc_p0_samples.h5", "w") as h5f:
        h5f.create_dataset("ndbd_samples", data=p0_samples)
        h5f.create_dataset("max_logp", data=p0_logp)

    # Mock task CSV and HDF5 for p1 (samples centered around 11.5)
    p1_samples = rng.normal(loc=11.5, scale=2.2, size=(n_samples, draws_per_fit))
    p1_logp = np.array([-42.0, -43.0, -41.5, -44.0])

    df1 = pd.DataFrame({
        "sample_idx": range(n_samples),
        "task_id": 1,
        "point_id": "p1",
        "co_shift": 0.8,
        "ndbd_shift": 0.0,
        "co_cut": "None",
        "ndbd_cutoff": "None",
        "stat_only": False,
        "ndbd_mode": [11.5] * n_samples,
        "ndbd_p16": [9.3] * n_samples,
        "ndbd_p90": [13.7] * n_samples,
        "max_logp": p1_logp,
        "ndbd_mode_bc": [np.nan] * n_samples,
        "ndbd_p16_bc": [np.nan] * n_samples,
        "ndbd_p90_bc": [np.nan] * n_samples,
    })
    df1.to_csv(tasks_dir / "syst_mc_p1.csv", index=False)
    with h5py.File(tasks_dir / "syst_mc_p1_samples.h5", "w") as h5f:
        h5f.create_dataset("ndbd_samples", data=p1_samples)
        h5f.create_dataset("max_logp", data=p1_logp)

    # Bias calibration dictionary
    bias_calib = {
        "p0": {"a": 1.0, "b": 1.0, "co_shift": 0.0},
        "p1": {"a": 1.0, "b": 1.5, "co_shift": 0.8},
    }

    # Nominal posterior: e.g. Gaussian centered at 9.8
    nom_post = rng.normal(loc=9.8, scale=1.9, size=2000)

    study = Study(
        name="test_bmmc_study",
        cfg=OmegaConf.create({"syst": {"category": "bmmc"}}),
        points=[p0, p1],
        study_dir=study_dir,
    )

    out = study.collect_syst_results(
        bias_calibration=bias_calib,
        nominal_post=nom_post,
        verbose=False,
    )

    results_df = out["results_df"]
    summaries = out["summary"]

    assert len(results_df) == 8
    assert (study_dir / "results.csv").is_file()
    assert (study_dir / "syst_summary.json").is_file()
    assert (study_dir / "test_bmmc_study.sqlite").is_file()

    assert "p0" in summaries
    assert "p1" in summaries

    p0_summary = summaries["p0"]
    assert p0_summary["category"] == "bmmc"
    assert p0_summary["n_samples"] == 4
    assert "delta_mode" in p0_summary
    assert "delta_right" in p0_summary
    assert "delta_left" in p0_summary
    assert p0_summary["delta_left"] <= 0.0
    assert p0_summary["delta_right"] >= 0.0

    # Verify per-systematic Format v2 HDF5 output
    syst_h5_p0 = study_dir / "syst_bmmc_p0.h5"
    assert syst_h5_p0.is_file()
    res_p0 = load_syst_result(syst_h5_p0)
    assert res_p0["name"] == "bmmc"
    assert len(res_p0["bc_post"]) > 0
    assert len(res_p0["raw_post"]) > 0

    # Verify diagnostic plot was generated
    assert (study_dir / "plots" / "bmmc_comparison_p0.png").is_file()
    assert (study_dir / "plots" / "bmmc_comparison_p1.png").is_file()


def test_study_check_and_print_syst_status(tmp_path):
    """Verify Study.check_syst and Study.print_syst_status format output cleanly."""
    study_dir = tmp_path / "bmmc_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    p0 = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)
    p1 = StudyPoint(point_id="p1", task_id=1, co_shift=0.8)

    study = Study(
        name="test_status_study",
        cfg=OmegaConf.create({"syst": {"category": "bmmc"}}),
        points=[p0, p1],
        study_dir=study_dir,
    )

    df_status = study.check_syst(n_samples=5)
    assert len(df_status) == 2
    assert "point_id" in df_status.columns
    assert "n_samples_completed" in df_status.columns
    assert "status" in df_status.columns
    assert (df_status["status"] == "PENDING").all()

    # Verify print_syst_status runs without error
    study.print_syst_status(df_status)
