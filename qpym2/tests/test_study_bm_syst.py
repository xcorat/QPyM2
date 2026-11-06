"""
Pure unit tests for Background Model (bm_threshold, bm_binning, bin_size)
systematics studies and master rate marginalisation in QPyM2 v0.6.

Does not require external NERSC datasets, live ROOT data, or PyROOT.
"""

from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from omegaconf import OmegaConf

from qpym2.io.bkg_model import read_bm_variant_modes, list_bm_variants
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


def test_read_bm_variant_modes_and_list_synthetic(tmp_path):
    """Verify read_bm_variant_modes and list_bm_variants on synthetic ResultsRaw.h5."""
    results_raw_path = tmp_path / "ResultsRaw.h5"

    df_thr250 = pd.DataFrame({
        "MC": ["co60_source", "u238_source", "th232_source"],
        "IsGaus": [True, True, False],
        "PDF_Mean": [105.0, 52.0, 12.0],
    })
    df_bin_5_20 = pd.DataFrame({
        "MC": ["co60_source", "u238_source", "th232_source"],
        "IsGaus": [True, True, True],
        "PDF_Mean": [98.0, 49.0, 15.0],
    })

    with pd.HDFStore(str(results_raw_path), mode="w") as store:
        store.put("Posterior_Raw_Threshold250", df_thr250)
        store.put("Posterior_Raw_MinBinSize_5_20", df_bin_5_20)
        store.put("Unrelated_Table", pd.DataFrame({"a": [1]}))

    # Test list_bm_variants
    variants = list_bm_variants(results_raw_path)
    assert variants == ["MinBinSize_5_20", "Threshold250"]

    # Test read_bm_variant_modes for Threshold250
    modes_thr250 = read_bm_variant_modes(results_raw_path, "Threshold250")
    assert modes_thr250["co60_source"] == 105.0
    assert modes_thr250["u238_source"] == 52.0
    # False IsGaus should be 0.0 (upper limit)
    assert modes_thr250["th232_source"] == 0.0

    # Test read_bm_variant_modes with Posterior_Raw_ prefix
    modes_prefix = read_bm_variant_modes(results_raw_path, "Posterior_Raw_Threshold250")
    pd.testing.assert_series_equal(modes_thr250, modes_prefix)

    # Test read_bm_variant_modes for MinBinSize_5_20
    modes_bin = read_bm_variant_modes(results_raw_path, "MinBinSize_5_20")
    assert modes_bin["co60_source"] == 98.0
    assert modes_bin["u238_source"] == 49.0
    assert modes_bin["th232_source"] == 15.0


def test_check_syst_bm_threshold_status(tmp_path):
    """Verify StudyPoint.check_syst_bm_threshold tracks completion properly."""
    study_dir = tmp_path / "thr_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    point = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)

    # 1. Pending (neither CSV nor H5 exists)
    status = point.check_syst_bm_threshold(study_dir, variant="Threshold250")
    assert status["is_completed"] is False
    assert status["status"] == "PENDING"
    assert status["n_samples_completed"] == 0

    # 2. In-progress (CSV exists with row but H5 is missing)
    task_csv = tasks_dir / "syst_bm_threshold_p0.csv"
    pd.DataFrame([{"variant": "Threshold250", "ndbd_mode": 10.0}]).to_csv(task_csv, index=False)
    status = point.check_syst_bm_threshold(study_dir, variant="Threshold250")
    assert status["is_completed"] is False
    assert status["status"] == "IN_PROGRESS"

    # 3. Completed (both CSV and non-empty H5 exist)
    task_h5 = tasks_dir / "syst_bm_threshold_p0_samples.h5"
    with h5py.File(task_h5, "w") as h5f:
        h5f.create_dataset("ndbd_samples", data=np.ones((1, 50)))
        h5f.create_dataset("max_logp", data=np.array([-40.0]))

    status = point.check_syst_bm_threshold(study_dir, variant="Threshold250")
    assert status["is_completed"] is True
    assert status["status"] == "COMPLETED"
    assert status["n_samples_completed"] == 1


def test_check_syst_bm_binning_status(tmp_path):
    """Verify StudyPoint.check_syst_bm_binning tracks multiple variants."""
    study_dir = tmp_path / "binning_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    point = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)
    variants = ["MinBinSize_5_20", "MinBinSize_25_60"]

    # 1. Pending
    status = point.check_syst_bm_binning(study_dir, variants=variants)
    assert status["is_completed"] is False
    assert status["status"] == "PENDING"
    assert status["n_samples_completed"] == 0

    # 2. In-progress (1 of 2 variants completed in CSV)
    task_csv = tasks_dir / "syst_bm_binning_p0.csv"
    pd.DataFrame([{"variant": "MinBinSize_5_20", "ndbd_mode": 10.0}]).to_csv(task_csv, index=False)
    status = point.check_syst_bm_binning(study_dir, variants=variants)
    assert status["is_completed"] is False
    assert status["status"] == "IN_PROGRESS"
    assert status["n_samples_completed"] == 1

    # 3. Completed (2 rows and H5 exists)
    task_h5 = tasks_dir / "syst_bm_binning_p0_samples.h5"
    pd.DataFrame([
        {"variant": "MinBinSize_5_20", "ndbd_mode": 10.0},
        {"variant": "MinBinSize_25_60", "ndbd_mode": 10.5},
    ]).to_csv(task_csv, index=False)
    with h5py.File(task_h5, "w") as h5f:
        h5f.create_dataset("ndbd_samples", data=np.ones((2, 50)))
        h5f.create_dataset("max_logp", data=np.array([-40.0, -41.0]))

    status = point.check_syst_bm_binning(study_dir, variants=variants)
    assert status["is_completed"] is True
    assert status["status"] == "COMPLETED"
    assert status["n_samples_completed"] == 2


def test_study_point_run_syst_bm_threshold_mock(tmp_path, monkeypatch):
    """Verify StudyPoint.run_syst_bm_threshold overrides norms and executes fit."""
    study_dir = tmp_path / "thr_study"
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

    # Mock read_bm_variant_modes
    monkeypatch.setattr(
        "qpym2.io.bkg_model.read_bm_variant_modes",
        lambda *args, **kwargs: pd.Series({"co60_source": 115.0}),
    )

    # Mock fit function
    def mock_fit(data, comps, xe, ye):
        samples = np.abs(np.random.normal(loc=8.5, scale=1.5, size=200))
        return MockTrace(samples, logp_samples=np.random.normal(loc=-42.0, scale=0.5, size=50))

    dummy_cfg = OmegaConf.create({})
    csv_path, h5_path = point.run_syst_bm_threshold(
        cfg=dummy_cfg,
        study_dir=study_dir,
        variant="Threshold250",
        model_key_syst="tech_note_250",
        results_raw_path=tmp_path / "ResultsRaw.h5",
        fit_fn=mock_fit,
        verbose=False,
    )

    assert Path(csv_path).is_file()
    assert Path(h5_path).is_file()

    df = pd.read_csv(csv_path)
    assert len(df) == 1
    assert df.loc[0, "variant"] == "Threshold250"
    assert df.loc[0, "model_key_syst"] == "tech_note_250"
    assert "ndbd_mode" in df.columns

    with h5py.File(h5_path, "r") as h5f:
        assert h5f["ndbd_samples"].shape == (1, 200)
        assert h5f.attrs["variant"] == "Threshold250"
        assert h5f.attrs["category"] == "bm_threshold"


def test_study_point_run_syst_bm_binning_mock(tmp_path, monkeypatch):
    """Verify StudyPoint.run_syst_bm_binning fits each binning variant."""
    study_dir = tmp_path / "binning_study"
    study_dir.mkdir()

    point = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)

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

    def mock_read_modes(result_file, fit_name):
        shift = 5.0 if "5_20" in fit_name else -5.0
        return pd.Series({"co60_source": 100.0 + shift})

    monkeypatch.setattr("qpym2.io.bkg_model.read_bm_variant_modes", mock_read_modes)

    def mock_fit(data, comps, xe, ye):
        samples = np.abs(np.random.normal(loc=11.0, scale=1.8, size=150))
        return MockTrace(samples, logp_samples=np.random.normal(loc=-41.0, scale=0.5, size=50))

    dummy_cfg = OmegaConf.create({})
    csv_path, h5_path = point.run_syst_bm_binning(
        cfg=dummy_cfg,
        study_dir=study_dir,
        variants=["MinBinSize_5_20", "MinBinSize_25_60"],
        results_raw_path=tmp_path / "ResultsRaw.h5",
        fit_fn=mock_fit,
        verbose=False,
    )

    assert Path(csv_path).is_file()
    assert Path(h5_path).is_file()

    df = pd.read_csv(csv_path)
    assert len(df) == 2
    assert list(df["variant"]) == ["MinBinSize_5_20", "MinBinSize_25_60"]

    with h5py.File(h5_path, "r") as h5f:
        assert h5f["ndbd_samples"].shape == (2, 150)
        assert h5f.attrs["category"] == "bm_binning"


def test_study_collect_syst_results_bm_threshold(tmp_path):
    """Verify Study.collect_syst_results computes deltas and saves Format v2 for bm_threshold."""
    study_dir = tmp_path / "thr_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    point = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)
    rng = np.random.default_rng(42)

    # Mock threshold fit output (samples centered at 8.84, nominal centered at 10.0 => delta ~ -1.16)
    thr_samples = rng.normal(loc=8.84, scale=1.8, size=(1, 20000))
    task_csv = tasks_dir / "syst_bm_threshold_p0.csv"
    task_h5 = tasks_dir / "syst_bm_threshold_p0_samples.h5"

    pd.DataFrame([{
        "variant": "Threshold250",
        "task_id": 0,
        "point_id": "p0",
        "model_key": "tech_note_200",
        "model_key_syst": "tech_note_250",
        "co_shift": 0.0,
        "ndbd_shift": 0.0,
        "ndbd_mode": 8.84,
        "ndbd_p16": 7.0,
        "ndbd_p90": 10.6,
        "max_logp": -42.0,
        "ndbd_mode_bc": np.nan,
        "ndbd_p16_bc": np.nan,
        "ndbd_p90_bc": np.nan,
    }]).to_csv(task_csv, index=False)

    with h5py.File(task_h5, "w") as h5f:
        h5f.create_dataset("ndbd_samples", data=thr_samples)
        h5f.create_dataset("max_logp", data=np.array([-42.0]))
        h5f.attrs["variant"] = "Threshold250"
        h5f.attrs["category"] = "bm_threshold"

    nom_post = rng.normal(loc=10.0, scale=1.8, size=20000)

    study = Study(
        name="test_bm_threshold_study",
        cfg=OmegaConf.create({"syst": {"category": "bm_threshold", "variant": "Threshold250"}}),
        points=[point],
        study_dir=study_dir,
    )

    out = study.collect_syst_results(nominal_post=nom_post, verbose=False)
    results_df = out["results_df"]
    summaries = out["summary"]

    assert len(results_df) == 1
    assert "p0" in summaries
    p0_summary = summaries["p0"]
    assert p0_summary["category"] == "bm_threshold"
    assert p0_summary["delta_left"] <= 0.0
    assert p0_summary["delta_right"] >= 0.0
    # Mode shifted from 10.0 to ~8.84 => delta_left should be ~ -1.16, delta_right = 0.0
    assert np.isclose(p0_summary["delta_left"], -1.16, atol=0.35)
    assert p0_summary["delta_right"] == 0.0

    # Format v2 HDF5 check
    syst_h5 = study_dir / "syst_bm_threshold_p0.h5"
    assert syst_h5.is_file()
    loaded = load_syst_result(syst_h5)
    assert loaded["name"] == "bm_threshold"
    assert len(loaded["raw_post"]) == 20000


def test_study_collect_syst_results_bm_binning(tmp_path):
    """Verify Study.collect_syst_results computes envelope across BM binning variants."""
    study_dir = tmp_path / "binning_study"
    tasks_dir = study_dir / "tasks"
    tasks_dir.mkdir(parents=True)

    point = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)
    rng = np.random.default_rng(42)

    # Variant 1 shifts mode down (-0.6), Variant 2 shifts mode up (+0.4)
    v1_samples = rng.normal(loc=9.4, scale=1.8, size=20000)
    v2_samples = rng.normal(loc=10.4, scale=1.8, size=20000)
    samples_mat = np.vstack([v1_samples, v2_samples])

    task_csv = tasks_dir / "syst_bm_binning_p0.csv"
    task_h5 = tasks_dir / "syst_bm_binning_p0_samples.h5"

    pd.DataFrame([
        {
            "variant": "MinBinSize_5_20",
            "task_id": 0,
            "point_id": "p0",
            "model_key": "tech_note_200",
            "co_shift": 0.0,
            "ndbd_shift": 0.0,
            "ndbd_mode": 9.4,
            "ndbd_p16": 7.6,
            "ndbd_p90": 11.2,
            "max_logp": -40.0,
            "ndbd_mode_bc": np.nan,
            "ndbd_p16_bc": np.nan,
            "ndbd_p90_bc": np.nan,
        },
        {
            "variant": "MinBinSize_25_60",
            "task_id": 0,
            "point_id": "p0",
            "model_key": "tech_note_200",
            "co_shift": 0.0,
            "ndbd_shift": 0.0,
            "ndbd_mode": 10.4,
            "ndbd_p16": 8.6,
            "ndbd_p90": 12.2,
            "max_logp": -41.0,
            "ndbd_mode_bc": np.nan,
            "ndbd_p16_bc": np.nan,
            "ndbd_p90_bc": np.nan,
        },
    ]).to_csv(task_csv, index=False)

    with h5py.File(task_h5, "w") as h5f:
        h5f.create_dataset("ndbd_samples", data=samples_mat)
        h5f.create_dataset("max_logp", data=np.array([-40.0, -41.0]))
        h5f.attrs["category"] = "bm_binning"

    nom_post = rng.normal(loc=10.0, scale=1.8, size=20000)

    study = Study(
        name="test_bm_binning_study",
        cfg=OmegaConf.create({"syst": {"category": "bm_binning"}}),
        points=[point],
        study_dir=study_dir,
    )

    out = study.collect_syst_results(nominal_post=nom_post, verbose=False)
    summaries = out["summary"]

    p0_summary = summaries["p0"]
    assert p0_summary["category"] == "bm_binning"
    # Envelope: delta_left should be ~ -0.60, delta_right should be ~ +0.40
    assert np.isclose(p0_summary["delta_left"], -0.60, atol=0.25)
    assert np.isclose(p0_summary["delta_right"], 0.40, atol=0.25)

    # Format v2 HDF5 check
    syst_h5 = study_dir / "syst_bm_binning_p0.h5"
    assert syst_h5.is_file()
    loaded = load_syst_result(syst_h5)
    assert loaded["name"] == "bm_binning"


def test_study_run_marginalization(tmp_path):
    """Verify Study.run_marginalization combines deltas, smears posteriors, and calculates limits."""
    study_dir = tmp_path / "rate_study"
    study_dir.mkdir()

    cfg = OmegaConf.create({
        "study": {"name": "test_rate_study", "type": "syst_rate"},
        "systematics": {
            "lss_dq": {"delta_right": 4.74, "delta_left": -3.59},
            "lss_sscl": {"delta_right": 1.20, "delta_left": -1.20},
            "bmmc": {"delta_right": 0.50, "delta_left": -1.50},
            "bm_threshold": {"delta_right": 0.00, "delta_left": -1.16},
            "bm_binning": {"delta_right": 0.00, "delta_left": -0.60},
            "ndbd_shift": {"delta_right": 0.00, "delta_left": -0.42},
            "bias_corr": {"delta_right": 0.30, "delta_left": -0.04},
        },
        "marginalisation": {
            "exposure": 1000.0,
            "eff_mc": 0.04,
            "psa_mean": 0.99,
            "psa_sig": 0.01,
            "n_smear": 5000,
            "rng_seed": 42,
            "extend_method": "skew_norm",
        },
    })

    point = StudyPoint(point_id="p0", task_id=0, co_shift=0.0)
    study = Study(
        name="test_rate_study",
        cfg=cfg,
        points=[point],
        study_dir=study_dir,
    )

    rng = np.random.default_rng(42)
    # Mock nominal posterior centered at 15.0 counts
    nom_post = np.clip(rng.normal(loc=15.0, scale=4.0, size=5000), 0.0, None)

    res = study.run_marginalization(nominal_post=nom_post, verbose=False)

    assert "sigma_left" in res
    assert "sigma_right" in res
    # Expected quadrature: sigma_left ~ sqrt(3.59^2 + 1.2^2 + 1.5^2 + 1.16^2 + 0.6^2 + 0.42^2 + 0.04^2) ~ 4.35
    # sigma_right ~ sqrt(4.74^2 + 1.2^2 + 0.5^2 + 0.3^2) ~ 4.92
    assert 4.0 <= abs(res["sigma_left"]) <= 4.7
    assert 4.6 <= res["sigma_right"] <= 5.3

    assert "half_life_limit_90_yr" in res
    assert res["half_life_limit_90_yr"] > 0.0

    # Verify generated artifacts
    assert (study_dir / "delta_table.csv").is_file()
    assert (study_dir / "syst_rate_summary.json").is_file()
    assert (study_dir / "combined_result.h5").is_file()
    assert (study_dir / "plots" / "delta_breakdown.png").is_file()
    assert (study_dir / "plots" / "smeared_posterior.png").is_file()
    assert (study_dir / "plots" / "rate_posterior.png").is_file()
