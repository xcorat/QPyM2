"""Pure unit tests for Co-60 energy shift (co_shift) marginalization in QPyM2 v0.6.

Does not require external NERSC datasets, live ROOT data, or PyROOT.
"""

import json
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from omegaconf import OmegaConf

from qpym2.studies.marginalize import run_marginalization
from qpym2.studies.point import StudyPoint
from qpym2.studies.study import Study
from qpym2.studies.syst_co_shift import (
    calculate_co_shift_weights,
    collect_syst_co_shift,
    fit_co_shift_parabola,
    plot_co_shift_parabola,
)
from qpym2.syst.io import load_syst_result


class MockPosteriorVar:
    def __init__(self, values):
        self.values = np.asarray(values)


class MockLogLikelihood:
    def __init__(self, pik_values):
        self.pik = MockPosteriorVar(pik_values)

    def __contains__(self, item):
        return item == "pik"

    def __getitem__(self, item):
        if item == "pik":
            return self.pik
        raise KeyError(item)


class MockTrace:
    def __init__(self, ndbd_samples, logp_samples=None):
        self.posterior = {
            "ndbd": MockPosteriorVar(ndbd_samples),
            "co60": MockPosteriorVar(ndbd_samples * 1.5),
        }
        default_logp = np.random.normal(loc=-100.0, scale=1.0, size=len(ndbd_samples))
        self.log_likelihood = MockLogLikelihood(
            logp_samples if logp_samples is not None else default_logp
        )


def test_fit_co_shift_parabola():
    """Verify quadratic parabola fitting and peak vertex finding."""
    shifts = [0.0, 0.8, 1.0, 1.2]
    # S(x) = -10 * (x - 0.9)^2 - 50 = -10*(x^2 - 1.8x + 0.81) - 50 = -10x^2 + 18x - 58.1
    scores = [-10.0 * (s - 0.9) ** 2 - 50.0 for s in shifts]

    res = fit_co_shift_parabola(shifts, scores)
    assert np.isclose(res["vertex_shift"], 0.9, atol=1e-3)
    assert np.isclose(res["vertex_score"], -50.0, atol=1e-3)
    assert res["sigma_shift"] > 0
    # sigma = 1 / sqrt(-2 * -10) = 1 / sqrt(20) ~ 0.2236
    assert np.isclose(res["sigma_shift"], 1.0 / np.sqrt(20.0), atol=1e-3)


def test_calculate_co_shift_weights():
    """Verify softmax weighting and draw-count normalization."""
    scores = np.array([-100.0, -98.0, -97.5, -99.0])
    weights = calculate_co_shift_weights(scores)

    assert len(weights) == 4
    assert np.isclose(np.sum(weights), 1.0)
    assert weights[2] > weights[1] > weights[3] > weights[0]

    # Test with varying draw counts
    draws = np.array([10000, 20000, 20000, 10000])
    weights_scaled = calculate_co_shift_weights(scores, draw_counts=draws)
    assert np.isclose(np.sum(weights_scaled), 1.0)


def test_plot_co_shift_parabola(tmp_path):
    """Verify diagnostic parabola figure generation."""
    shifts = [0.0, 0.8, 1.0, 1.2]
    scores = [-10.0 * (s - 0.9) ** 2 - 50.0 for s in shifts]
    parabola = fit_co_shift_parabola(shifts, scores)

    fig, ax = plot_co_shift_parabola(shifts, scores, parabola)
    out_file = tmp_path / "test_parabola.png"
    fig.savefig(out_file)
    plt.close(fig)

    assert out_file.exists()
    assert out_file.stat().st_size > 0


def test_collect_syst_co_shift_end_to_end(tmp_path):
    """Verify end-to-end Co-shift marginalization with synthetic traces and bias calibration."""
    study_dir = tmp_path / "test_co_shift_study"
    study_dir.mkdir(parents=True)
    plots_dir = study_dir / "plots"
    plots_dir.mkdir(parents=True)

    # Setup synthetic bias calibration
    calib_data = {
        "p0_co_shift_0.0": {"a": 1.0, "b": 0.0},
        "p1_co_shift_0.8": {"a": 1.02, "b": -0.5},
        "p2_co_shift_1.0": {"a": 1.01, "b": -0.2},
        "p3_co_shift_1.2": {"a": 0.99, "b": 0.3},
    }
    calib_file = study_dir / "bias_calibration.json"
    with open(calib_file, "w") as f:
        json.dump(calib_data, f)

    # Create synthetic StudyPoints
    points = [
        StudyPoint(point_id="p0_co_shift_0.0", task_id=0, co_shift=0.0),
        StudyPoint(point_id="p1_co_shift_0.8", task_id=1, co_shift=0.8),
        StudyPoint(point_id="p2_co_shift_1.0", task_id=2, co_shift=1.0),
        StudyPoint(point_id="p3_co_shift_1.2", task_id=3, co_shift=1.2),
    ]

    cfg = OmegaConf.create(
        {
            "study": {"name": "test_co_shift_study", "type": "syst"},
            "reference": {"co_shift": 1.0, "co_cut": "gamma10_v150", "ndbd_cutoff": "u1740"},
            "syst": {
                "category": "co_shift",
                "reference_co_shift": 1.0,
                "bias_calibration": str(calib_file),
                "estimator": "mode",
                "lp_var": "pik",
            },
            "infrastructure": {
                "data_root": str(tmp_path),
                "output_paths": {"fit_data_dir": str(tmp_path / "fit_data")},
            },
        }
    )

    study = Study(name="test_co_shift_study", study_dir=study_dir, points=points, cfg=cfg)

    # Generate synthetic traces
    rng = np.random.default_rng(42)
    mock_traces = {
        0.0: MockTrace(
            rng.normal(loc=14.0, scale=4.0, size=10000),
            logp_samples=rng.normal(loc=-110.0, scale=1.0, size=1000),
        ),
        0.8: MockTrace(
            rng.normal(loc=15.0, scale=4.0, size=10000),
            logp_samples=rng.normal(loc=-102.0, scale=1.0, size=1000),
        ),
        1.0: MockTrace(
            rng.normal(loc=15.5, scale=4.0, size=10000),
            logp_samples=rng.normal(loc=-100.0, scale=1.0, size=1000),
        ),
        1.2: MockTrace(
            rng.normal(loc=15.2, scale=4.0, size=10000),
            logp_samples=rng.normal(loc=-105.0, scale=1.0, size=1000),
        ),
    }

    def mock_from_netcdf(path_str):
        for s, tr in mock_traces.items():
            if f"cs_{s}" in path_str or f"shift_{s}" in path_str or f"{s}" in path_str:
                return tr
        return mock_traces[1.0]

    # Create dummy NetCDF paths so check passes
    dummy_nc_files = []
    for p in points:
        nc_f = study_dir / f"fit_dummy_{p.point_id}.nc"
        nc_f.write_text("fake nc")
        dummy_nc_files.append(nc_f)

    def mock_get_fit_hash(self, point_cfg):
        return {
            "output_uuid": f"dummy_{self.point_id}",
            "fit_prefix": study_dir / f"fit_dummy_{self.point_id}",
        }

    with (
        patch("arviz.from_netcdf", side_effect=mock_from_netcdf),
        patch.object(StudyPoint, "get_fit_hash", mock_get_fit_hash),
    ):
        res = collect_syst_co_shift(study, verbose=True)

    summary = res["summary"]
    assert "nominal" in summary
    assert "co_shift_marginal" in summary

    nom_summary = summary["nominal"]
    assert nom_summary["category"] == "co_shift"
    assert nom_summary["reference_co_shift"] == 1.0
    assert "delta_right" in nom_summary
    assert "delta_left" in nom_summary
    assert "delta_mode" in nom_summary
    assert nom_summary["delta_right"] >= 0.0
    assert nom_summary["delta_left"] <= 0.0

    # Verify JSON summary written to disk
    summary_file = study_dir / "syst_summary.json"
    assert summary_file.exists()
    with open(summary_file) as f:
        loaded_json = json.load(f)
    assert "nominal" in loaded_json

    # Verify Format v2 HDF5 written to disk
    h5_file = study_dir / "syst_co_shift_marginal.h5"
    assert h5_file.exists()
    h5_data = load_syst_result(h5_file)
    assert h5_data["name"] == "co_shift"
    assert len(h5_data["bc_post"]) > 0

    # Verify diagnostic plots
    assert (plots_dir / "co_shift_logp_parabola.png").exists()
    assert (plots_dir / "co_shift_posterior_overlay.png").exists()
    assert (plots_dir / "co_shift_marginal_comparison.png").exists()


def test_co_shift_integration_with_marginalize(tmp_path):
    """Verify that run_marginalization seamlessly ingests the co_shift summary."""
    co_shift_dir = tmp_path / "example_co_shift"
    co_shift_dir.mkdir(parents=True)
    rate_study_dir = tmp_path / "example_syst_rate"
    rate_study_dir.mkdir(parents=True)

    # Write syst_summary.json in co_shift study
    co_shift_summary = {
        "nominal": {
            "category": "co_shift",
            "reference_co_shift": 1.0,
            "delta_right": 0.45,
            "delta_left": -0.55,
            "delta_mode": -0.10,
        }
    }
    with open(co_shift_dir / "syst_summary.json", "w") as f:
        json.dump(co_shift_summary, f)

    # Configure master rate study referencing the co_shift study
    rate_cfg = OmegaConf.create(
        {
            "study": {"name": "example_syst_rate", "type": "syst_rate"},
            "systematics": {
                "co_shift": {
                    "source": str(co_shift_dir),
                    "delta_right": 9.99,  # Should be overridden by summary file
                    "delta_left": -9.99,
                },
                "bias_corr": {"delta_right": 0.30, "delta_left": -0.04},
            },
            "marginalisation": {
                "exposure": 1038.4,
                "eff_mc": 0.04,
                "psa_mean": 0.99,
                "psa_sig": 0.01,
                "n_smear": 5000,
                "rng_seed": 42,
            },
        }
    )

    rate_study = Study(
        name="example_syst_rate",
        study_dir=rate_study_dir,
        points=[StudyPoint(point_id="nominal", task_id=0, co_shift=1.0)],
        cfg=rate_cfg,
    )

    rng = np.random.default_rng(42)
    nom_post = np.clip(rng.normal(loc=15.0, scale=3.0, size=5000), 0.0, None)

    results = run_marginalization(rate_study, nominal_post=nom_post, verbose=False)

    assert "sigma_right" in results
    assert "sigma_left" in results
    assert results["half_life_limit_90_yr"] > 0

    # Verify delta table contains the values from co_shift syst_summary.json
    delta_df = pd.read_csv(rate_study_dir / "delta_table.csv")
    co_row = delta_df[delta_df["source"] == "co_shift"]
    assert not co_row.empty
    assert np.isclose(co_row["delta_right"].iloc[0], 0.45)
    assert np.isclose(co_row["delta_left"].iloc[0], -0.55)
