"""
Pure unit tests for QPyM2 storage fit routines (qpym2/storage/fits.py).

Does not require external NERSC datasets or PyROOT.
"""

import numpy as np
import pandas as pd
import pytest
from omegaconf import OmegaConf

from qpym2.storage.backend import fit_exists, read_fit, write_fit
from qpym2.storage.fits import (
    assemble_fit_components,
    build_fit_hash,
    resolve_fit_template_uuids,
)


@pytest.fixture
def mock_cfg(tmp_path):
    """Minimal test configuration fixture."""
    data_root = tmp_path / "data"
    data_root.mkdir()

    cfg_dict = {
        "staging": {"version": "v1"},
        "infrastructure": {
            "data_root": str(data_root),
            "output_paths": {
                "fit_data_dir": "fit_data",
                "staging_dir": "staging",
            },
            "data_files": {
                "bkg_model_jags": "jags_bkg.h5",
                "bkg_syst_corrmat": "corrmat.csv",
            },
            "naming_conventions": {
                "treename": "uvtree",
                "data_fname": "data.root",
            },
        },
        "cuts": {
            "co_cuts": {"gamma10_v150": "v>150"},
            "ndbd_cutoff": {"u1740": "u>1740"},
        },
        "analysis": {
            "components": [
                {"name": "co60", "filter_regex": "60Co"},
                {"name": "u238", "filter_regex": "238U|210Pb"},
                {"name": "th232", "filter_regex": "232Th|208Pb"},
                {"name": "others", "filter_regex_inverse": "60Co|238U|210Pb|232Th|208Pb"},
            ],
        },
        "models": {
            "tech_note_200": {
                "roi": {
                    "esum_min": 2400.0,
                    "esum_max": 2570.0,
                    "threshold": 200.0,
                },
                "binning": {
                    "bin_size": 2.5,
                },
                "coord_type": "h2uv",
                "smooth": [1.2, 1.2],
            }
        },
    }
    return OmegaConf.create(cfg_dict)


def test_parameter_set_uuid_referencing(mock_cfg):
    """Test parameter-set UUID resolution and linkage to components."""
    # 1. Unshifted baseline
    tpl_uuids_0, base_0, co_0, ndbd_0 = resolve_fit_template_uuids(
        mock_cfg,
        model_key="tech_note_200",
        co_cut="gamma10_v150",
        co_shift=0.0,
    )
    # When co_shift == 0.0, co UUID equals base UUID
    assert co_0 == base_0 == ndbd_0
    assert tpl_uuids_0["co60"] == base_0
    assert tpl_uuids_0["u238"] == base_0
    assert tpl_uuids_0["ndbd"] == base_0

    # 2. Shifted Co-60
    tpl_uuids_8, base_8, co_8, ndbd_8 = resolve_fit_template_uuids(
        mock_cfg,
        model_key="tech_note_200",
        co_cut="gamma10_v150",
        co_shift=0.8,
    )
    # Invariant base and ndbd UUIDs must match the unshifted run
    assert base_8 == base_0
    assert ndbd_8 == ndbd_0
    # Co-60 UUID must differ due to the shift
    assert co_8 != base_8
    assert tpl_uuids_8["co60"] == co_8
    assert tpl_uuids_8["u238"] == base_8
    assert tpl_uuids_8["ndbd"] == ndbd_8


def test_build_fit_hash_determinism():
    """Test output UUID determinism and dependency invalidation."""
    templates = {
        "co60": "uuid_co_a1",
        "u238": "uuid_base_b2",
        "th232": "uuid_base_b2",
        "others": "uuid_base_b2",
        "ndbd": "uuid_base_b2",
    }
    mcmc = {"chains": 8, "draws": 20000, "tune": 1000}
    priors = {"prior_width": "5sig", "sigma_scale": 5.0}
    smoothing = {"enabled": True, "method": "smooth_nph2"}

    # Canonical hash
    p1, uuid1 = build_fit_hash("m2", "data", templates, smoothing, priors, None, mcmc)
    p2, uuid2 = build_fit_hash("m2", "data", templates, smoothing, priors, None, mcmc)
    assert uuid1 == uuid2
    assert len(uuid1) == 12

    # Decoupled MCMC: Changing draws/chains preserves output UUID (same mathematical target)
    mcmc_mod = {"chains": 4, "draws": 40000, "tune": 2000}
    _, uuid_mod_draws = build_fit_hash("m2", "data", templates, smoothing, priors, None, mcmc_mod)
    assert uuid_mod_draws == uuid1

    # Changing Co-60 input UUID changes output UUID
    templates_shifted = dict(templates, co60="uuid_co_different")
    _, uuid_shifted = build_fit_hash("m2", "data", templates_shifted, smoothing, priors, None, mcmc)
    assert uuid_shifted != uuid1

    # Changing priors changes output UUID
    priors_mod = {"prior_width": "syst", "sigma_scale": 1.0}
    _, uuid_mod_priors = build_fit_hash("m2", "data", templates, smoothing, priors_mod, None, mcmc)
    assert uuid_mod_priors != uuid1


def test_fit_template_assembly_mock(tmp_path):
    """Test assembling mock source templates into likelihood components."""
    # Create mock JAGS table
    jags_table = pd.DataFrame(
        {
            "mode": [100.0, 50.0, 30.0],
            "mean": [100.0, 50.0, 30.0],
            "std": [10.0, 5.0, 3.0],
            "fname": ["60Co.root", "238U.root", "232Th.root"],
        },
        index=["60Co", "238U", "232Th"],
    )

    def mock_read_jagsh5(fpath):
        return jags_table, pd.DataFrame()

    # Create mock template reader returning synthetic 2D Gaussian
    edges = [np.linspace(0, 10, 11), np.linspace(0, 10, 11)]
    hist_sample = np.ones((10, 10), dtype=float)

    def mock_read_template(path):
        return hist_sample.copy(), edges, {"meta": "mock"}

    comp_cfg = [
        {"name": "co60", "filter_regex": "60Co"},
        {"name": "u238", "filter_regex": "238U"},
        {"name": "th232", "filter_regex": "232Th"},
    ]

    class MockHM:
        smooth = None

    input_dir = tmp_path / "mock_input"
    input_dir.mkdir()
    # Create dummy files for existence checks
    for s in ["60Co", "238U", "232Th", "ndbd"]:
        (input_dir / f"{s}_uuid.h5").touch()

    mctable, sig_df, edges_out, comps_df = assemble_fit_components(
        input_dir=input_dir,
        jags_fpath="dummy.h5",
        comp_cfg=comp_cfg,
        hm=MockHM(),
        uuid_base="uuid",
        uuid_co="uuid",
        uuid_ndbd="uuid",
        read_jagsh5_fn=mock_read_jagsh5,
        read_template_fn=mock_read_template,
    )

    assert "mchist" in mctable.columns
    assert len(mctable) == 3
    assert len(sig_df) == 1
    assert "ndbd" in sig_df.index
    assert "pdf_hist" in comps_df.columns
    assert "integral" in comps_df.columns
    assert set(comps_df.index) == {"co60", "u238", "th232", "ndbd"}


def test_fit_write_and_read_roundtrip(tmp_path):
    """Test writing and reading fit outputs with and without trace."""
    out_dir = tmp_path / "output"
    prefix = out_dir / "fit_test123456"

    comps_df = pd.DataFrame(
        {"integral": [120.0, 45.0], "prior": [None, None]},
        index=["co60", "ndbd"],
    )
    data_events = np.array([2450.5, 2520.1, 2530.0])
    metadata = {"model_name": "m2", "output_uuid": "test123456"}

    # Write without trace
    h5_p, nc_p = write_fit(
        prefix,
        table=comps_df,
        trace=None,
        data=data_events,
        metadata=metadata,
        save_trace=False,
    )
    assert h5_p.exists()
    assert nc_p is None
    assert fit_exists(prefix)
    assert not fit_exists(prefix, check_trace=True)

    # Read back
    table_read, trace_read, data_read, meta_read = read_fit(prefix)
    assert set(table_read.index) == {"co60", "ndbd"}
    assert np.allclose(table_read["integral"].values, [120.0, 45.0])
    assert len(data_read) == 3
    assert meta_read["output_uuid"] == "test123456"


def test_setup_storage_fit_suffix_and_cache_quality(mock_cfg, tmp_path):
    """Test suffix support and quality-aware caching in setup_storage_fit."""
    from qpym2.storage.fits import setup_storage_fit

    # 1. Suffix creates distinct prefix but identical output_uuid
    res_nosuffix = setup_storage_fit(mock_cfg, model_key="tech_note_200", dry_run=True)
    res_suffix = setup_storage_fit(mock_cfg, model_key="tech_note_200", suffix="d40k", dry_run=True)

    assert res_nosuffix["output_uuid"] == res_suffix["output_uuid"]
    assert str(res_suffix["fit_prefix"]).endswith(f"fit_{res_suffix['output_uuid']}_d40k")
    assert str(res_nosuffix["fit_prefix"]).endswith(f"fit_{res_nosuffix['output_uuid']}")
    assert res_suffix["fit_params"]["suffix"] == "d40k"

    # 2. Quality-aware caching: write existing fit with 1000 draws
    out_dir = res_nosuffix["output_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = res_nosuffix["fit_prefix"]

    comps_df = pd.DataFrame({"integral": [10.0]}, index=["ndbd"])
    meta = {
        "output_uuid": res_nosuffix["output_uuid"],
        "mcmc": {"draws": 1000, "chains": 4},
    }
    write_fit(prefix, table=comps_df, metadata=meta)

    # Requesting 2000 draws should detect lower cached quality -> is_cached=False
    check_high = setup_storage_fit(
        mock_cfg, model_key="tech_note_200", draws=2000, chains=4, dry_run=True
    )
    assert check_high["is_cached"] is False

    # Requesting 500 draws should detect sufficient quality -> is_cached=True
    check_low = setup_storage_fit(
        mock_cfg, model_key="tech_note_200", draws=500, chains=4, dry_run=True
    )
    assert check_low["is_cached"] is True
