import pytest
from qpym2.cfg.hydra_adapter_new import get_config, resolve_data_path
from qpym2.cfg.config_builders_new import (
    build_staging_config,
    build_input_config,
    get_h2uv_model,
    get_comp_cfg,
    get_artifact_store_dir,
    get_staging_dir,
    get_fit_dir,
)


def test_default_config_composition():
    cfg = get_config()
    assert cfg.models.tech_note_200.roi.esum_min == 2400
    assert cfg.models.tech_note_200.roi.esum_max == 2570
    assert cfg.models.tech_note_200.roi.threshold == 200


def test_staging_config_decoupled_from_roi():
    cfg = get_config()
    staging_cfg = build_staging_config(cfg)
    assert "TotalEnergy >= 2390 && TotalEnergy <= 2570" in staging_cfg['cuts']
    assert staging_cfg['treename'] == "uvtree"
    assert staging_cfg['outpath'] == "/data/cuore/Analysis/MC/staging"


def test_build_input_config():
    cfg = get_config()
    input_cfg = build_input_config(cfg, model_key="tech_note_200")
    assert input_cfg['outpath'] == "/data/cuore/Analysis/MC/staging"
    assert input_cfg['treename'] == "uvtree"
    assert input_cfg['roi']['threshold'] == 200
    assert input_cfg['binning']['bin_size'] == 2.5


def test_get_h2uv_model_defaults():
    cfg = get_config()
    model = get_h2uv_model(cfg)
    assert model.name == "hmu_tech_note"
    assert model.esum_min == 2400
    assert model.esum_max == 2570
    assert model.threshold == 200
    assert model.bin_size == 2.5
    assert model.binsize == pytest.approx(2.5 / (2 ** 0.5))
    assert model.umin == pytest.approx(2400 / (2 ** 0.5))
    assert model.umax == pytest.approx(2570 / (2 ** 0.5))


def test_get_h2uv_model_named():
    cfg = get_config()
    model = get_h2uv_model(cfg, model_key="tech_note_150")
    assert model.name == "hmu_tech_note_150"
    assert model.threshold == 150
    assert model.bin_size == 1.0


def test_component_grouping_custom_model():
    cfg = get_config()
    # Default components: 4 groups (co60, u238, th232, others)
    default_comps = get_comp_cfg(cfg)
    assert len(default_comps) == 4
    comp_names = [c['name'] for c in default_comps]
    assert "u238" in comp_names
    assert "th232" in comp_names

    # Custom 3-component grouping on tech_note_3comp
    m3_comps = get_comp_cfg(cfg, model_key="tech_note_3comp")
    assert len(m3_comps) == 3
    m3_names = [c['name'] for c in m3_comps]
    assert "u_th" in m3_names
    assert "u238" not in m3_names


def test_path_resolution():
    cfg = get_config()
    assert resolve_data_path(cfg, "/absolute/path") == "/absolute/path"
    assert resolve_data_path(cfg, "relative/dir") == "/data/cuore/Analysis/MC/relative/dir"
    assert get_artifact_store_dir(cfg) == "/data/cuore/Analysis/MC/artifact_store"
    assert get_staging_dir(cfg) == "/data/cuore/Analysis/MC/staging"
    assert get_fit_dir(cfg) == "/data/cuore/Analysis/MC/tech_note/fits"


def test_nersc_environment_paths():
    cfg = get_config(environment="nersc")
    assert get_staging_dir(cfg) == "/global/cfs/projectdirs/cuore/scratch/xcorat/data/MC/staging"
