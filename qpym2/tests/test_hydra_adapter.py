from qpym2.cfg.hydra_adapter import HydraConfigManager, get_config


def test_list_options_returns_expected_keys():
    opts = HydraConfigManager.list_options()
    assert isinstance(opts, dict)
    for key in ("environments", "staging", "analysis", "models"):
        assert key in opts
        assert isinstance(opts[key], list)
        assert len(opts[key]) > 0


def test_get_default_selections():
    defaults = HydraConfigManager.get_default_selections()
    assert isinstance(defaults, dict)
    assert defaults.get('environment') == 'local'
    assert defaults.get('staging') == 'default'
    assert defaults.get('analysis') == 'base'
    assert defaults.get('models') == 'tech_note'


def test_get_config_defaults():
    cfg = get_config()
    assert cfg.pipeline.version == "0.6.0"
    assert "tech_note_200" in cfg.models
    assert cfg.staging.version == 1
    assert "TotalEnergy >= 2390 && TotalEnergy <= 2570" in cfg.staging.cuts.broad_preselection


def test_get_config_overrides():
    cfg = get_config(environment="nersc")
    assert cfg.infrastructure.data_root == "/global/cfs/projectdirs/cuore/scratch/xcorat/data"
