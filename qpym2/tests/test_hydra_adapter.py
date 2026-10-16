import pytest
from qpym2.cfg.hydra_adapter import HydraConfigManager


def test_list_options_returns_expected_keys():
    opts = HydraConfigManager.list_options()
    assert isinstance(opts, dict)
    for key in ("environments", "analysis", "hyperparams"):
        assert key in opts
        assert isinstance(opts[key], list)


def test_list_options_has_some_entries():
    opts = HydraConfigManager.list_options()
    # Ensure at least one of the sections is non-empty (repository may vary)
    assert any(len(v) > 0 for v in opts.values())


def test_get_default_selections():
    defaults = HydraConfigManager.get_default_selections()
    assert isinstance(defaults, dict)
    # These are the defaults declared in config/config.yaml for this repo
    assert defaults.get('environment') == 'base'
    assert defaults.get('analysis') == 'models_base'
    assert defaults.get('hyperparams') == 'base'
