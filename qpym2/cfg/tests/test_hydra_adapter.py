"""Unit tests for `qpym2.cfg.hydra_adapter`."""

from __future__ import annotations

from pathlib import Path

import pytest
from omegaconf import DictConfig, OmegaConf

from qpym2.cfg import ConfigLoadError, HydraConfigManager, get_config


# ----------------------------------------------------------------------
# Basic config loading
# ----------------------------------------------------------------------
class TestBasicLoading:
    def test_loads_default_config(self, use_test_configs):
        cfg = HydraConfigManager.get_config()
        assert isinstance(cfg, DictConfig)
        assert cfg.pipeline.name == "test"

    def test_default_composition_picks_base_env(self, use_test_configs):
        cfg = HydraConfigManager.get_config()
        # base.yaml sets data_dir = "data"
        assert cfg.infrastructure.output_paths.data_dir == "data"

    def test_default_composition_picks_base_analysis(self, use_test_configs):
        cfg = HydraConfigManager.get_config()
        assert cfg.models.base.threshold == 1.0

    def test_get_config_module_shortcut(self, use_test_configs):
        cfg = get_config()
        assert cfg.pipeline.name == "test"


# ----------------------------------------------------------------------
# Override behaviour
# ----------------------------------------------------------------------
class TestOverrides:
    def test_environment_group_override(self, use_test_configs):
        cfg = HydraConfigManager.get_config(environment="nersc")
        assert cfg.infrastructure.output_paths.data_dir == "/data/nersc"

    def test_analysis_group_override(self, use_test_configs):
        cfg = HydraConfigManager.get_config(analysis="m2")
        assert cfg.models.base.threshold == 0.5
        assert cfg.models.base.binning == 2

    def test_hyperparams_group_override(self, use_test_configs):
        cfg = HydraConfigManager.get_config(hyperparams="base_params")
        # top-level keys are LISTS of available values; selected.* is active
        assert list(cfg.hyperparams.co_shift) == [0.0, 0.5, 1.0]
        assert cfg.hyperparams.selected.co_shift == 0.0
        assert cfg.hyperparams.blinded is True

    def test_multiple_group_overrides(self, use_test_configs):
        cfg = HydraConfigManager.get_config(
            environment="nersc",
            analysis="m2",
            hyperparams="base_params",
        )
        assert cfg.infrastructure.output_paths.data_dir == "/data/nersc"
        assert cfg.models.base.threshold == 0.5
        assert cfg.hyperparams.selected.co_shift == 0.0

    def test_selected_value_override(self, use_test_configs):
        # The active value lives at hyperparams.selected.*; override that.
        cfg = HydraConfigManager.get_config(
            hyperparams="base_params",
            overrides=["hyperparams.selected.co_shift=0.5"],
        )
        assert cfg.hyperparams.selected.co_shift == 0.5
        # the list of available values is untouched
        assert list(cfg.hyperparams.co_shift) == [0.0, 0.5, 1.0]

    def test_cli_style_value_override(self, use_test_configs):
        cfg = HydraConfigManager.get_config(
            overrides=["models.base.threshold=42"]
        )
        assert cfg.models.base.threshold == 42

    def test_combined_group_and_value_overrides(self, use_test_configs):
        cfg = HydraConfigManager.get_config(
            analysis="m2",
            overrides=["models.base.binning=99"],
        )
        assert cfg.models.base.threshold == 0.5
        assert cfg.models.base.binning == 99


# ----------------------------------------------------------------------
# Output path defaults
# ----------------------------------------------------------------------
class TestOutputPathDefaults:
    def test_derives_missing_subdirs_from_data_dir(self, use_test_configs):
        cfg = HydraConfigManager.get_config()
        op = cfg.infrastructure.output_paths
        assert op.data_dir == "data"
        assert op.fits_dir == "data/fits"
        assert op.logs_dir == "data/logs"
        assert op.plots_dir == "data/plots"

    def test_respects_explicitly_set_subdirs(self, use_test_configs):
        # nersc env sets fits_dir explicitly
        cfg = HydraConfigManager.get_config(environment="nersc")
        op = cfg.infrastructure.output_paths
        assert op.fits_dir == "/data/nersc/fits"
        # other dirs derived
        assert op.logs_dir == "/data/nersc/logs"
        assert op.plots_dir == "/data/nersc/plots"


# ----------------------------------------------------------------------
# Discovery helpers
# ----------------------------------------------------------------------
class TestListOptions:
    def test_returns_sorted_options(self, use_test_configs):
        opts = HydraConfigManager.list_options()
        assert opts["environments"] == ["base", "nersc"]
        assert opts["analysis"] == ["m2", "models_base"]
        assert opts["hyperparams"] == ["base_params"]

    def test_keys_present_even_when_subdir_missing(self, tmp_path):
        # Empty config dir: subdirs absent -> still returns all keys
        HydraConfigManager.set_config_dir(tmp_path)
        try:
            opts = HydraConfigManager.list_options()
            assert opts == {"environments": [], "analysis": [], "hyperparams": []}
        finally:
            HydraConfigManager.reset_config_dir()

    def test_returns_lists_of_strings(self, use_test_configs):
        opts = HydraConfigManager.list_options()
        for key, value in opts.items():
            assert isinstance(value, list)
            assert all(isinstance(v, str) for v in value), key


# ----------------------------------------------------------------------
# Dict / YAML helpers
# ----------------------------------------------------------------------
class TestConversionHelpers:
    def test_get_config_dict_returns_plain_dict(self, use_test_configs):
        d = HydraConfigManager.get_config_dict()
        assert isinstance(d, dict)
        assert not isinstance(d, DictConfig)
        assert d["pipeline"]["name"] == "test"

    def test_get_config_dict_resolves_nested_structure(self, use_test_configs):
        d = HydraConfigManager.get_config_dict()
        assert isinstance(d["infrastructure"], dict)
        assert isinstance(d["infrastructure"]["output_paths"], dict)

    def test_print_config(self, use_test_configs, capsys):
        HydraConfigManager.print_config()
        captured = capsys.readouterr().out
        assert "pipeline" in captured
        assert "name: test" in captured


# ----------------------------------------------------------------------
# Error handling
# ----------------------------------------------------------------------
class TestErrorHandling:
    def test_missing_config_directory_raises(self, tmp_path):
        missing = tmp_path / "does_not_exist"
        HydraConfigManager.set_config_dir(missing)
        try:
            with pytest.raises(ConfigLoadError) as exc_info:
                HydraConfigManager.get_config()
            assert exc_info.value.config_dir == missing
        finally:
            HydraConfigManager.reset_config_dir()

    def test_invalid_config_name_raises_config_load_error(self, use_test_configs):
        with pytest.raises(ConfigLoadError) as exc_info:
            HydraConfigManager.get_config(config_name="nonexistent")
        assert exc_info.value.config_name == "nonexistent"
        assert exc_info.value.__cause__ is not None

    def test_invalid_override_raises_config_load_error(self, use_test_configs):
        with pytest.raises(ConfigLoadError):
            HydraConfigManager.get_config(overrides=["this.is.not=valid+++"])

    def test_invalid_group_name_raises_config_load_error(self, use_test_configs):
        with pytest.raises(ConfigLoadError):
            HydraConfigManager.get_config(environment="nope_does_not_exist")


# ----------------------------------------------------------------------
# Hydra session isolation
# ----------------------------------------------------------------------
class TestHydraIsolation:
    def test_consecutive_loads_do_not_collide(self, use_test_configs):
        cfg1 = HydraConfigManager.get_config()
        cfg2 = HydraConfigManager.get_config(environment="nersc")
        cfg3 = HydraConfigManager.get_config()
        assert cfg1.infrastructure.output_paths.data_dir == "data"
        assert cfg2.infrastructure.output_paths.data_dir == "/data/nersc"
        assert cfg3.infrastructure.output_paths.data_dir == "data"

    def test_session_cleared_after_error(self, use_test_configs):
        with pytest.raises(ConfigLoadError):
            HydraConfigManager.get_config(config_name="nonexistent")
        # Should still be able to load after a failure
        cfg = HydraConfigManager.get_config()
        assert cfg.pipeline.name == "test"


# ----------------------------------------------------------------------
# Smoke test against real repo config
# ----------------------------------------------------------------------
class TestRealRepoConfig:
    def test_real_default_config_loads(self, use_real_configs):
        cfg = HydraConfigManager.get_config()
        assert isinstance(cfg, DictConfig)
        assert "infrastructure" in cfg

    def test_real_list_options_includes_known_files(self, use_real_configs):
        opts = HydraConfigManager.list_options()
        assert "base" in opts["environments"]
        assert "models_base" in opts["analysis"]
