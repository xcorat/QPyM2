"""Pytest fixtures for QPyM2 config tests."""

from __future__ import annotations

from pathlib import Path

import pytest
from hydra.core.global_hydra import GlobalHydra

from qpym2.cfg.hydra_adapter import HydraConfigManager


FIXTURES_DIR = Path(__file__).parent / "fixtures" / "test_configs"


@pytest.fixture(autouse=True)
def _clear_hydra():
    """Ensure no stale GlobalHydra state between tests."""
    GlobalHydra.instance().clear()
    yield
    GlobalHydra.instance().clear()


@pytest.fixture
def test_config_dir() -> Path:
    """Path to the bundled minimal test config directory."""
    return FIXTURES_DIR


@pytest.fixture
def use_test_configs(test_config_dir: Path):
    """Point HydraConfigManager at the test fixture directory."""
    HydraConfigManager.set_config_dir(test_config_dir)
    yield test_config_dir
    HydraConfigManager.reset_config_dir()


@pytest.fixture
def use_real_configs():
    """Point HydraConfigManager at the real repo `config/` directory."""
    real = Path(__file__).resolve().parents[3] / "config"
    HydraConfigManager.set_config_dir(real)
    yield real
    HydraConfigManager.reset_config_dir()
