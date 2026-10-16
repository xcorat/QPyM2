"""
QPyM2 Configuration System

Hydra-based configuration management for QPyM2.
Provides hierarchical config composition with three tiers:
  - Tier 1: Infrastructure (environments/)
  - Tier 2: Analysis models (analysis/)
  - Tier 3: Hyperparameters (hyperparams/)

Primary interface: HydraConfigManager from hydra_adapter module.
Legacy config code archived in qpym2/cfg_old/ for reference only.

Toy configurations (config/toys/) use simple YAML loading via toy_config module.
Follows the same pattern: base.yaml is default, named variants override.
"""

from .hydra_adapter import HydraConfigManager, get_config
from .config_builders import (
    get_comp_cfg, get_h2uv_model,
    get_data_dir, get_plots_dir, get_syst_dir, get_toy_results_dir,
)
from .toy_config import load_toy_config

__version__ = "2.0.0"
__all__ = [
    "HydraConfigManager", "get_config", "get_comp_cfg", "get_h2uv_model",
    "load_toy_config",
    "get_data_dir", "get_plots_dir", "get_syst_dir", "get_toy_results_dir",
]
