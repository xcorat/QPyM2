"""
QPyM2 Configuration System (v0.6)

Hydra-based configuration management for QPyM2.
Provides decoupled hierarchical config composition across four tiers:
  - Tier 1: Infrastructure (environments/)
  - Tier 2: Staging data contract (staging/)
  - Tier 3: Analysis baseline and physics cuts (analysis/)
  - Tier 4: Analysis models, ROI, and binning (models/)

Primary interface: HydraConfigManager and get_config from hydra_adapter.
"""

from .hydra_adapter import (
    HydraConfigManager,
    get_config,
    resolve_data_path,
    resolve_path,
    resolve_repo_path,
    get_data_root,
    get_cuore_data_root,
    get_repo_root,
)
from .config_builders import (
    build_staging_config,
    build_input_config,
    build_mc_config,
    build_ndbd_config,
    build_data_config,
    build_write_config,
    get_comp_cfg,
    get_h2uv_model,
    get_jags_path,
    get_data_names,
    get_staging_dir,
    get_bkg_mc_dir,
    get_ndbd_mc_dir,
    get_fit_dir,
    get_data_dir,
    get_plots_dir,
    get_syst_dir,
    get_toy_results_dir,
    get_artifact_store_dir,
    get_fit_data_dir,
    get_model_storage_coords,
    get_input_fname,
    get_corrmat_path,
)
from .toy_config import load_toy_config

__version__ = "0.6.0"
__all__ = [
    "HydraConfigManager",
    "get_config",
    "resolve_data_path",
    "resolve_path",
    "resolve_repo_path",
    "get_data_root",
    "get_cuore_data_root",
    "get_repo_root",
    "build_staging_config",
    "build_input_config",
    "build_mc_config",
    "build_ndbd_config",
    "build_data_config",
    "build_write_config",
    "get_comp_cfg",
    "get_h2uv_model",
    "get_jags_path",
    "get_data_names",
    "get_staging_dir",
    "get_bkg_mc_dir",
    "get_ndbd_mc_dir",
    "get_fit_dir",
    "get_data_dir",
    "get_plots_dir",
    "get_syst_dir",
    "get_toy_results_dir",
    "get_artifact_store_dir",
    "get_fit_data_dir",
    "get_model_storage_coords",
    "get_input_fname",
    "get_corrmat_path",
    "load_toy_config",
]
