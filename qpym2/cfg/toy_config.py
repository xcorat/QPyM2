"""
Toy configuration loader for QPyM2.

Uses Hydra compose API to load toy fit configurations from config/toys/.
Provides bias ranges, output settings, and batch parameter grids.

Follows the same pattern as other config tiers (environments/, analysis/, 
hyperparams/): base.yaml is the default, other files are named variants.
"""

from pathlib import Path
from typing import Optional
from hydra import compose, initialize_config_dir
from hydra.core.global_hydra import GlobalHydra
from omegaconf import OmegaConf


def _get_config_dir() -> Path:
    """Get the config directory path."""
    module_dir = Path(__file__).parent
    pkg_root = module_dir.parent.parent
    return pkg_root / "config"


def load_toy_config(config_name: Optional[str] = None) -> dict:
    """
    Load a toy configuration using Hydra compose API.
    
    Args:
        config_name: Name of the config file (without .yaml extension).
                    If None, loads 'base.yaml' (the default).
                    Examples: 'batch', 'tech_note', etc.
    
    Returns:
        Dictionary with toy configuration parameters:
        - bias: {min, max} - signal injection ranges
        - output: {csv, logs, append} - output settings  
        - model_grid: {co_cuts, ndbd_cutoff, co_shift} - parameter grids
    
    Examples:
        >>> cfg = load_toy_config()  # loads base.yaml
        >>> print(cfg['bias']['min'])
        5
        >>> batch_cfg = load_toy_config('batch')  # loads batch.yaml
        >>> print(batch_cfg['model_grid']['co_shift'])
        [0.0, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5]
    """
    if config_name is None:
        config_name = "base"
    
    config_dir = _get_config_dir()
    toys_dir = config_dir / "toys"
    
    # Validate toy config exists
    config_file = toys_dir / f"{config_name}.yaml"
    if not config_file.exists():
        available = [f.stem for f in toys_dir.glob('*.yaml')]
        raise FileNotFoundError(
            f"Toy config file not found: {config_file}\n"
            f"Available configs: {', '.join(available)}"
        )
    
    # Clear any existing Hydra instance
    GlobalHydra.instance().clear()
    
    # Load via Hydra compose with toys override
    with initialize_config_dir(
        version_base=None,
        config_dir=str(config_dir.absolute())
    ):
        cfg = compose(
            config_name="config",
            overrides=[f"toys={config_name}"]
        )
        # Extract and return just the toys section as plain dict
        result = OmegaConf.to_container(cfg.toys, resolve=True)
        # Resolve output paths relative to project root (package root)
        pkg_root = _get_config_dir().parent
        for key in ('csv', 'logs'):
            if key in result.get('output', {}):
                result['output'][key] = str(pkg_root / result['output'][key])
        return result
