"""
Hydra Configuration Adapter for QPyM2.

Provides a simple interface for loading Hydra configs in existing code.
Uses Hydra's compose API for programmatic config loading without decorators.

Usage:
    from qpym2.cfg.hydra_adapter import HydraConfigManager
    
    # Load default config
    cfg = HydraConfigManager.get_config()
    
    # Load with overrides
    cfg = HydraConfigManager.get_config(
        environment="nersc",
        analysis="tech_note",
        hyperparams="bias_study"
    )
    
    # Access config values
    staging_dir = cfg.infrastructure.output_paths.staging_dir
    threshold = cfg.models.tech_note_200.roi.threshold
    co_shift = cfg.hyperparams.co_shift
"""

from pathlib import Path
from typing import Optional, Any, List, Union
from contextlib import contextmanager
import logging

from hydra import compose, initialize_config_dir
from hydra.core.global_hydra import GlobalHydra
from omegaconf import DictConfig, OmegaConf

logger = logging.getLogger(__name__)


class HydraConfigManager:
    """Load Hydra configurations for QPyM2."""
    
    _config_dir: Optional[Path] = None
    
    @classmethod
    def _get_config_dir(cls) -> Path:
        """Get the config directory path."""
        if cls._config_dir is None:
            # Default: config/ relative to project root (3 levels up from this file)
            cls._config_dir = Path(__file__).parent.parent.parent / "config"
        return cls._config_dir
    
    @classmethod
    def set_config_dir(cls, config_dir: Union[str, Path]) -> None:
        """Override the config directory path."""
        cls._config_dir = Path(config_dir)
    
    @classmethod
    def _apply_output_paths_defaults(cls, cfg: DictConfig) -> None:
        """
        Apply defaults to infrastructure.output_paths.
        
        If plots_dir is not set (None or empty string), set it to fit_dir + '/plots'.
        
        Args:
            cfg: Hydra DictConfig object to modify in-place
        """
        if not hasattr(cfg, 'infrastructure') or not hasattr(cfg.infrastructure, 'output_paths'):
            return
        
        output_paths = cfg.infrastructure.output_paths
        
        # Check if plots_dir needs a default
        plots_dir = output_paths.get('plots_dir', None) if output_paths else None
        fit_dir = output_paths.get('fit_dir', None) if output_paths else None
        
        # If plots_dir is not set and fit_dir exists, use fit_dir + '/plots'
        if (plots_dir is None or plots_dir == '') and fit_dir:
            output_paths.plots_dir = f"{fit_dir}/plots"
            logger.debug(f"Applied default plots_dir: {output_paths.plots_dir}")
    
    @classmethod
    def get_config(
        cls,
        environment: Optional[str] = None,
        analysis: Optional[str] = None,
        hyperparams: Optional[str] = None,
        overrides: Optional[List[str]] = None,
        config_name: str = "config"
    ) -> DictConfig:
        """
        Load configuration with optional overrides.
        
        Args:
            environment: Environment config name (base, development, nersc, tech_note, production)
                        If None, uses default from config.yaml
            analysis: Analysis config name (models_base, tech_note, threshold_350)
                     If None, uses default from config.yaml
            hyperparams: Hyperparameter config name (base, bias_study)
                        If None, uses default from config.yaml
            overrides: Additional CLI-style overrides (e.g., ["hyperparams.co_shift=0.5"])
            config_name: Master config file name (default: "config")
        
        Returns:
            OmegaConf DictConfig object with merged configuration
        
        Example:
            cfg = HydraConfigManager.get_config(
                environment="nersc",
                analysis="tech_note",
                overrides=["hyperparams.co_shift=0.5", "hyperparams.co_cuts=v150"]
            )
        """
        
        config_dir = cls._get_config_dir()
        
        if not config_dir.exists():
            raise FileNotFoundError(f"Config directory not found: {config_dir}")
        
        # Build overrides list - only add if explicitly specified
        override_list = []
        if environment is not None:
            override_list.append(f"environments={environment}")
        if analysis is not None:
            override_list.append(f"analysis={analysis}")
        if hyperparams is not None:
            override_list.append(f"hyperparams={hyperparams}")
        
        if overrides:
            override_list.extend(overrides)
        
        # Clear any existing Hydra instance
        GlobalHydra.instance().clear()
        
        try:
            with initialize_config_dir(
                version_base=None, 
                config_dir=str(config_dir.absolute())
            ):
                cfg = compose(config_name=config_name, overrides=override_list)
                logger.debug(f"Loaded config with overrides: {override_list}")
                
                # Apply default for plots_dir if not explicitly set
                cls._apply_output_paths_defaults(cfg)
                
                return cfg
        except Exception as e:
            logger.error(f"Failed to load config: {e}")
            raise

    @classmethod
    def list_options(cls) -> dict:
        """
        List available configuration option names for environments, analysis and hyperparams.

        Returns:
            Dict with keys 'environments', 'analysis', 'hyperparams' mapping to lists of names (without extension).
        """
        config_dir = cls._get_config_dir()
        options = {}
        for sub in ["environments", "analysis", "hyperparams"]:
            subdir = config_dir / sub
            names = []
            if subdir.exists():
                for f in sorted(subdir.iterdir()):
                    if f.is_file() and f.suffix in ('.yaml', '.yml'):
                        names.append(f.stem)
            options[sub] = names
        return options

    @classmethod
    def get_default_selections(cls) -> dict:
        """
        Return the selections defined in `config/config.yaml` defaults.

        Returns:
            Dict with keys 'environment', 'analysis', 'hyperparams' and values as the default name strings or None.
        """
        config_dir = cls._get_config_dir()
        config_file = config_dir / "config.yaml"
        result = {"environment": None, "analysis": None, "hyperparams": None}
        if not config_file.exists():
            return result

        try:
            from omegaconf import OmegaConf
            cfg = OmegaConf.load(str(config_file))
            defaults = cfg.get('defaults', [])
            for item in defaults:
                # defaults entries can be dicts like {'environments': 'base'} or strings; handle dicts
                if isinstance(item, dict):
                    for k, v in item.items():
                        if k == 'environments':
                            result['environment'] = v
                        elif k == 'analysis':
                            result['analysis'] = v
                        elif k == 'hyperparams':
                            result['hyperparams'] = v
            return result
        except Exception:
            return result
    
    @classmethod
    def get_config_dict(cls, **kwargs) -> dict:
        """
        Load configuration and convert to plain Python dict.
        
        Same arguments as get_config().
        
        Returns:
            Plain Python dictionary (not OmegaConf)
        """
        cfg = cls.get_config(**kwargs)
        return OmegaConf.to_container(cfg, resolve=True)
    
    @classmethod
    def print_config(cls, cfg: Optional[DictConfig] = None, **kwargs) -> None:
        """
        Print configuration in YAML format.
        
        Args:
            cfg: Existing config to print, or None to load fresh
            **kwargs: Arguments passed to get_config() if cfg is None
        """
        if cfg is None:
            cfg = cls.get_config(**kwargs)
        print(OmegaConf.to_yaml(cfg))

    @classmethod
    def list_options(cls) -> dict:
        """
        List available config options from the config directory.

        Returns:
            A dictionary with keys 'environments', 'analysis', and 'hyperparams',
            each containing a sorted list of available config names (without file extension).
        """
        config_dir = cls._get_config_dir()

        def _list_subdir(subdir_name: str):
            path = config_dir / subdir_name
            if not path.exists() or not path.is_dir():
                return []
            opts = [p.stem for p in sorted(path.glob('*.yaml'))]
            return opts

        return {
            'environments': _list_subdir('environments'),
            'analysis': _list_subdir('analysis'),
            'hyperparams': _list_subdir('hyperparams')
        }

    @classmethod
    def get_default_selections(cls) -> dict:
        """
        Read the base `config.yaml` defaults and return the default selections used by Hydra.

        Returns:
            A dict with keys 'environment', 'analysis', 'hyperparams' mapping to their default names
            (or None if not specified).
        """
        config_file = cls._get_config_dir() / "config.yaml"
        defaults = {'environment': None, 'analysis': None, 'hyperparams': None}
        if not config_file.exists():
            return defaults

        try:
            cfg = OmegaConf.load(str(config_file))
            for entry in cfg.get('defaults', []):
                # entries may be plain dicts or OmegaConf DictConfig
                if isinstance(entry, (dict, DictConfig)):
                    for k, v in dict(entry).items():
                        if k == 'environments':
                            defaults['environment'] = v
                        elif k == 'analysis':
                            defaults['analysis'] = v
                        elif k == 'hyperparams':
                            defaults['hyperparams'] = v
        except Exception:
            pass

        return defaults


# Convenience function for quick access
def get_config(**kwargs) -> DictConfig:
    """Shortcut for HydraConfigManager.get_config()."""
    return HydraConfigManager.get_config(**kwargs)
