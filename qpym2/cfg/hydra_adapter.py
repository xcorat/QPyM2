"""Hydra Configuration Adapter for QPyM2 (v0.6).

Provides a simple interface for loading Hydra configs across the unified 4-tier
hierarchy: environments, staging, analysis, and models.
Uses Hydra's compose API for programmatic config loading without decorators.

Usage:
    from qpym2.cfg.hydra_adapter import HydraConfigManager, get_config
    
    # Load default config
    cfg = get_config()
    
    # Load with overrides
    cfg = get_config(
        environment="nersc",
        staging="default",
        analysis="base",
        models="tech_note"
    )
"""

import os
import logging
from pathlib import Path
from typing import Optional, List, Tuple, Union

from hydra import compose, initialize_config_dir
from hydra.core.global_hydra import GlobalHydra
from omegaconf import DictConfig, OmegaConf

logger = logging.getLogger(__name__)


def get_repo_root() -> Path:
    """Return the repository root directory (3 levels up from this file)."""
    return Path(__file__).resolve().parent.parent.parent


def resolve_path(path_val: Union[str, Path, None], base: Union[str, Path, None] = None) -> Optional[str]:
    """Resolve a path value against a base directory.

    Rules:
        * ``None`` or ``""`` resolves to ``None``.
        * An absolute path (``Path.is_absolute()``, after ``expanduser``) is
          returned unchanged.
        * Otherwise the path is joined onto ``base`` (or the repository root
          when ``base`` is ``None``). A relative ``base`` is resolved against
          the repository root first.

    Returns the absolute path as a string, or ``None``.
    """
    if path_val is None or path_val == "":
        return None
    p = Path(str(path_val)).expanduser()
    if p.is_absolute():
        return str(p)
    root = Path(str(base)).expanduser() if base is not None else get_repo_root()
    if not root.is_absolute():
        root = get_repo_root() / root
    return str((root / p).resolve())


def resolve_repo_path(path_val: Union[str, Path, None]) -> Optional[str]:
    """Resolve a path relative to the repository root (absolute paths pass through)."""
    return resolve_path(path_val)


def resolve_config_dir(
    config_dir: Union[str, Path, None] = None,
    data_dir: Union[str, Path, None] = None,
) -> Tuple[Path, str]:
    """Resolve which directory holds the Hydra configuration.

    Precedence:
        1. explicit ``config_dir`` argument
        2. ``QPYM2_CONFIG_DIR`` environment variable
        3. ``data_dir`` argument or ``QPYM2_DATA_DIR`` -> ``<dir>/config``
        4. ``./config`` (current working directory); emits a warning
        5. repository ``config/`` (packaged example)

    Returns a ``(path, source)`` tuple where ``source`` is a short label
    describing which rule matched.
    """
    if config_dir is not None:
        return Path(config_dir).expanduser().resolve(), 'flag'

    env_config_dir = os.environ.get('QPYM2_CONFIG_DIR')
    if env_config_dir:
        return Path(env_config_dir).expanduser().resolve(), 'QPYM2_CONFIG_DIR'

    resolved_data_dir = data_dir or os.environ.get('QPYM2_DATA_DIR')
    if resolved_data_dir:
        candidate = Path(resolved_data_dir).expanduser().resolve() / 'config'
        if candidate.is_dir():
            return candidate, 'data_dir'

    cwd_candidate = Path.cwd() / 'config'
    if cwd_candidate.is_dir():
        return cwd_candidate.resolve(), 'cwd'

    return (get_repo_root() / 'config').resolve(), 'repo'


def get_data_root(cfg: DictConfig) -> str:
    """Return the absolute data root from config (or ``QPYM2_DATA_DIR``).

    Raises:
        ValueError: If no data root is set, or it is not an absolute path.
    """
    raw = None
    infrastructure = getattr(cfg, 'infrastructure', None)
    if infrastructure is not None:
        raw = getattr(infrastructure, 'data_root', None)
    if raw is None or raw == "":
        raw = os.environ.get('QPYM2_DATA_DIR')
    if raw is None or raw == "":
        raise ValueError(
            "No data root configured. Set infrastructure.data_root to an absolute "
            "path, or set the QPYM2_DATA_DIR environment variable."
        )
    p = Path(str(raw)).expanduser()
    if not p.is_absolute():
        raise ValueError(
            f"infrastructure.data_root must be an absolute path (got {raw!r}); "
            "all relative paths are resolved against it."
        )
    return str(p)


def resolve_data_path(cfg: DictConfig, path_val: Union[str, Path, None]) -> Optional[str]:
    """Resolve a folder path relative to the configured data root."""
    if path_val is None or path_val == "":
        return None
    return resolve_path(path_val, base=get_data_root(cfg))


class HydraConfigManager:
    """Load Hydra configurations for QPyM2."""
    
    _config_dir: Optional[Path] = None
    
    @classmethod
    def _get_config_dir(cls) -> Path:
        """Get the config directory path."""
        if cls._config_dir is None:
            cls._config_dir = get_repo_root() / "config"
        return cls._config_dir
    
    @classmethod
    def set_config_dir(cls, config_dir: Union[str, Path]) -> None:
        """Override the config directory path."""
        cls._config_dir = Path(config_dir)
    
    @classmethod
    def get_config(
        cls,
        environment: Optional[str] = None,
        staging: Optional[str] = None,
        analysis: Optional[str] = None,
        models: Optional[str] = None,
        model: Optional[str] = None,
        hyperparams: Optional[str] = None,
        overrides: Optional[List[str]] = None,
        config_name: str = "config",
        config_dir: Union[str, Path, None] = None,
        data_dir: Union[str, Path, None] = None,
    ) -> DictConfig:
        """
        Load configuration with optional overrides.

        The configuration directory is resolved in this order:
            1. ``config_dir`` argument
            2. ``QPYM2_CONFIG_DIR`` environment variable
            3. ``data_dir`` / ``QPYM2_DATA_DIR`` -> ``<dir>/config``
            4. ``./config`` (current working directory); emits a warning
            5. repository ``config/`` (packaged example)
        """
        if cls._config_dir is not None:
            resolved_config_dir = Path(cls._config_dir).expanduser().resolve()
            config_source = 'set_config_dir'
        else:
            resolved_config_dir, config_source = resolve_config_dir(
                config_dir=config_dir, data_dir=data_dir
            )
            if config_source == 'cwd':
                logger.warning(
                    "Using './config' as the configuration directory. Pass "
                    "--config-dir or set QPYM2_CONFIG_DIR to silence this warning."
                )

        config_dir = resolved_config_dir

        if not config_dir.exists():
            raise FileNotFoundError(f"Config directory not found: {config_dir}")
        
        override_list = []
        if environment is not None:
            override_list.append(f"environments={environment}")
        if staging is not None:
            override_list.append(f"staging={staging}")
        if analysis is not None:
            override_list.append(f"analysis={analysis}")
        resolved_models = models or model
        if resolved_models is not None:
            override_list.append(f"models={resolved_models}")
        if hyperparams is not None:
            logger.debug("hyperparams argument is deprecated and ignored.")
        
        if overrides:
            override_list.extend(overrides)
        
        GlobalHydra.instance().clear()
        
        try:
            with initialize_config_dir(
                version_base=None, 
                config_dir=str(config_dir.absolute())
            ):
                cfg = compose(config_name=config_name, overrides=override_list)
                logger.debug(f"Loaded config with overrides: {override_list}")
                return cfg
        except Exception as e:
            logger.error(f"Failed to load config: {e}")
            raise

    @classmethod
    def get_config_dict(cls, **kwargs) -> dict:
        """Load configuration and convert to plain Python dict."""
        cfg = cls.get_config(**kwargs)
        return OmegaConf.to_container(cfg, resolve=True)
    
    @classmethod
    def print_config(cls, cfg: Optional[DictConfig] = None, **kwargs) -> None:
        """Print configuration in YAML format."""
        if cfg is None:
            cfg = cls.get_config(**kwargs)
        print(OmegaConf.to_yaml(cfg))

    @classmethod
    def list_options(cls) -> dict:
        """List available config options from the config directory."""
        config_dir = cls._get_config_dir()

        def _list_subdir(subdir_name: str):
            path = config_dir / subdir_name
            if not path.exists() or not path.is_dir():
                return []
            opts = [p.stem for p in sorted(path.glob('*.yaml'))]
            return opts

        return {
            'environments': _list_subdir('environments'),
            'staging': _list_subdir('staging'),
            'analysis': _list_subdir('analysis'),
            'models': _list_subdir('models'),
        }

    @classmethod
    def get_default_selections(cls) -> dict:
        """Read the base `config.yaml` defaults and return default selections."""
        config_file = cls._get_config_dir() / "config.yaml"
        defaults = {
            'environment': None,
            'staging': None,
            'analysis': None,
            'models': None,
        }
        if not config_file.exists():
            return defaults

        try:
            cfg = OmegaConf.load(str(config_file))
            for entry in cfg.get('defaults', []):
                if isinstance(entry, (dict, DictConfig)):
                    for k, v in dict(entry).items():
                        if k == 'environments':
                            defaults['environment'] = v
                        elif k == 'staging':
                            defaults['staging'] = v
                        elif k == 'analysis':
                            defaults['analysis'] = v
                        elif k == 'models':
                            defaults['models'] = v
        except Exception:
            pass

        return defaults


def get_config(
    environment: Optional[str] = None,
    staging: Optional[str] = None,
    analysis: Optional[str] = None,
    models: Optional[str] = None,
    model: Optional[str] = None,
    hyperparams: Optional[str] = None,
    overrides: Optional[List[str]] = None,
    config_name: str = "config",
    config_dir: Union[str, Path, None] = None,
    data_dir: Union[str, Path, None] = None,
) -> DictConfig:
    """Shortcut for HydraConfigManager.get_config()."""
    return HydraConfigManager.get_config(
        environment=environment,
        staging=staging,
        analysis=analysis,
        models=models or model,
        hyperparams=hyperparams,
        overrides=overrides,
        config_name=config_name,
        config_dir=config_dir,
        data_dir=data_dir,
    )
