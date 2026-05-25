"""
Hydra Configuration Adapter for QPyM2.

Provides a simple interface for loading Hydra configs in existing code.
Uses Hydra's compose API for programmatic config loading without decorators.

Example:
    Basic usage::

        from qpym2.cfg import HydraConfigManager, get_config

        # Load default config
        cfg = HydraConfigManager.get_config()

        # Load with overrides. NOTE: top-level keys in
        # config/hyperparams/base_params.yaml (co_shift, co_cuts,
        # ndbd_cutoff, ndbd_shift) are *lists of available values*. The
        # active scalar lives under `hyperparams.selected.*` (see
        # config/hyperparams/m2uv.yaml). Override the active value via:
        cfg = HydraConfigManager.get_config(
            environment="nersc",
            analysis="m2",
            hyperparams="m2uv",
            overrides=["hyperparams.selected.co_shift=0.5"],
        )

        # Access config values
        staging_dir = cfg.infrastructure.output_paths.staging_dir
        active_co_shift = cfg.hyperparams.selected.co_shift
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional, Union

from hydra import compose, initialize_config_dir
from hydra.core.global_hydra import GlobalHydra
from omegaconf import DictConfig, OmegaConf

logger = logging.getLogger(__name__)


__all__ = [
    "HydraConfigManager",
    "ConfigLoadError",
    "get_config",
]


class ConfigLoadError(Exception):
    """Raised when a Hydra configuration cannot be loaded.

    Wraps the underlying exception with contextual information (config
    directory, config name, overrides) so callers can diagnose failures
    without inspecting log output.
    """

    def __init__(
        self,
        message: str,
        *,
        config_dir: Optional[Path] = None,
        config_name: Optional[str] = None,
        overrides: Optional[list[str]] = None,
        cause: Optional[BaseException] = None,
    ) -> None:
        details = []
        if config_name is not None:
            details.append(f"config_name={config_name!r}")
        if config_dir is not None:
            details.append(f"config_dir={str(config_dir)!r}")
        if overrides is not None:
            details.append(f"overrides={overrides!r}")
        suffix = f" ({', '.join(details)})" if details else ""
        super().__init__(f"{message}{suffix}")
        self.config_dir = config_dir
        self.config_name = config_name
        self.overrides = overrides
        self.__cause__ = cause


class HydraConfigManager:
    """Load Hydra configurations for QPyM2.

    Uses Hydra's ``compose`` API for programmatic config loading without
    decorators. Supports composition of multiple config groups
    (``environments``, ``analysis``, ``hyperparams``) with CLI-style
    overrides.
    """

    _config_dir: Optional[Path] = None

    # ------------------------------------------------------------------
    # Config directory management
    # ------------------------------------------------------------------
    @classmethod
    def _get_config_dir(cls) -> Path:
        """Return the configured config directory (defaults to ``<repo>/config``)."""
        if cls._config_dir is None:
            cls._config_dir = Path(__file__).parent.parent.parent / "config"
        return cls._config_dir

    @classmethod
    def set_config_dir(cls, config_dir: Union[str, Path]) -> None:
        """Override the config directory path.

        Args:
            config_dir: Path to a directory containing ``config.yaml`` and
                config groups (``environments/``, ``analysis/``,
                ``hyperparams/``).
        """
        cls._config_dir = Path(config_dir)

    @classmethod
    def reset_config_dir(cls) -> None:
        """Reset the config directory to the default (``<repo>/config``)."""
        cls._config_dir = None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    @classmethod
    def _apply_output_paths_defaults(cls, cfg: DictConfig) -> None:
        """Fill in derived ``output_paths`` defaults from ``data_dir``.

        If ``fits_dir``, ``logs_dir``, or ``plots_dir`` is missing or empty
        and ``data_dir`` is set, derive them as ``<data_dir>/<name>``.

        Args:
            cfg: Hydra ``DictConfig`` to modify in-place.
        """
        if not hasattr(cfg, "infrastructure") or not hasattr(
            cfg.infrastructure, "output_paths"
        ):
            return

        output_paths = cfg.infrastructure.output_paths
        data_dir = output_paths.get("data_dir", None)
        if not data_dir:
            return

        # Allow inserting new keys even when the parent node is in
        # struct mode (e.g. when fits_dir/logs_dir/plots_dir are absent
        # from the loaded environment file).
        prev_struct = OmegaConf.is_struct(output_paths)
        try:
            OmegaConf.set_struct(output_paths, False)
            for key in ("fits_dir", "logs_dir", "plots_dir"):
                current = output_paths.get(key, None)
                if current is None or current == "":
                    # fits_dir -> fits, logs_dir -> logs, plots_dir -> plots
                    derived = f"{data_dir}/{key.removesuffix('_dir')}"
                    output_paths[key] = derived
                    logger.debug("Applied default %s: %s", key, derived)
        finally:
            OmegaConf.set_struct(output_paths, prev_struct)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    @classmethod
    def get_config(
        cls,
        environment: Optional[str] = None,
        analysis: Optional[str] = None,
        hyperparams: Optional[str] = None,
        overrides: Optional[list[str]] = None,
        config_name: str = "config",
    ) -> DictConfig:
        """Load configuration with optional overrides.

        Args:
            environment: Environment config name (e.g. ``base``, ``nersc``).
                If ``None``, uses the default from ``config.yaml``.
            analysis: Analysis config name (e.g. ``models_base``, ``m2``).
                If ``None``, uses the default from ``config.yaml``.
            hyperparams: Hyperparameter config name (e.g. ``base_params``).
                If ``None``, uses the default from ``config.yaml``.
            overrides: Additional CLI-style overrides
                (e.g. ``["hyperparams.selected.co_shift=0.5"]``).
            config_name: Master config file name (default ``"config"``).

        Returns:
            Merged ``OmegaConf`` ``DictConfig``.

        Raises:
            ConfigLoadError: If the config directory is missing or the
                composition fails (invalid config name, malformed YAML,
                bad override, etc.).
        """
        config_dir = cls._get_config_dir()

        if not config_dir.exists():
            raise ConfigLoadError(
                "Config directory not found",
                config_dir=config_dir,
                config_name=config_name,
            )

        override_list: list[str] = []
        if environment is not None:
            override_list.append(f"environments={environment}")
        if analysis is not None:
            override_list.append(f"analysis={analysis}")
        if hyperparams is not None:
            # `hyperparams` has no entry in the master defaults list, so use
            # Hydra's append syntax (+) rather than the override syntax.
            override_list.append(f"+hyperparams={hyperparams}")
        if overrides:
            override_list.extend(overrides)

        try:
            GlobalHydra.instance().clear()
            with initialize_config_dir(
                version_base=None,
                config_dir=str(config_dir.absolute()),
            ):
                cfg = compose(config_name=config_name, overrides=override_list)
            logger.debug("Loaded config with overrides: %s", override_list)
            cls._apply_output_paths_defaults(cfg)
            return cfg
        except ConfigLoadError:
            raise
        except Exception as exc:
            logger.error(
                "Failed to load config '%s' from %s with overrides %s: %s",
                config_name,
                config_dir,
                override_list,
                exc,
                exc_info=True,
            )
            raise ConfigLoadError(
                "Failed to compose Hydra configuration",
                config_dir=config_dir,
                config_name=config_name,
                overrides=override_list,
                cause=exc,
            ) from exc

    @classmethod
    def get_config_dict(cls, **kwargs: Any) -> dict[str, Any]:
        """Load configuration and convert to a plain Python ``dict``.

        Interpolations are resolved.

        Args:
            **kwargs: Forwarded to :meth:`get_config`.

        Returns:
            Plain Python ``dict`` (not ``OmegaConf``).
        """
        cfg = cls.get_config(**kwargs)
        return OmegaConf.to_container(cfg, resolve=True)  # type: ignore[return-value]

    @classmethod
    def print_config(
        cls,
        cfg: Optional[DictConfig] = None,
        **kwargs: Any,
    ) -> None:
        """Print configuration in YAML format.

        Args:
            cfg: Existing config to print, or ``None`` to load fresh.
            **kwargs: Forwarded to :meth:`get_config` if ``cfg`` is ``None``.
        """
        if cfg is None:
            cfg = cls.get_config(**kwargs)
        print(OmegaConf.to_yaml(cfg))

    @classmethod
    def list_options(cls) -> dict[str, list[str]]:
        """List available config option names per config group.

        Returns:
            Dict with keys ``"environments"``, ``"analysis"``,
            ``"hyperparams"``, each mapping to a sorted list of config
            names (without file extension).
        """
        config_dir = cls._get_config_dir()

        def _list_subdir(subdir_name: str) -> list[str]:
            path = config_dir / subdir_name
            if not path.exists() or not path.is_dir():
                return []
            return [p.stem for p in sorted(path.glob("*.yaml"))]

        return {
            "environments": _list_subdir("environments"),
            "analysis": _list_subdir("analysis"),
            "hyperparams": _list_subdir("hyperparams"),
        }


def get_config(**kwargs: Any) -> DictConfig:
    """Shortcut for :meth:`HydraConfigManager.get_config`."""
    return HydraConfigManager.get_config(**kwargs)
