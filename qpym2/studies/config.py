"""
Study configuration loader for QPyM2 (v0.6).

Discovers study directories, loads study-specific configuration files,
inherits default configurations from the 4-tier Hydra system, and
merges study-level overrides using OmegaConf.
"""

from pathlib import Path
from typing import List, Optional, Tuple, Union
import logging
from omegaconf import DictConfig, OmegaConf

from qpym2.cfg.hydra_adapter import HydraConfigManager, get_repo_root
from qpym2.cfg.config_builders import get_fit_data_dir

logger = logging.getLogger(__name__)

STUDY_CONFIG_NAMES = (
    "study.yaml",
    "study.yml",
    "config.yaml",
    "config.yml",
)


def resolve_study_dir(study_dir_or_name: Union[str, Path]) -> Path:
    """
    Resolve a study directory path from a path or named identifier.

    Candidate search locations:
    1. Direct path (absolute or relative to current working directory)
    2. Path relative to repository root
    3. <repo_root>/studies/<name>
    4. <repo_root>/config/studies/<name>
    5. <fit_data_root>/studies/<name>

    Args:
        study_dir_or_name: Path or name of the study.

    Returns:
        Resolved absolute Path to study directory.

    Raises:
        FileNotFoundError: If the study directory cannot be located.
    """
    p = Path(study_dir_or_name).expanduser()

    # 1. Direct path exists
    if p.is_dir():
        return p.resolve()

    # 2. Relative to repo root
    repo_root = get_repo_root()
    repo_rel = repo_root / p
    if repo_rel.is_dir():
        return repo_rel.resolve()

    name = str(study_dir_or_name).strip()

    # 3. repo/studies/<name>
    cand_studies = repo_root / "studies" / name
    if cand_studies.is_dir():
        return cand_studies.resolve()

    # 4. repo/config/studies/<name>
    cand_cfg_studies = repo_root / "config" / "studies" / name
    if cand_cfg_studies.is_dir():
        return cand_cfg_studies.resolve()

    # 5. fit_data_dir/studies/<name>
    try:
        base_cfg = HydraConfigManager.get_config()
        fit_data_root = Path(get_fit_data_dir(base_cfg))
        cand_fit_studies = fit_data_root / "studies" / name
        if cand_fit_studies.is_dir():
            return cand_fit_studies.resolve()
    except Exception:
        pass

    raise FileNotFoundError(
        f"Study directory not found for '{study_dir_or_name}'. "
        f"Checked: '{p}', '{repo_rel}', '{cand_studies}', '{cand_cfg_studies}'."
    )


def find_study_config_file(study_dir: Path) -> Path:
    """
    Locate the primary study YAML file within a study directory.

    Args:
        study_dir: Path to study directory.

    Returns:
        Path to primary study YAML file.

    Raises:
        FileNotFoundError: If no configuration file is found.
    """
    for name in STUDY_CONFIG_NAMES:
        cand = study_dir / name
        if cand.is_file():
            return cand

    # Check for any YAML file in study_dir if none of standard names matched
    yaml_files = sorted(list(study_dir.glob("*.yaml")) + list(study_dir.glob("*.yml")))
    if yaml_files:
        return yaml_files[0]

    raise FileNotFoundError(
        f"No study configuration file ({', '.join(STUDY_CONFIG_NAMES)}) found in '{study_dir}'."
    )


def load_study_config(
    study_dir_or_name: Union[str, Path],
    overrides: Optional[List[str]] = None,
) -> Tuple[DictConfig, Path]:
    """
    Load a study configuration, inheriting defaults and merging overrides.

    Workflow:
    1. Resolves the study directory and locates primary configuration file.
    2. Reads the study config YAML.
    3. Inspects 'defaults' section to determine base configuration tiers
       (environment, staging, analysis, models).
    4. Loads base configuration using HydraConfigManager.
    5. Disables Hydra struct locking (OmegaConf.set_struct) to allow custom keys.
    6. Loads any auxiliary YAML files in study folder as intermediate overrides.
    7. Merges primary study configuration on top of base defaults.
    8. Applies optional command-line dotlist overrides.

    Args:
        study_dir_or_name: Directory path or registered study name.
        overrides: Optional list of dotlist overrides (e.g. ['reference.co_shift=0.5']).

    Returns:
        Tuple of (merged_DictConfig, resolved_study_dir).
    """
    study_dir = resolve_study_dir(study_dir_or_name)
    primary_cfg_file = find_study_config_file(study_dir)

    raw_study_cfg = OmegaConf.load(primary_cfg_file)

    # 1. Extract base configuration defaults
    defaults = raw_study_cfg.get("defaults", {})
    env_name = defaults.get("environment") if defaults else None
    staging_name = defaults.get("staging") if defaults else None
    analysis_name = defaults.get("analysis") if defaults else None
    models_name = defaults.get("models") or (defaults.get("model") if defaults else None)

    # 2. Load base Hydra configuration
    base_cfg = HydraConfigManager.get_config(
        environment=env_name,
        staging=staging_name,
        analysis=analysis_name,
        models=models_name,
    )

    # Allow custom keys in base configuration
    OmegaConf.set_struct(base_cfg, False)

    # 3. Check for any auxiliary YAML files in study_dir (excluding primary)
    for yaml_path in sorted(study_dir.glob("*.yaml")):
        if yaml_path.resolve() != primary_cfg_file.resolve():
            aux_cfg = OmegaConf.load(yaml_path)
            base_cfg = OmegaConf.merge(base_cfg, aux_cfg)

    # Also check study_dir / "config" subfolder if it exists
    aux_cfg_dir = study_dir / "config"
    if aux_cfg_dir.is_dir():
        for yaml_path in sorted(aux_cfg_dir.glob("*.yaml")):
            aux_cfg = OmegaConf.load(yaml_path)
            base_cfg = OmegaConf.merge(base_cfg, aux_cfg)

    # 4. Merge primary study config on top
    merged_cfg = OmegaConf.merge(base_cfg, raw_study_cfg)

    # 5. Apply any CLI/runtime overrides
    if overrides:
        cli_overrides = OmegaConf.from_dotlist(overrides)
        merged_cfg = OmegaConf.merge(merged_cfg, cli_overrides)

    # Ensure study metadata fields exist
    if "study" not in merged_cfg:
        merged_cfg.study = OmegaConf.create({})
    if "name" not in merged_cfg.study or not merged_cfg.study.name:
        merged_cfg.study.name = study_dir.name
    merged_cfg.study.study_dir = str(study_dir)

    return merged_cfg, study_dir
