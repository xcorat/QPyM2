"""
Study assembly routines for QPyM2 (v0.6).

Builds a complete Study instance by expanding configuration parameters,
resolving parameter sweeps (Cartesian products), and determining
StudyPoint specifications.
"""

import itertools
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from omegaconf import OmegaConf

from qpym2.studies.config import load_study_config
from qpym2.studies.point import StudyPoint


_KEY_ALIASES = {
    "co_shifts": "co_shift",
    "co_cuts": "co_cut",
    "ndbd_cutoffs": "ndbd_cutoff",
    "ndbd_shifts": "ndbd_shift",
    "model_keys": "model_key",
    "models": "model_key",
    "model": "model_key",
}


def _canonical_param_key(key: str) -> str:
    """Map plural or alternate parameter names to canonical singular names."""
    return _KEY_ALIASES.get(key, key)


def _format_point_id(task_id: int, swept_params: Dict[str, Any]) -> str:
    """Format a clean, descriptive point_id based on swept parameters."""
    if not swept_params:
        return "nominal"

    tokens = []
    for k, v in swept_params.items():
        if k == "co_shift":
            tokens.append(f"co_shift_{v}")
        elif k == "ndbd_shift":
            tokens.append(f"ndbd_shift_{v}")
        elif k == "co_cut":
            tokens.append(f"co_{v}")
        elif k == "ndbd_cutoff":
            tokens.append(f"ndbd_{v}")
        elif k == "model_key":
            tokens.append(f"model_{v}")
        else:
            tokens.append(f"{k}_{v}")

    return f"p{task_id}_" + "_".join(tokens)


def expand_study_points(
    reference: Dict[str, Any],
    sweep: Optional[Dict[str, Any]] = None,
    default_model_key: str = "default",
) -> List[StudyPoint]:
    """
    Expand reference parameters and parameter sweep into a list of StudyPoints.

    Args:
        reference: Dictionary of baseline parameter values.
        sweep: Optional dictionary mapping parameter names to lists of values to sweep.
        default_model_key: Default model identifier if not in reference.

    Returns:
        List of initialized StudyPoint instances.
    """
    ref_dict = dict(reference or {})

    base_model_key = str(ref_dict.get("model_key", default_model_key))
    base_co_cut = ref_dict.get("co_cut")
    base_ndbd_cutoff = ref_dict.get("ndbd_cutoff")
    base_co_shift = float(ref_dict.get("co_shift", 0.0))
    base_ndbd_shift = float(ref_dict.get("ndbd_shift", 0.0))
    base_blinded = bool(ref_dict.get("blinded", False))

    if not sweep:
        # Single-point baseline study
        return [
            StudyPoint(
                point_id="nominal",
                task_id=0,
                model_key=base_model_key,
                co_cut=base_co_cut,
                ndbd_cutoff=base_ndbd_cutoff,
                co_shift=base_co_shift,
                ndbd_shift=base_ndbd_shift,
                blinded=base_blinded,
            )
        ]

    # Normalize sweep keys and ensure lists
    sweep_axes: List[Tuple[str, List[Any]]] = []
    for raw_k, vals in sweep.items():
        k = _canonical_param_key(raw_k)
        if isinstance(vals, (list, tuple)):
            sweep_axes.append((k, list(vals)))
        else:
            sweep_axes.append((k, [vals]))

    if not sweep_axes:
        return [
            StudyPoint(
                point_id="nominal",
                task_id=0,
                model_key=base_model_key,
                co_cut=base_co_cut,
                ndbd_cutoff=base_ndbd_cutoff,
                co_shift=base_co_shift,
                ndbd_shift=base_ndbd_shift,
                blinded=base_blinded,
            )
        ]

    axis_names = [ax[0] for ax in sweep_axes]
    axis_values = [ax[1] for ax in sweep_axes]

    points: List[StudyPoint] = []
    for task_id, combo in enumerate(itertools.product(*axis_values)):
        combo_dict = dict(zip(axis_names, combo))

        m_key = str(combo_dict.get("model_key", base_model_key))
        c_cut = combo_dict.get("co_cut", base_co_cut)
        n_cutoff = combo_dict.get("ndbd_cutoff", base_ndbd_cutoff)
        c_shift = float(combo_dict.get("co_shift", base_co_shift))
        n_shift = float(combo_dict.get("ndbd_shift", base_ndbd_shift))
        blind = bool(combo_dict.get("blinded", base_blinded))

        point_id = _format_point_id(task_id, combo_dict)

        points.append(
            StudyPoint(
                point_id=point_id,
                task_id=task_id,
                model_key=m_key,
                co_cut=c_cut,
                ndbd_cutoff=n_cutoff,
                co_shift=c_shift,
                ndbd_shift=n_shift,
                blinded=blind,
                params=combo_dict,
            )
        )

    return points


def assemble_study(
    study_dir_or_name: Union[str, Path],
    overrides: Optional[List[str]] = None,
) -> Any:
    """
    Assemble a Study instance from a study directory and configuration.

    Workflow:
    1. Loads merged Hydra + study configuration.
    2. Resolves baseline/reference parameters and sweep axes.
    3. Expands the parameter combinations into StudyPoints.
    4. Instantiates and returns the Study object.

    Args:
        study_dir_or_name: Path or registered name of study directory.
        overrides: Optional dotlist overrides.

    Returns:
        Study instance.
    """
    # Lazy import to avoid circular dependency
    from qpym2.studies.study import Study

    cfg, study_dir = load_study_config(study_dir_or_name, overrides=overrides)

    study_sec = cfg.get("study", {})
    study_name = str(study_sec.get("name", study_dir.name))
    description = str(study_sec.get("description", ""))

    defaults_sec = cfg.get("defaults", {})
    default_model_key = str(
        defaults_sec.get("model_key") or defaults_sec.get("models") or "default"
    )

    reference = OmegaConf.to_container(cfg.get("reference", {}), resolve=True) if "reference" in cfg else {}
    sweep = OmegaConf.to_container(cfg.get("sweep", {}), resolve=True) if "sweep" in cfg else None

    points = expand_study_points(
        reference=reference,
        sweep=sweep,
        default_model_key=default_model_key,
    )

    return Study(
        name=study_name,
        cfg=cfg,
        points=points,
        study_dir=study_dir,
        description=description,
    )
