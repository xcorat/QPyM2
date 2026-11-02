"""
StudyPoint definition and input readiness routines for QPyM2 studies (v0.6).

Represents a single concrete evaluation point in a parameter sweep or study.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Optional
from omegaconf import DictConfig

from qpym2.storage.inputs import create_storage_inputs_for_point


@dataclass
class StudyPoint:
    """
    Representation of an individual parameter point within a study.

    Attributes:
        point_id: Unique string label for this point (e.g. 'shift_0.8' or 'p0').
        task_id: Sequential integer index.
        model_key: Named model key (e.g. 'tech_note_200').
        co_cut: Co-60 cut key or expression (e.g. 'gamma10_v150').
        ndbd_cutoff: NDBD cutoff key or expression (e.g. 'u1740').
        co_shift: Co-60 energy shift in keV.
        ndbd_shift: NDBD energy shift in keV.
        blinded: Whether data is blinded.
        params: Optional dictionary of arbitrary extra parameters.
    """
    point_id: str
    task_id: int = 0
    model_key: str = "tech_note_200"
    co_cut: Optional[str] = None
    ndbd_cutoff: Optional[str] = None
    co_shift: float = 0.0
    ndbd_shift: float = 0.0
    blinded: bool = False
    params: Dict[str, Any] = field(default_factory=dict)

    def get_required_inputs(self, cfg: DictConfig) -> Dict[str, Any]:
        """
        Inspect and resolve all required templates and data without modifying disk.

        Args:
            cfg: Merged Hydra configuration.

        Returns:
            Point manifest dictionary with expected paths, UUIDs, and cache flags.
        """
        return create_storage_inputs_for_point(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
            blinded=self.blinded,
            dry_run=True,
        )

    def check_inputs(self, cfg: DictConfig) -> Dict[str, Any]:
        """
        Check existence of all required inputs on disk.

        Args:
            cfg: Merged Hydra configuration.

        Returns:
            Dictionary with readiness metrics:
            - point_id, task_id, model_key, co_shift, etc.
            - total_templates, cached_templates, missing_templates
            - missing_template_names, missing_template_paths
            - data_path, data_cached
            - is_ready (True if all templates and data exist)
        """
        manifest = self.get_required_inputs(cfg)
        templates = manifest.get("templates", {})
        data_info = manifest.get("data", {})

        total_templates = len(templates)
        cached_templates = sum(1 for t in templates.values() if t.get("cached", False))
        missing_templates = total_templates - cached_templates

        missing_names = [k for k, t in templates.items() if not t.get("cached", False)]
        missing_paths = [t.get("path") for t in templates.values() if not t.get("cached", False)]

        data_cached = bool(data_info.get("cached", False))
        data_path = str(data_info.get("path", ""))

        is_ready = (missing_templates == 0) and data_cached

        return {
            "point_id": self.point_id,
            "task_id": self.task_id,
            "model_key": self.model_key,
            "co_cut": self.co_cut,
            "ndbd_cutoff": self.ndbd_cutoff,
            "co_shift": self.co_shift,
            "ndbd_shift": self.ndbd_shift,
            "blinded": self.blinded,
            "total_templates": total_templates,
            "cached_templates": cached_templates,
            "missing_templates": missing_templates,
            "missing_template_names": missing_names,
            "missing_template_paths": missing_paths,
            "data_path": data_path,
            "data_cached": data_cached,
            "is_ready": is_ready,
            "status": "READY" if is_ready else "MISSING_INPUTS",
            "input_dir": manifest.get("input_dir", ""),
            "param_dir": manifest.get("param_dir", ""),
        }

    def generate_inputs(
        self,
        cfg: DictConfig,
        overwrite: bool = False,
        hist_cache: Optional[Dict[Any, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Materialize missing input templates and data on disk.

        Args:
            cfg: Merged Hydra configuration.
            overwrite: If True, regenerate and overwrite existing files.
            hist_cache: Optional shared in-memory histogram cache.

        Returns:
            Manifest dictionary returned by create_storage_inputs_for_point.
        """
        return create_storage_inputs_for_point(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
            blinded=self.blinded,
            overwrite=overwrite,
            dry_run=False,
            hist_cache=hist_cache,
        )
