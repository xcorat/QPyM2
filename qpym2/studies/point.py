"""StudyPoint definition and input readiness routines for QPyM2 studies (v0.6).

Represents a single concrete evaluation point in a parameter sweep or study.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from omegaconf import DictConfig

from qpym2.storage.fits import build_fit_hash, resolve_fit_template_uuids
from qpym2.storage.inputs import create_storage_inputs_for_point, resolve_storage_coords
from qpym2.storage.paths import get_output_dir


@dataclass
class StudyPoint:
    """Representation of an individual parameter point within a study.

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
    model_key: str = "default"
    co_cut: Optional[str] = None
    ndbd_cutoff: Optional[str] = None
    co_shift: float = 0.0
    ndbd_shift: float = 0.0
    blinded: bool = False
    params: Dict[str, Any] = field(default_factory=dict)

    def get_required_inputs(self, cfg: DictConfig) -> Dict[str, Any]:
        """Inspect and resolve all required templates and data without modifying disk."""
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
        """Check existence of all required inputs on disk."""
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
        """Materialize missing input templates and data on disk."""
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

    def get_fit_hash(
        self,
        cfg: DictConfig,
        prior_width: str = "5sig",
        mkchain_samples: int = 100,
        chains: int = 8,
        draws: int = 20000,
        tune: int = 1000,
        target_accept: float = 0.9,
        study_context: Optional[Dict[str, Any]] = None,
    ) -> Tuple[Dict[str, Any], str, Path, Path]:
        """Compute output UUID and target fit paths without performing MCMC sampling."""
        coords, hm, param_dir, _ = resolve_storage_coords(cfg, model_key=self.model_key)
        out_dir = get_output_dir(param_dir)

        tpl_uuids, _, _, _ = resolve_fit_template_uuids(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
        )

        smoothing_cfg = {
            "enabled": bool(getattr(hm, "smooth", None)),
            "method": "smooth_nph2",
            "smooth": list(hm.smooth) if getattr(hm, "smooth", None) else None,
        }

        if prior_width in ("syst", "2syst"):
            sigma_scale = 1.0 if prior_width == "syst" else 2.0
            use_sample_std = True
        else:
            sigma_scale = 5.0 if prior_width == "5sig" else 10.0
            use_sample_std = False

        prior_cfg = {
            "prior_width": prior_width,
            "sigma_scale": sigma_scale,
            "use_sample_std": use_sample_std,
            "mkchain_samples": mkchain_samples,
        }

        mcmc_cfg = {
            "chains": chains,
            "draws": draws,
            "tune": tune,
            "target_accept": target_accept,
        }
        data_source_name = "data_blinded" if self.blinded else "data"

        skip_roi = self.params.get("skip_roi")
        roi_cfg = {"skip_roi": list(skip_roi) if skip_roi else None}

        fit_params, output_uuid = build_fit_hash(
            model_name="m2",
            data_source=data_source_name,
            template_uuids=tpl_uuids,
            smoothing=smoothing_cfg,
            priors=prior_cfg,
            roi_settings=roi_cfg,
            mcmc=mcmc_cfg,
            study_context=study_context,
        )

        suffix = self.params.get("suffix")
        metadata = dict(fit_params)
        metadata["mcmc"] = mcmc_cfg
        metadata["output_uuid"] = output_uuid
        if suffix:
            metadata["suffix"] = suffix
            fit_prefix = out_dir / f"fit_{output_uuid}_{suffix}"
        else:
            fit_prefix = out_dir / f"fit_{output_uuid}"

        return metadata, output_uuid, fit_prefix, out_dir

    def check_fit(
        self,
        cfg: DictConfig,
        prior_width: str = "5sig",
        mkchain_samples: int = 100,
        chains: int = 8,
        draws: int = 20000,
        tune: int = 1000,
    ) -> Dict[str, Any]:
        """Check if fit output files already exist on disk for this point."""
        _, output_uuid, fit_prefix, out_dir = self.get_fit_hash(
            cfg=cfg,
            prior_width=prior_width,
            mkchain_samples=mkchain_samples,
            chains=chains,
            draws=draws,
            tune=tune,
        )
        h5_path = fit_prefix.with_suffix(".h5")
        nc_path = fit_prefix.with_suffix(".nc")
        is_fitted = h5_path.exists()

        return {
            "point_id": self.point_id,
            "task_id": self.task_id,
            "model_key": self.model_key,
            "co_shift": self.co_shift,
            "ndbd_shift": self.ndbd_shift,
            "co_cut": self.co_cut,
            "output_uuid": output_uuid,
            "is_fitted": is_fitted,
            "has_trace": nc_path.exists(),
            "fit_prefix": str(fit_prefix),
            "h5_path": str(h5_path),
            "nc_path": str(nc_path),
            "output_dir": str(out_dir),
            "status": "FITTED" if is_fitted else "PENDING_FIT",
        }

    def run_fit(
        self,
        cfg: DictConfig,
        save_trace: bool = False,
        force: bool = False,
        chains: int = 8,
        draws: int = 20000,
        tune: int = 1000,
        prior_width: str = "5sig",
        mkchain_samples: int = 100,
        single_core: bool = False,
        study_dir: Optional[Union[str, Path]] = None,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """Execute Bayesian fit for this study point."""
        from qpym2.run.run_single_storage_fit import run_single_storage_fit

        return run_single_storage_fit(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
            blinded=self.blinded,
            prior_width=prior_width,
            mkchain_samples=mkchain_samples,
            chains=chains,
            draws=draws,
            tune=tune,
            single_core=single_core,
            save_trace=save_trace,
            study_dir=study_dir,
            force=force,
            dry_run=False,
            verbose=verbose,
        )

    # --- Modular Task Delegates ---

    def check_bias_toys(
        self,
        study_dir: Union[str, Path],
        n_toys: int = 200,
    ) -> Dict[str, Any]:
        """Check completion status of bias toy fits for this study point."""
        from qpym2.studies.toys import check_bias_toys

        return check_bias_toys(self, study_dir=study_dir, n_toys=n_toys)

    def run_bias_toys(
        self,
        cfg: DictConfig,
        study_dir: Union[str, Path],
        n_toys: int = 200,
        bias_min: float = 5.0,
        bias_max: float = 200.0,
        seed: int = 42,
        chains: int = 4,
        draws: int = 2000,
        tune: int = 500,
        target_accept: float = 0.9,
        single_core: bool = True,
        prior_width: str = "5sig",
        mkchain_samples: int = 100,
        fit_fn: Optional[Any] = None,
        force: bool = False,
        verbose: bool = True,
    ) -> Path:
        """Execute a batch of bias toy fits for this study point."""
        from qpym2.studies.toys import run_bias_toys

        return run_bias_toys(
            point=self,
            cfg=cfg,
            study_dir=study_dir,
            n_toys=n_toys,
            bias_min=bias_min,
            bias_max=bias_max,
            seed=seed,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            single_core=single_core,
            prior_width=prior_width,
            mkchain_samples=mkchain_samples,
            fit_fn=fit_fn,
            force=force,
            verbose=verbose,
        )

    def check_sens_toys(
        self,
        study_dir: Union[str, Path],
        n_toys: int = 200,
    ) -> Dict[str, Any]:
        """Check completion status of sensitivity toy fits for this study point."""
        from qpym2.studies.toys import check_sens_toys

        return check_sens_toys(self, study_dir=study_dir, n_toys=n_toys)

    def run_sens_toys(
        self,
        cfg: DictConfig,
        study_dir: Union[str, Path],
        n_toys: int = 200,
        seed: int = 42,
        chains: int = 4,
        draws: int = 2000,
        tune: int = 500,
        target_accept: float = 0.9,
        single_core: bool = True,
        prior_width: str = "5sig",
        mkchain_samples: int = 100,
        fit_fn: Optional[Any] = None,
        force: bool = False,
        verbose: bool = True,
    ) -> Path:
        """Execute a batch of sensitivity toy fits (background-only) for this study point."""
        from qpym2.studies.toys import run_sens_toys

        return run_sens_toys(
            point=self,
            cfg=cfg,
            study_dir=study_dir,
            n_toys=n_toys,
            seed=seed,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            single_core=single_core,
            prior_width=prior_width,
            mkchain_samples=mkchain_samples,
            fit_fn=fit_fn,
            force=force,
            verbose=verbose,
        )

    def check_syst_mc(
        self,
        study_dir: Union[str, Path],
        n_samples: int = 100,
    ) -> Dict[str, Any]:
        """Check completion status of MC-stats (bmmc) systematic fits for this study point."""
        from qpym2.studies.syst_mc import check_syst_mc

        return check_syst_mc(self, study_dir=study_dir, n_samples=n_samples)

    def run_syst_mc(
        self,
        cfg: DictConfig,
        study_dir: Union[str, Path],
        n_samples: int = 100,
        seed: int = 42,
        chains: int = 2,
        draws: int = 2500,
        tune: int = 500,
        target_accept: float = 0.9,
        stat_only: bool = False,
        single_core: bool = True,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        fit_fn: Optional[Any] = None,
        save_traces: bool = False,
        force: bool = False,
        verbose: bool = True,
    ) -> Tuple[Path, Path]:
        """Execute a batch of MC-stats (bmmc) systematic fits for this study point."""
        from qpym2.studies.syst_mc import run_syst_mc

        return run_syst_mc(
            point=self,
            cfg=cfg,
            study_dir=study_dir,
            n_samples=n_samples,
            seed=seed,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            stat_only=stat_only,
            single_core=single_core,
            bias_calibration=bias_calibration,
            fit_fn=fit_fn,
            save_traces=save_traces,
            force=force,
            verbose=verbose,
        )

    def check_syst_lss(
        self,
        study_dir: Union[str, Path],
        category: str = "lss_sscl",
        datasets: Optional[List[int]] = None,
        n_samples: Optional[int] = None,
        dq_mode: str = "pm1",
    ) -> Dict[str, Any]:
        """Check completion status of Lineshape Systematic (LSS: sscl or dq) fits for this point."""
        from qpym2.studies.syst_lss import check_syst_lss

        return check_syst_lss(
            point=self,
            study_dir=study_dir,
            category=category,
            datasets=datasets,
            n_samples=n_samples,
            dq_mode=dq_mode,
        )

    def run_syst_lss(
        self,
        cfg: DictConfig,
        study_dir: Union[str, Path],
        category: str = "lss_sscl",
        datasets: Optional[List[int]] = None,
        dq_mode: str = "pm1",
        n_samples: Optional[int] = None,
        lss_data_path: Optional[Union[str, Path]] = None,
        lss_params: Optional[Dict[str, Any]] = None,
        evdata: Optional[Dict[str, Any]] = None,
        seed: int = 42,
        chains: int = 2,
        draws: int = 2500,
        tune: int = 500,
        target_accept: float = 0.9,
        single_core: bool = True,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        fit_fn: Optional[Any] = None,
        save_traces: bool = False,
        force: bool = False,
        verbose: bool = True,
    ) -> Tuple[Path, Path]:
        """Execute Lineshape Systematic (LSS: SSCL or dQ) fits for this study point."""
        from qpym2.studies.syst_lss import run_syst_lss

        return run_syst_lss(
            point=self,
            cfg=cfg,
            study_dir=study_dir,
            category=category,
            datasets=datasets,
            dq_mode=dq_mode,
            n_samples=n_samples,
            lss_data_path=lss_data_path,
            lss_params=lss_params,
            evdata=evdata,
            seed=seed,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            single_core=single_core,
            bias_calibration=bias_calibration,
            fit_fn=fit_fn,
            save_traces=save_traces,
            force=force,
            verbose=verbose,
        )

    def check_syst_bm_threshold(
        self,
        study_dir: Union[str, Path],
        variant: str = "Threshold250",
    ) -> Dict[str, Any]:
        """Check completion status of BM threshold systematic fit for this study point."""
        from qpym2.studies.syst_bm import check_syst_bm_threshold

        return check_syst_bm_threshold(self, study_dir=study_dir, variant=variant)

    def run_syst_bm_threshold(
        self,
        cfg: DictConfig,
        study_dir: Union[str, Path],
        variant: str = "Threshold250",
        model_key_syst: str = "tech_note_250",
        results_raw_path: Optional[Union[str, Path]] = None,
        seed: int = 42,
        chains: int = 2,
        draws: int = 2500,
        tune: int = 500,
        target_accept: float = 0.9,
        single_core: bool = True,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        dedicated_bias_csv: Optional[Union[str, Path]] = None,
        fit_fn: Optional[Any] = None,
        save_traces: bool = False,
        force: bool = False,
        verbose: bool = True,
    ) -> Tuple[Path, Path]:
        """Execute BM threshold systematic fit (e.g. Threshold250 vs Threshold200)."""
        from qpym2.studies.syst_bm import run_syst_bm_threshold

        return run_syst_bm_threshold(
            point=self,
            cfg=cfg,
            study_dir=study_dir,
            variant=variant,
            model_key_syst=model_key_syst,
            results_raw_path=results_raw_path,
            seed=seed,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            single_core=single_core,
            bias_calibration=bias_calibration,
            dedicated_bias_csv=dedicated_bias_csv,
            fit_fn=fit_fn,
            save_traces=save_traces,
            force=force,
            verbose=verbose,
        )

    def check_syst_bm_binning(
        self,
        study_dir: Union[str, Path],
        variants: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """Check completion status of BM binning systematic fits for this study point."""
        from qpym2.studies.syst_bm import check_syst_bm_binning

        return check_syst_bm_binning(self, study_dir=study_dir, variants=variants)

    def run_syst_bm_binning(
        self,
        cfg: DictConfig,
        study_dir: Union[str, Path],
        variants: Optional[List[str]] = None,
        results_raw_path: Optional[Union[str, Path]] = None,
        seed: int = 42,
        chains: int = 2,
        draws: int = 2500,
        tune: int = 500,
        target_accept: float = 0.9,
        single_core: bool = True,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        fit_fn: Optional[Any] = None,
        save_traces: bool = False,
        force: bool = False,
        verbose: bool = True,
    ) -> Tuple[Path, Path]:
        """Execute BM binning systematic fits across variants (MinBinSize, IntraAlpha)."""
        from qpym2.studies.syst_bm import run_syst_bm_binning

        return run_syst_bm_binning(
            point=self,
            cfg=cfg,
            study_dir=study_dir,
            variants=variants,
            results_raw_path=results_raw_path,
            seed=seed,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            single_core=single_core,
            bias_calibration=bias_calibration,
            fit_fn=fit_fn,
            save_traces=save_traces,
            force=force,
            verbose=verbose,
        )

    def check_syst_template_binning(
        self,
        study_dir: Union[str, Path],
        bin_sizes: Optional[List[float]] = None,
    ) -> Dict[str, Any]:
        """Check completion status of template binning (bin_size) systematic fits."""
        from qpym2.studies.syst_template import check_syst_template_binning

        return check_syst_template_binning(self, study_dir=study_dir, bin_sizes=bin_sizes)

    def run_syst_template_binning(
        self,
        cfg: DictConfig,
        study_dir: Union[str, Path],
        bin_sizes: Optional[List[float]] = None,
        seed: int = 42,
        chains: int = 2,
        draws: int = 2500,
        tune: int = 500,
        target_accept: float = 0.9,
        single_core: bool = True,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        fit_fn: Optional[Any] = None,
        save_traces: bool = False,
        force: bool = False,
        verbose: bool = True,
    ) -> Tuple[Path, Path]:
        """Execute template likelihood binning fits for coarser bin sizes."""
        from qpym2.studies.syst_template import run_syst_template_binning

        return run_syst_template_binning(
            point=self,
            cfg=cfg,
            study_dir=study_dir,
            bin_sizes=bin_sizes,
            seed=seed,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            single_core=single_core,
            bias_calibration=bias_calibration,
            fit_fn=fit_fn,
            save_traces=save_traces,
            force=force,
            verbose=verbose,
        )
