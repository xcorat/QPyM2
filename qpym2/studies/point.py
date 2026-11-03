"""
StudyPoint definition and input readiness routines for QPyM2 studies (v0.6).

Represents a single concrete evaluation point in a parameter sweep or study.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union
from omegaconf import DictConfig

from qpym2.storage.fits import build_fit_hash, resolve_fit_template_uuids
from qpym2.storage.inputs import create_storage_inputs_for_point, resolve_storage_coords
from qpym2.storage.paths import get_output_dir


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
        """
        Compute output UUID and target fit paths without performing MCMC sampling.

        Returns:
            Tuple of (fit_params_dict, output_uuid, fit_prefix, output_dir).
        """
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

        fit_prefix = out_dir / f"fit_{output_uuid}"
        return fit_params, output_uuid, fit_prefix, out_dir

    def check_fit(
        self,
        cfg: DictConfig,
        prior_width: str = "5sig",
        mkchain_samples: int = 100,
        chains: int = 8,
        draws: int = 20000,
        tune: int = 1000,
    ) -> Dict[str, Any]:
        """
        Check if fit output files already exist on disk for this point.

        Returns:
            Dictionary with fit existence, output UUID, and output paths.
        """
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
        """
        Execute Bayesian fit for this study point.
        """
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

    def check_bias_toys(
        self,
        study_dir: Union[str, Path],
        n_toys: int = 200,
    ) -> Dict[str, Any]:
        """
        Check completion status of bias toy fits for this study point.
        """
        sdir = Path(study_dir)
        task_csv = sdir / "tasks" / f"toys_bias_{self.point_id}.csv"
        is_done = False
        row_count = 0
        if task_csv.exists() and task_csv.stat().st_size > 0:
            import pandas as pd
            try:
                df = pd.read_csv(task_csv)
                row_count = len(df)
                is_done = row_count >= n_toys
            except Exception:
                pass

        return {
            "point_id": self.point_id,
            "task_id": self.task_id,
            "model_key": self.model_key,
            "co_shift": self.co_shift,
            "ndbd_shift": self.ndbd_shift,
            "co_cut": self.co_cut,
            "n_toys_target": n_toys,
            "n_toys_completed": row_count,
            "is_completed": is_done,
            "task_csv": str(task_csv),
            "status": "COMPLETED" if is_done else ("IN_PROGRESS" if row_count > 0 else "PENDING"),
        }

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
        """
        Execute a batch of bias toy fits for this study point.

        Loads precomputed templates via setup_storage_fit, generates pseudo-experiments
        with uniformly distributed injected signal, fits each toy, and saves the task CSV.
        """
        import contextlib
        import numpy as np
        import pandas as pd

        sdir = Path(study_dir)
        tasks_dir = sdir / "tasks"
        tasks_dir.mkdir(parents=True, exist_ok=True)
        task_csv = tasks_dir / f"toys_bias_{self.point_id}.csv"

        if not force:
            status = self.check_bias_toys(study_dir, n_toys=n_toys)
            if status["is_completed"]:
                if verbose:
                    print(f"[{self.point_id}] Bias toys already completed ({status['n_toys_completed']} rows). Skipping.")
                return task_csv

        # Setup fit objects, templates, and edges from storage
        from qpym2.storage.fits import setup_storage_fit
        fit_setup = setup_storage_fit(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
            prior_width=prior_width,
            mkchain_samples=mkchain_samples,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            dry_run=False,
        )

        fit_comps = fit_setup["fit_comps"]
        xedges, yedges = fit_setup["edges"]
        component_names = list(fit_comps.index)
        if "ndbd" not in component_names:
            raise ValueError(f"Signal component 'ndbd' not found in components: {component_names}")
        isig = component_names.index("ndbd")

        rng = np.random.default_rng(seed + self.task_id * 10000)
        injections = rng.uniform(bias_min, bias_max, size=n_toys)
        base_norms = np.array([
            float(fit_comps.loc[c, "integral"]) if "integral" in fit_comps.columns else 1.0
            for c in component_names
        ])

        # Prevent fork deadlocks in PyMC parallel sampler if inside threaded runtime
        try:
            import pymc.sampling.parallel as _parallel
            _parallel.threadpool_limits = lambda **kwargs: contextlib.nullcontext()
        except Exception:
            pass

        from qpym2.utils.hists import get_sum_hist, sample_hist
        from qpym2.models import m2
        from qpym2.fit.results import _get_stats
        from qpym2.utils.mathy import get_mode
        import pymc as pm

        rows = []
        for i, inj in enumerate(injections):
            norms = base_norms.copy()
            norms[isig] = float(inj)
            hsum = get_sum_hist(fit_comps, norms)
            data_i = sample_hist(hsum, xedges, yedges)

            if fit_fn is not None:
                trace = fit_fn(data_i, fit_comps, xedges, yedges)
            else:
                fik = m2.make_fik(fit_comps, data_i, xedges, yedges)
                model = m2.m2_model(fit_comps, fik)
                with model:
                    trace = pm.sample(
                        draws=draws,
                        tune=tune,
                        chains=chains,
                        cores=1 if single_core else chains,
                        return_inferencedata=True,
                        target_accept=target_accept,
                        progressbar=False,
                        idata_kwargs={"log_likelihood": True},
                        random_seed=int(seed + self.task_id * 10000 + i + 1),
                    )

            ndbd_post = trace.posterior["ndbd"].values.flatten()
            stats = _get_stats(ndbd_post)
            mode_val = float(stats["mode"])
            upp34_val = float(stats["upp34"])
            uerr = upp34_val - mode_val
            pull = (mode_val - float(inj)) / uerr if abs(uerr) > 1e-9 else 0.0

            logp_val = 0.0
            if hasattr(trace, "log_likelihood") and trace.log_likelihood is not None:
                ll = trace.log_likelihood
                if hasattr(ll, "pik"):
                    arr = np.asarray(ll.pik.values).flatten()
                elif isinstance(ll, dict) and "pik" in ll:
                    arr = np.asarray(ll["pik"]).flatten()
                else:
                    arr = np.array([])
                if len(arr) > 1:
                    logp_val = float(get_mode(arr))
                elif len(arr) == 1:
                    logp_val = float(arr[0])

            row = {
                "toy_id": i,
                "task_id": self.task_id,
                "point_id": self.point_id,
                "model_key": self.model_key,
                "co_shift": self.co_shift,
                "ndbd_shift": self.ndbd_shift,
                "co_cut": str(self.co_cut),
                "ndbd_cutoff": str(self.ndbd_cutoff),
                "type_": "bias",
                "ndbd_inj_": float(inj),
                "logp_mode_": logp_val,
                "pull": float(pull),
            }

            for var in component_names:
                var_post = trace.posterior[var].values.flatten()
                v_stats = _get_stats(var_post)
                for sname, sval in v_stats.items():
                    row[f"{var}_{sname}"] = float(sval)

            rows.append(row)
            if verbose and ((i + 1) % 10 == 0 or (i + 1) == n_toys):
                print(f"[{self.point_id}] Completed {i + 1}/{n_toys} bias toys", flush=True)

        df_res = pd.DataFrame(rows)
        temp_csv = tasks_dir / f".tmp_{task_csv.name}"
        df_res.to_csv(temp_csv, index=False)
        temp_csv.replace(task_csv)

        if verbose:
            print(f"[{self.point_id}] Saved {len(df_res)} bias toys: {task_csv}")

        return task_csv

