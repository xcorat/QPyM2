"""
StudyPoint definition and input readiness routines for QPyM2 studies (v0.6).

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
    model_key: str = "default"
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
                    print(
                        f"[{self.point_id}] Bias toys already completed ({status['n_toys_completed']} rows). Skipping."
                    )
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
        base_norms = np.array(
            [
                float(fit_comps.loc[c, "integral"]) if "integral" in fit_comps.columns else 1.0
                for c in component_names
            ]
        )

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

    def check_sens_toys(
        self,
        study_dir: Union[str, Path],
        n_toys: int = 200,
    ) -> Dict[str, Any]:
        """
        Check completion status of sensitivity toy fits for this study point.
        """
        sdir = Path(study_dir)
        task_csv = sdir / "tasks" / f"toys_sens_{self.point_id}.csv"
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
            "ndbd_cutoff": self.ndbd_cutoff,
            "n_toys_target": n_toys,
            "n_toys_completed": row_count,
            "is_completed": is_done,
            "task_csv": str(task_csv),
            "status": "COMPLETED" if is_done else ("IN_PROGRESS" if row_count > 0 else "PENDING"),
        }

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
        """
        Execute a batch of sensitivity toy fits (background-only hypothesis) for this study point.

        Loads precomputed templates via setup_storage_fit, generates pseudo-experiments
        with zero signal injection (ndbd=0), fits each toy, and saves the task CSV.
        """
        import contextlib
        import numpy as np
        import pandas as pd

        sdir = Path(study_dir)
        tasks_dir = sdir / "tasks"
        tasks_dir.mkdir(parents=True, exist_ok=True)
        task_csv = tasks_dir / f"toys_sens_{self.point_id}.csv"

        if not force:
            status = self.check_sens_toys(study_dir, n_toys=n_toys)
            if status["is_completed"]:
                if verbose:
                    print(
                        f"[{self.point_id}] Sensitivity toys already completed ({status['n_toys_completed']} rows). Skipping."
                    )
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

        base_norms = np.array(
            [
                float(fit_comps.loc[c, "integral"]) if "integral" in fit_comps.columns else 1.0
                for c in component_names
            ]
        )

        # Enforce background-only hypothesis (zero signal)
        norms = base_norms.copy()
        norms[isig] = 0.0

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

        hsum = get_sum_hist(fit_comps, norms)

        rows = []
        for i in range(n_toys):
            # Seed per toy for reproducibility
            toy_seed = int(seed + self.task_id * 10000 + i + 1)
            np.random.seed(toy_seed)

            # Sample pseudo-data under H0
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
                        random_seed=toy_seed,
                    )

            ndbd_post = trace.posterior["ndbd"].values.flatten()
            stats = _get_stats(ndbd_post)
            mode_val = float(stats["mode"])
            upp34_val = float(stats["upp34"])
            uerr = upp34_val - mode_val
            pull = (mode_val - 0.0) / uerr if abs(uerr) > 1e-9 else 0.0

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
                "type_": "sens",
                "ndbd_inj_": 0.0,
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
                print(f"[{self.point_id}] Completed {i + 1}/{n_toys} sensitivity toys", flush=True)

        df_res = pd.DataFrame(rows)
        temp_csv = tasks_dir / f".tmp_{task_csv.name}"
        df_res.to_csv(temp_csv, index=False)
        temp_csv.replace(task_csv)

        if verbose:
            print(f"[{self.point_id}] Saved {len(df_res)} sensitivity toys: {task_csv}")

        return task_csv

    def check_syst_mc(
        self,
        study_dir: Union[str, Path],
        n_samples: int = 100,
    ) -> Dict[str, Any]:
        """
        Check completion status of MC-stats (bmmc) systematic fits for this study point.
        """
        sdir = Path(study_dir)
        task_csv = sdir / "tasks" / f"syst_mc_{self.point_id}.csv"
        task_h5 = sdir / "tasks" / f"syst_mc_{self.point_id}_samples.h5"
        is_done = False
        row_count = 0
        if task_csv.exists() and task_csv.stat().st_size > 0:
            import pandas as pd

            try:
                df = pd.read_csv(task_csv)
                row_count = len(df)
                is_done = row_count >= n_samples and task_h5.exists() and task_h5.stat().st_size > 0
            except Exception:
                pass

        return {
            "point_id": self.point_id,
            "task_id": self.task_id,
            "model_key": self.model_key,
            "co_shift": self.co_shift,
            "ndbd_shift": self.ndbd_shift,
            "co_cut": self.co_cut,
            "ndbd_cutoff": self.ndbd_cutoff,
            "n_samples_target": n_samples,
            "n_samples_completed": row_count,
            "is_completed": is_done,
            "task_csv": str(task_csv),
            "task_h5": str(task_h5),
            "status": "COMPLETED" if is_done else ("IN_PROGRESS" if row_count > 0 else "PENDING"),
        }

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
        """
        Execute a batch of MC-stats (bmmc) systematic fits for this study point.

        Perturbs background source normalisation vectors according to the background
        model covariance matrix (or draws unperturbed raw chain if stat_only=True),
        refits the observed data for each sample, and records summary metrics and
        ndbd posterior samples.
        """
        import contextlib
        import json
        import h5py
        import numpy as np
        import pandas as pd

        sdir = Path(study_dir)
        tasks_dir = sdir / "tasks"
        tasks_dir.mkdir(parents=True, exist_ok=True)
        task_csv = tasks_dir / f"syst_mc_{self.point_id}.csv"
        task_h5 = tasks_dir / f"syst_mc_{self.point_id}_samples.h5"

        if not force:
            status = self.check_syst_mc(study_dir, n_samples=n_samples)
            if status["is_completed"]:
                if verbose:
                    print(
                        f"[{self.point_id}] MC-stats systematic fits already completed ({status['n_samples_completed']} rows). Skipping."
                    )
                return task_csv, task_h5

        # 1. Resolve storage inputs, templates, data, and nominal fit components
        from qpym2.storage.fits import (
            setup_storage_fit,
            assemble_fit_components,
            resolve_storage_coords,
            resolve_fit_template_uuids,
        )
        from qpym2.cfg.config_builders import get_comp_cfg, get_jags_path, get_corrmat_path
        from qpym2.io.bkg_model import get_mcsyst_samples

        coords, hm, param_dir, input_dir = resolve_storage_coords(cfg, model_key=self.model_key)
        jags_fpath = get_jags_path(cfg)
        corrmat_fpath = get_corrmat_path(cfg) if not stat_only else None
        comp_cfg_raw = [dict(c) for c in get_comp_cfg(cfg, model_key=self.model_key)]

        # Run storage fit setup to resolve nominal components, priors, and observed data
        fit_setup = setup_storage_fit(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
            blinded=self.blinded,
            prior_width="5sig",
            mkchain_samples=100,
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            dry_run=False,
        )

        data = fit_setup["data"]
        xedges, yedges = fit_setup["edges"]
        nominal_comps = fit_setup["fit_comps"]
        nominal_priors = list(nominal_comps["prior"])

        # Also get mctable with individual background sources for normalization variation
        tpl_uuids, uuid_base, uuid_co, uuid_ndbd = resolve_fit_template_uuids(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
        )

        mctable, signal_df, _, _ = assemble_fit_components(
            input_dir=input_dir,
            jags_fpath=jags_fpath,
            comp_cfg=comp_cfg_raw,
            hm=hm,
            uuid_base=uuid_base,
            uuid_co=uuid_co,
            uuid_ndbd=uuid_ndbd,
        )

        # 2. Draw n_samples perturbed background normalisation vectors
        if verbose:
            mode_str = "stat-only" if stat_only else "covariance-perturbed"
            print(
                f"[{self.point_id}] Generating {n_samples} {mode_str} normalisation vectors (seed={seed})..."
            )

        mc_samples = get_mcsyst_samples(
            chain=jags_fpath,
            corrmat_fpath=corrmat_fpath,
            n=n_samples,
            seed=seed,
            stat_only=stat_only,
        )

        # 3. Resolve bias calibration if provided
        popt = None
        if isinstance(bias_calibration, dict):
            calib = bias_calibration.get(self.point_id) or bias_calibration
            if "a" in calib and "b" in calib:
                popt = (float(calib["a"]), float(calib["b"]))
        elif bias_calibration is not None:
            cpath = Path(bias_calibration)
            if not cpath.is_absolute():
                cpath = sdir / cpath if (sdir / cpath).exists() else Path(cpath)
            if cpath.exists():
                with open(cpath, "r", encoding="utf-8") as f:
                    calib_dict = json.load(f)
                calib = calib_dict.get(self.point_id) or calib_dict
                if "a" in calib and "b" in calib:
                    popt = (float(calib["a"]), float(calib["b"]))

        # Optional trace saving directory
        traces_dir = sdir / "traces"
        if save_traces:
            traces_dir.mkdir(parents=True, exist_ok=True)

        try:
            import pymc.sampling.parallel as _parallel

            _parallel.threadpool_limits = lambda **kwargs: contextlib.nullcontext()
        except Exception:
            pass

        from qpym2.models import m2
        from qpym2.fit.bias import apply_bias_correction
        from qpym2.fit.results import _get_stats
        from qpym2.utils.mathy import get_mode
        import pymc as pm

        rows = []
        all_ndbd_samples = []
        all_max_logp = []

        for i in range(n_samples):
            # Seed per sample
            sample_seed = int(seed + self.task_id * 10000 + i + 1)
            np.random.seed(sample_seed)

            # Update mctable with sample normalisation
            mt_i = mctable.copy()
            mt_i["norm"] = pd.Series(mc_samples.iloc[i].to_dict()).reindex(mt_i.index)

            # Rebuild grouped background components with perturbed norms
            bkg_comps_i = m2.make_bkgcomps_table(comp_cfg_raw, mt_i, smooth=hm.smooth)
            comps_i = pd.concat([bkg_comps_i, signal_df], join="inner")
            comps_i["prior"] = nominal_priors
            comps_i["varname"] = comps_i.index

            if fit_fn is not None:
                trace = fit_fn(data, comps_i, xedges, yedges)
            else:
                fik_i = m2.make_fik(comps_i, data, xedges, yedges)
                model_i = m2.m2_model(comps_i, fik_i)
                with model_i:
                    trace = pm.sample(
                        draws=draws,
                        tune=tune,
                        chains=chains,
                        cores=1 if single_core else chains,
                        return_inferencedata=True,
                        target_accept=target_accept,
                        progressbar=False,
                        idata_kwargs={"log_likelihood": True},
                        random_seed=sample_seed,
                    )

            if save_traces and hasattr(trace, "to_netcdf"):
                trace_path = traces_dir / f"trace_mc_{self.point_id}_s{i:04d}.nc"
                trace.to_netcdf(str(trace_path))

            ndbd_post = trace.posterior["ndbd"].values.flatten()
            all_ndbd_samples.append(ndbd_post)

            stats = _get_stats(ndbd_post)
            mode_val = float(stats["mode"])
            p16_val = float(stats["p16"])
            p90_val = float(stats["p90"])

            # Log-likelihood mode score
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
                    try:
                        logp_val = float(get_mode(arr))
                    except Exception:
                        logp_val = float(np.median(arr))
                elif len(arr) == 1:
                    logp_val = float(arr[0])
            all_max_logp.append(logp_val)

            mode_bc = p16_bc = p90_bc = np.nan
            if popt is not None:
                corrected = apply_bias_correction(ndbd_post, popt)
                mode_bc = float(get_mode(corrected))
                p16_bc = float(np.percentile(corrected, 16))
                p90_bc = float(np.percentile(corrected, 90))

            row = {
                "sample_idx": i,
                "task_id": self.task_id,
                "point_id": self.point_id,
                "model_key": self.model_key,
                "co_shift": self.co_shift,
                "ndbd_shift": self.ndbd_shift,
                "co_cut": str(self.co_cut),
                "ndbd_cutoff": str(self.ndbd_cutoff),
                "stat_only": stat_only,
                "ndbd_mode": mode_val,
                "ndbd_p16": p16_val,
                "ndbd_p90": p90_val,
                "max_logp": logp_val,
                "ndbd_mode_bc": mode_bc,
                "ndbd_p16_bc": p16_bc,
                "ndbd_p90_bc": p90_bc,
            }

            # Also record mode for background components
            for cname in comps_i.index:
                if cname != "ndbd" and cname in trace.posterior:
                    c_post = trace.posterior[cname].values.flatten()
                    row[f"{cname}_mode"] = float(get_mode(c_post))

            rows.append(row)
            if verbose and ((i + 1) % 10 == 0 or (i + 1) == n_samples):
                print(
                    f"[{self.point_id}] Completed {i + 1}/{n_samples} MC-stats systematic fits",
                    flush=True,
                )

        # Write HDF5 with all posterior samples
        temp_h5 = tasks_dir / f".tmp_{task_h5.name}"
        with h5py.File(temp_h5, "w") as h5f:
            h5f.create_dataset(
                "ndbd_samples",
                data=np.array(all_ndbd_samples, dtype=float),
                compression="gzip",
                compression_opts=4,
            )
            h5f.create_dataset("max_logp", data=np.array(all_max_logp, dtype=float))
            h5f.attrs["point_id"] = self.point_id
            h5f.attrs["n_samples"] = n_samples
            h5f.attrs["stat_only"] = stat_only
        temp_h5.replace(task_h5)

        # Write CSV
        df_res = pd.DataFrame(rows)
        temp_csv = tasks_dir / f".tmp_{task_csv.name}"
        df_res.to_csv(temp_csv, index=False)
        temp_csv.replace(task_csv)

        if verbose:
            print(f"[{self.point_id}] Saved {len(df_res)} MC-stats fits: {task_csv} and {task_h5}")

        return task_csv, task_h5

    def check_syst_lss(
        self,
        study_dir: Union[str, Path],
        category: str = "lss_sscl",
        datasets: Optional[List[int]] = None,
        n_samples: Optional[int] = None,
        dq_mode: str = "pm1",
    ) -> Dict[str, Any]:
        """
        Check completion status of Lineshape Systematic (LSS: sscl or dq) fits for this study point.
        """
        sdir = Path(study_dir)
        cat_norm = "lss_dq" if category in ("dq", "lss_dq") else "lss_sscl"
        task_csv = sdir / "tasks" / f"syst_{cat_norm}_{self.point_id}.csv"
        task_h5 = sdir / "tasks" / f"syst_{cat_norm}_{self.point_id}_samples.h5"

        _ds = datasets if datasets is not None else list(range(3601, 3616))
        if cat_norm == "lss_sscl":
            target_n = len(_ds)
        elif dq_mode == "pm1":
            target_n = 2 * len(_ds)
        else:
            target_n = (n_samples or 100) * len(_ds)

        is_done = False
        row_count = 0
        if task_csv.exists() and task_csv.stat().st_size > 0:
            import pandas as pd

            try:
                df = pd.read_csv(task_csv)
                row_count = len(df)
                is_done = row_count >= target_n
            except Exception:
                pass

        return {
            "point_id": self.point_id,
            "task_id": self.task_id,
            "model_key": self.model_key,
            "co_shift": self.co_shift,
            "ndbd_shift": self.ndbd_shift,
            "co_cut": self.co_cut,
            "ndbd_cutoff": self.ndbd_cutoff,
            "category": cat_norm,
            "n_samples_target": target_n,
            "n_samples_completed": row_count,
            "is_completed": is_done,
            "task_csv": str(task_csv),
            "task_h5": str(task_h5),
            "status": "COMPLETED" if is_done else ("IN_PROGRESS" if row_count > 0 else "PENDING"),
        }

    def run_syst_lss(
        self,
        cfg: DictConfig,
        study_dir: Union[str, Path],
        category: str = "lss_sscl",
        datasets: Optional[List[int]] = None,
        n_samples: Optional[int] = None,
        dq_mode: str = "pm1",
        seed: int = 42,
        chains: int = 2,
        draws: int = 2500,
        tune: int = 500,
        target_accept: float = 0.9,
        single_core: bool = True,
        bias_calibration: Optional[Union[str, Path, Dict[str, Any]]] = None,
        fit_fn: Optional[Any] = None,
        evdata: Optional[Dict[str, Any]] = None,
        lss_params: Optional[Dict[str, Any]] = None,
        lss_data_path: Optional[str] = None,
        save_traces: bool = False,
        force: bool = False,
        verbose: bool = True,
    ) -> Tuple[Path, Path]:
        """
        Execute Lineshape Systematic (LSS: sscl or dq) fits across datasets for this study point.
        """
        import contextlib
        import json
        import h5py
        import numpy as np
        import pandas as pd

        sdir = Path(study_dir)
        tasks_dir = sdir / "tasks"
        tasks_dir.mkdir(parents=True, exist_ok=True)

        cat_norm = "lss_dq" if category in ("dq", "lss_dq") else "lss_sscl"
        task_csv = tasks_dir / f"syst_{cat_norm}_{self.point_id}.csv"
        task_h5 = tasks_dir / f"syst_{cat_norm}_{self.point_id}_samples.h5"

        _ds_list = datasets if datasets is not None else list(range(3601, 3616))

        if not force:
            status = self.check_syst_lss(
                study_dir=study_dir,
                category=cat_norm,
                datasets=_ds_list,
                n_samples=n_samples,
                dq_mode=dq_mode,
            )
            if status["is_completed"]:
                if verbose:
                    print(
                        f"[{self.point_id}] LSS systematic fits ({cat_norm}) already completed ({status['n_samples_completed']} rows). Skipping."
                    )
                return task_csv, task_h5

        # 1. Nominal storage setup (data, components, priors, edges)
        from qpym2.storage.fits import setup_storage_fit, resolve_storage_coords
        from qpym2.cfg.config_builders import get_comp_cfg

        coords, hm, param_dir, input_dir = resolve_storage_coords(cfg, model_key=self.model_key)
        comp_cfg_raw = [dict(c) for c in get_comp_cfg(cfg, model_key=self.model_key)]

        fit_setup = setup_storage_fit(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
            blinded=self.blinded,
            prior_width="5sig",
            chains=chains,
            draws=draws,
            tune=tune,
            target_accept=target_accept,
            dry_run=False,
        )

        data = fit_setup["data"]
        nominal_comps = fit_setup["fit_comps"]
        nominal_priors = list(nominal_comps["prior"])

        # 2. Resolve LSS parameters
        if lss_params is None:
            from qpym2.io.qls_scaling import load_lss_params

            resolved_lss_path = lss_data_path
            if not resolved_lss_path and hasattr(cfg, "syst"):
                resolved_lss_path = getattr(cfg.syst, "lss_data_path", None)
            lss_params = load_lss_params(cfg, n_ds=len(_ds_list), lss_data_path=resolved_lss_path)

        ds_start = _ds_list[0]

        # 3. Resolve staged event data
        if evdata is None:
            from qpym2.cfg.config_builders import get_jags_path, build_input_config
            from qpym2.io.bkg_model import read_jagsh5
            from qpym2.fit.create_evtrees import create_evtrees
            import gc

            jags_fpath = get_jags_path(cfg)
            input_cfg = build_input_config(cfg, model_key=self.model_key)
            mcfit_table, _ = read_jagsh5(jags_fpath)
            comp_fnames = dict(zip(mcfit_table.index, mcfit_table.fname))
            comp_fnames["ndbd"] = "ndbd.root"
            evtrees = create_evtrees(input_cfg, comp_fnames)

            _NP_COLS = ["u0", "v0", "ch1", "ch2", "ds"]
            evdata = {}
            for comp, rdf in evtrees.items():
                if rdf is None:
                    evdata[comp] = None
                else:
                    evdata[comp] = {
                        c: np.array(a) for c, a in rdf.AsNumpy(columns=_NP_COLS).items()
                    }
            del evtrees
            gc.collect()

        # Build mcfit_table for histogramming
        from qpym2.fit.syst_lss import (
            build_sscl_mctable,
            build_dq_mctable,
            run_lss_fit,
            extract_posterior_stats,
        )
        from qpym2.fit.bias import apply_bias_correction
        from qpym2.utils.mathy import get_mode
        from qpym2.fit.params import model_cfg

        mcfg = model_cfg(
            hist_model=hm,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
            co_cuts_dict=dict(getattr(cfg.cuts, "co_cuts", {})) if hasattr(cfg, "cuts") else None,
            ndbd_cuts_dict=dict(getattr(cfg.cuts, "ndbd_cutoff", {}))
            if hasattr(cfg, "cuts")
            else None,
        )
        cuts = mcfg.get("cuts", {})
        cuts_mc = [f for f in [cuts.get("common"), cuts.get("mc_only")] if f]
        cuts_ndbd = [f for f in [cuts.get("common"), cuts.get("ndbd_only")] if f]

        # Background fit table for templates
        from qpym2.storage.fits import assemble_fit_components, resolve_fit_template_uuids

        tpl_uuids, uuid_base, uuid_co, uuid_ndbd = resolve_fit_template_uuids(
            cfg=cfg,
            model_key=self.model_key,
            co_cut=self.co_cut,
            ndbd_cutoff=self.ndbd_cutoff,
            co_shift=self.co_shift,
            ndbd_shift=self.ndbd_shift,
        )
        from qpym2.cfg.config_builders import get_jags_path

        jags_fpath = get_jags_path(cfg)
        mctable, _, _, _ = assemble_fit_components(
            input_dir=input_dir,
            jags_fpath=jags_fpath,
            comp_cfg=comp_cfg_raw,
            hm=hm,
            uuid_base=uuid_base,
            uuid_co=uuid_co,
            uuid_ndbd=uuid_ndbd,
        )

        # 4. Bias calibration
        popt = None
        if isinstance(bias_calibration, dict):
            calib = bias_calibration.get(self.point_id) or bias_calibration
            if "a" in calib and "b" in calib:
                popt = (float(calib["a"]), float(calib["b"]))
        elif bias_calibration is not None:
            cpath = Path(bias_calibration)
            if not cpath.is_absolute():
                cpath = sdir / cpath if (sdir / cpath).exists() else Path(cpath)
            if cpath.exists():
                with open(cpath, "r", encoding="utf-8") as f:
                    calib_dict = json.load(f)
                calib = calib_dict.get(self.point_id) or calib_dict
                if "a" in calib and "b" in calib:
                    popt = (float(calib["a"]), float(calib["b"]))

        # Optional trace saving directory
        traces_dir = sdir / "traces"
        if save_traces:
            traces_dir.mkdir(parents=True, exist_ok=True)

        # Prevent fork deadlocks in PyMC parallel sampler if inside threaded runtime
        try:
            import pymc.sampling.parallel as _parallel

            _parallel.threadpool_limits = lambda **kwargs: contextlib.nullcontext()
        except Exception:
            pass

        # 5. Build fit jobs list: (ds, variant_name, transform_args)
        jobs = []
        if cat_norm == "lss_sscl":
            for ds in _ds_list:
                ids = (
                    ds - ds_start
                    if (ds - ds_start) < len(lss_params["delta_p1_1r"])
                    else _ds_list.index(ds)
                )
                ch_sigma0 = lss_params["sigma0_2d"][ids]
                delta_p1 = lss_params["delta_p1_1r"][ids]
                jobs.append((ds, "sscale_1r", {"ch_sigma0": ch_sigma0, "delta_p1": delta_p1}))
        elif dq_mode == "pm1":
            for ds in _ds_list:
                ids = (
                    ds - ds_start
                    if (ds - ds_start) < len(lss_params["q_pars_arr"])
                    else _ds_list.index(ds)
                )
                q_mean = lss_params["q_pars_arr"][ids]
                q_1r = lss_params["q_pars_sample_1r"][ids]
                q_1l = lss_params["q_pars_sample_1l"][ids]
                jobs.append((ds, "dqres_1r", {"q_pars_mean": q_mean, "q_pars_sample": q_1r}))
                jobs.append((ds, "dqres_1l", {"q_pars_mean": q_mean, "q_pars_sample": q_1l}))
        else:
            # Random parameter draws
            from qpym2.io.qls_scaling import sample_dq_random

            n_dq = n_samples or 100
            dq_random_samples = sample_dq_random(cfg, n_samples=n_dq, seed=seed)
            for si in range(n_dq):
                sname = f"dq_r{si:04d}"
                qps = dq_random_samples[si]
                for ds in _ds_list:
                    ids = (
                        ds - ds_start
                        if (ds - ds_start) < len(lss_params["q_pars_arr"])
                        else _ds_list.index(ds)
                    )
                    q_mean = lss_params["q_pars_arr"][ids]
                    jobs.append((ds, sname, {"q_pars_mean": q_mean, "q_pars_sample": qps[ids]}))

        if verbose:
            print(
                f"[{self.point_id}] Running {len(jobs)} {cat_norm} fits across {len(_ds_list)} datasets..."
            )

        rows = []
        all_ndbd_samples = []

        for idx, (ds, variant, par_dict) in enumerate(jobs):
            fit_seed = int(seed + self.task_id * 10000 + idx + 1)
            rng = np.random.default_rng(fit_seed)

            # Build perturbed mctable
            if cat_norm == "lss_sscl":
                mt_ds, sig_ds = build_sscl_mctable(
                    mc_events_dict=evdata,
                    mcfit_table=mctable,
                    hm=hm,
                    target_ds=ds,
                    ch_sigma0=par_dict["ch_sigma0"],
                    delta_p1=par_dict["delta_p1"],
                    cuts_mc=cuts_mc,
                    cuts_ndbd=cuts_ndbd,
                    rng=rng,
                )
            else:
                mt_ds, sig_ds = build_dq_mctable(
                    mc_events_dict=evdata,
                    mcfit_table=mctable,
                    hm=hm,
                    target_ds=ds,
                    q_pars_mean=par_dict["q_pars_mean"],
                    q_pars_sample=par_dict["q_pars_sample"],
                    cuts_mc=cuts_mc,
                    cuts_ndbd=cuts_ndbd,
                )

            # Fit
            if fit_fn is not None:
                trace = fit_fn(mt_ds, sig_ds, data, hm, comp_cfg_raw, nominal_priors)
            else:
                trace = run_lss_fit(
                    mctable=mt_ds,
                    signal_df=sig_ds,
                    data=data,
                    hm=hm,
                    comp_cfg=comp_cfg_raw,
                    nominal_priors=nominal_priors,
                    chains=chains,
                    draws=draws,
                )

            if save_traces and hasattr(trace, "to_netcdf"):
                trace_path = traces_dir / f"trace_{cat_norm}_{variant}_ds{ds}_{self.point_id}.nc"
                trace.to_netcdf(str(trace_path))

            stats = extract_posterior_stats(trace)
            ndbd_post = stats["ndbd_samples"]
            all_ndbd_samples.append(ndbd_post)

            mode_val = float(stats["component_modes"].get("ndbd", get_mode(ndbd_post)))
            p16_val = float(stats["p16"])
            p90_val = float(stats["p90"])

            mode_bc = p16_bc = p90_bc = np.nan
            if popt is not None:
                corrected = apply_bias_correction(ndbd_post, popt)
                mode_bc = float(get_mode(corrected))
                p16_bc = float(np.percentile(corrected, 16))
                p90_bc = float(np.percentile(corrected, 90))

            row = {
                "ds": ds,
                "variant": variant,
                "category": cat_norm,
                "point_id": self.point_id,
                "task_id": self.task_id,
                "model_key": self.model_key,
                "co_shift": self.co_shift,
                "ndbd_shift": self.ndbd_shift,
                "co_cut": str(self.co_cut),
                "ndbd_cutoff": str(self.ndbd_cutoff),
                "ndbd_mode": mode_val,
                "ndbd_p16": p16_val,
                "ndbd_p90": p90_val,
                "ndbd_mode_bc": mode_bc,
                "ndbd_p16_bc": p16_bc,
                "ndbd_p90_bc": p90_bc,
            }

            for cname, cmode in stats["component_modes"].items():
                if cname != "ndbd":
                    row[f"{cname}_mode"] = cmode

            rows.append(row)
            if verbose and ((idx + 1) % 5 == 0 or (idx + 1) == len(jobs)):
                print(
                    f"[{self.point_id}] Completed {idx + 1}/{len(jobs)} {cat_norm} fits", flush=True
                )

        # Write HDF5 with all posterior samples
        temp_h5 = tasks_dir / f".tmp_{task_h5.name}"
        with h5py.File(temp_h5, "w") as h5f:
            h5f.create_dataset(
                "ndbd_samples",
                data=np.array(all_ndbd_samples, dtype=float),
                compression="gzip",
                compression_opts=4,
            )
            h5f.attrs["point_id"] = self.point_id
            h5f.attrs["category"] = cat_norm
            h5f.attrs["n_fits"] = len(rows)
        temp_h5.replace(task_h5)

        # Write CSV
        df_res = pd.DataFrame(rows)
        temp_csv = tasks_dir / f".tmp_{task_csv.name}"
        df_res.to_csv(temp_csv, index=False)
        temp_csv.replace(task_csv)

        if verbose:
            print(
                f"[{self.point_id}] Saved {len(df_res)} {cat_norm} fits: {task_csv} and {task_h5}"
            )

        return task_csv, task_h5
