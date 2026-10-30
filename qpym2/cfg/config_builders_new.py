"""
Configuration builder functions for QPyM2 (v0.6).

Converts Hydra DictConfig objects into dict structures expected by
runner scripts (run_create_staging.py, run_create_input.py, run_fits.py).

These builders provide the bridge between the new Hydra YAML config system
and the existing function signatures that expect specific dict structures.
"""

import math
from typing import Dict, Any, List, Tuple, Optional
from omegaconf import DictConfig, OmegaConf
from pathlib import Path

from qpym2.cfg.hydra_adapter_new import resolve_data_path


def _get_model_key(cfg: DictConfig, model_key: Optional[str]) -> str:
    """Return model_key if given, otherwise return the first key in cfg.models."""
    if model_key is not None:
        return model_key
    if 'models' not in cfg or not cfg.models:
        raise ValueError('No models defined in config and no model_key provided')
    return next(iter(cfg.models))


def _merge_cuts(*cuts: Optional[str]) -> str:
    """Join non-empty cuts with && while skipping trivial true literals."""
    merged: List[str] = []
    for cut in cuts:
        if cut is None:
            continue
        cut_str = str(cut).strip()
        if not cut_str:
            continue
        if cut_str.lower() == 'true':
            continue
        merged.append(cut_str)
    if not merged:
        return 'true'
    return ' && '.join(merged)


def build_staging_config(cfg: DictConfig, model_key: Optional[str] = None) -> Dict[str, Any]:
    """
    Build staging output configuration dict from Hydra config.

    Staging is decoupled from model ROI/selection: the staged files use the
    common transform contract defined in ``cfg.staging`` (versioned, shared).
    Model-specific cuts (co_cut, ndbd_cutoff) belong downstream in input
    construction, not staging.

    Args:
        cfg: Hydra DictConfig object
        model_key: Unused; kept for backward compatibility.

    Returns:
        Dict with keys: outpath, treename, out_cols, cuts
    """
    s_cuts = cfg.staging.cuts
    cuts = _merge_cuts(
        s_cuts.get('base'),
        s_cuts.get('multiplicity'),
        s_cuts.get('data_quality'),
        s_cuts.get('broad_preselection'),
    )

    raw_outpath = cfg.infrastructure.output_paths.staging_dir

    return {
        'outpath': resolve_data_path(cfg, raw_outpath),
        'treename': cfg.qroot.output_tree_name,
        'out_cols': list(cfg.qroot.output_columns),
        'cuts': cuts
    }


def build_mc_config(cfg: DictConfig) -> Dict[str, Any]:
    """
    Build background MC configuration dict from Hydra config.

    Args:
        cfg: Hydra DictConfig object

    Returns:
        Dict with keys: inpath, treename, evar, esum_var, ch1, ch2, ds, t, extra_cuts
    """
    multiplicity_cut = (
        cfg.analysis_cuts.get('multiplicity', 'Multiplicity == 2')
        if hasattr(cfg, 'analysis_cuts') else 'Multiplicity == 2'
    )
    base_cuts = cfg.cuts.get('base_cuts', 'true') if hasattr(cfg, 'cuts') else 'true'
    mc_cuts = cfg.cuts.get('mc_cuts', None) if hasattr(cfg, 'cuts') else None

    raw_inpath = cfg.infrastructure.directories.bkg_mc_dir

    return {
        'inpath': resolve_data_path(cfg, raw_inpath),
        'treename': cfg.qroot.tree_name,
        'evar': cfg.qroot.energy_var,
        'esum_var': cfg.qroot.esum_var,
        'ch1': cfg.qroot.channel_vars[0],
        'ch2': cfg.qroot.channel_vars[1],
        'ds': 'Dataset',  # Standard ROOT variable names
        't': 'Time',
        'extra_cuts': _merge_cuts(multiplicity_cut, base_cuts, mc_cuts)
    }


def build_ndbd_config(cfg: DictConfig, model_key: Optional[str] = None) -> Dict[str, Any]:
    """
    Build NDBD (0νββ signal) configuration dict from Hydra config.

    The NDBD staging job applies only the common staging cuts. The
    model-specific ``ndbd_cutoff`` selection is applied downstream in input
    construction (it feeds the input identity, not the staging identity).

    Args:
        cfg: Hydra DictConfig object
        model_key: Unused; kept for backward compatibility.

    Returns:
        Dict with keys: inpath, fname, out_fname, treename, evar, esum_var,
                       ch1, ch2, ds, t, extra_cuts
    """
    multiplicity_cut = (
        cfg.analysis_cuts.get('multiplicity', 'Multiplicity == 2')
        if hasattr(cfg, 'analysis_cuts') else 'Multiplicity == 2'
    )
    base_cuts = cfg.cuts.get('base_cuts', 'true') if hasattr(cfg, 'cuts') else 'true'
    mc_cuts = cfg.cuts.get('mc_cuts', None) if hasattr(cfg, 'cuts') else None

    ndbd = cfg.staging.sources.ndbd_mc
    raw_inpath = cfg.infrastructure.directories.ndbd_mc_dir

    return {
        'inpath': resolve_data_path(cfg, raw_inpath),
        'fname': ndbd.pattern,
        'out_fname': ndbd.out_fname,
        'treename': cfg.qroot.tree_name,
        'evar': cfg.qroot.energy_var,
        'esum_var': cfg.qroot.esum_var,
        'ch1': cfg.qroot.channel_vars[0],
        'ch2': cfg.qroot.channel_vars[1],
        'ds': 'Dataset',
        't': 'Time',
        'extra_cuts': _merge_cuts(multiplicity_cut, base_cuts, mc_cuts)
    }


def build_data_config(
    cfg: DictConfig, 
    blinded: bool = False,
    model_key: Optional[str] = None
) -> Dict[str, Any]:
    """
    Build experimental data configuration dict from Hydra config.

    Args:
        cfg: Hydra DictConfig object
        blinded: If True, use the blinded output name
        model_key: Model configuration key (unused; reserved for future use)

    Returns:
        Dict with keys: inpath, fname, out_fname, treename, evar, esum_var,
                       ch1, ch2, ds, t, extra_cuts
    """
    multiplicity_cut = (
        cfg.analysis_cuts.get('multiplicity', 'Multiplicity == 2')
        if hasattr(cfg, 'analysis_cuts') else 'Multiplicity == 2'
    )
    base_cuts = cfg.cuts.get('base_cuts', 'true') if hasattr(cfg, 'cuts') else 'true'
    data_cuts = cfg.cuts.get('data_cuts', None) if hasattr(cfg, 'cuts') else None

    s_data = cfg.staging.sources.data
    raw_inpath = cfg.infrastructure.directories.bm_dir
    reduced_fname = cfg.infrastructure.data_files.reduced_data

    out_fname = s_data.out_blinded_fname if blinded else s_data.out_fname

    return {
        'inpath': resolve_data_path(cfg, raw_inpath),
        'fname': reduced_fname,
        'out_fname': out_fname,
        'treename': cfg.qroot.tree_name,
        'evar': cfg.qroot.energy_var,
        'esum_var': cfg.qroot.esum_var,
        'ch1': cfg.qroot.channel_vars[0],
        'ch2': cfg.qroot.channel_vars[1],
        'ds': 'Dataset',
        't': 'Time',
        'extra_cuts': _merge_cuts(multiplicity_cut, base_cuts, data_cuts)
    }


def get_jags_path(cfg: DictConfig) -> str:
    """Extract JAGS HDF5 file path from Hydra config."""
    return resolve_data_path(cfg, cfg.infrastructure.data_files.jags_output)


def get_data_names(cfg: DictConfig) -> Tuple[str, str]:
    """Get data naming conventions from Hydra config."""
    naming = cfg.infrastructure.naming_conventions
    return (naming.data_name, naming.data_blinded_name)


def get_staging_dir(cfg: DictConfig) -> str:
    """Extract staging directory path from Hydra config."""
    return resolve_data_path(cfg, cfg.infrastructure.output_paths.staging_dir)


def get_bkg_mc_dir(cfg: DictConfig) -> str:
    """Extract background MC directory from Hydra config."""
    return resolve_data_path(cfg, cfg.infrastructure.directories.bkg_mc_dir)


def get_ndbd_mc_dir(cfg: DictConfig) -> str:
    """Extract NDBD MC directory from Hydra config."""
    return resolve_data_path(cfg, cfg.infrastructure.directories.ndbd_mc_dir)


def build_input_config(cfg: DictConfig, model_key: Optional[str] = None) -> Dict[str, Any]:
    """Build fit input creation configuration dict from Hydra config."""
    staging_dir = resolve_data_path(cfg, cfg.infrastructure.output_paths.staging_dir)
    model = cfg.models[_get_model_key(cfg, model_key)]

    return {
        'outpath': staging_dir,
        'treename': cfg.qroot.output_tree_name,
        'data_fname': f"{cfg.infrastructure.naming_conventions.data_name}.root",
        'roi': OmegaConf.to_container(model.roi, resolve=True),
        'binning': OmegaConf.to_container(model.binning, resolve=True),
        'data_blinded_name': cfg.infrastructure.naming_conventions.data_blinded_name
    }


def build_write_config(cfg: DictConfig, model_key: Optional[str] = None) -> Dict[str, Any]:
    """Build fit output write configuration dict from Hydra config."""
    fit_dir = resolve_data_path(cfg, cfg.infrastructure.output_paths.fit_dir)
    return {
        'outdir': fit_dir,
        'overwrite': True,
    }


def get_h2uv_model(cfg: DictConfig, model_key: Optional[str] = None):
    """Build a TempConfig histogram model from Hydra config."""
    from qpym2.utils import TempConfig
    from qpym2.utils.hists import HistType

    resolved_key = _get_model_key(cfg, model_key)
    model = cfg.models[resolved_key]
    s2 = math.sqrt(2)

    esum_min = model.roi.esum_min
    esum_max = model.roi.esum_max
    threshold = model.roi.threshold
    if 'u_range' in model.binning or 'v_range' in model.binning:
        raise ValueError(
            f"Model '{resolved_key}' uses deprecated binning.u_range/v_range. "
            "Active h2uv configs derive rotated-axis ranges from roi/esum and threshold. "
            "Please remove u_range/v_range and keep only binning.bin_size (defined in e1/e2 space)."
        )

    binsize = model.binning.bin_size / s2

    umin = esum_min / s2
    umax = esum_max / s2
    vmin = 0
    vmax = (esum_max - 2 * threshold) / s2

    if hasattr(cfg, 'hyperparams') and hasattr(cfg.hyperparams, 'smooth'):
        smooth_enabled = cfg.hyperparams.smooth
    elif hasattr(cfg, 'base_model') and hasattr(cfg.base_model, 'smoothing'):
        smooth_enabled = cfg.base_model.smoothing.enabled
    else:
        smooth_enabled = True

    if smooth_enabled:
        kernel = getattr(cfg.base_model.smoothing, 'kernel', 'k5b') if hasattr(cfg, 'base_model') else 'k5b'
        iterations = getattr(cfg.base_model.smoothing, 'iterations', 1) if hasattr(cfg, 'base_model') else 1
        smooth = (iterations, kernel)
    else:
        smooth = False

    skip_roi = getattr(cfg.hyperparams, 'skip_roi', None) if hasattr(cfg, 'hyperparams') else None
    name = getattr(model, 'name', f'hmu_{resolved_key}')

    return TempConfig(
        esum_min=esum_min,
        esum_max=esum_max,
        threshold=threshold,
        bin_size=float(model.binning.bin_size),
        binsize=binsize,
        umin=umin,
        umax=umax,
        vmin=vmin,
        vmax=vmax,
        smooth=smooth,
        skip_roi=skip_roi,
        hist_type=HistType.h2uv,
        name=name,
    )


def get_comp_cfg(cfg: DictConfig, model_key: Optional[str] = None) -> List[Dict[str, str]]:
    """Read background component grouping configuration from Hydra config."""
    if model_key is not None and 'models' in cfg and model_key in cfg.models:
        model = cfg.models[model_key]
        if hasattr(model, 'components') and model.components:
            return [OmegaConf.to_container(c, resolve=True) for c in model.components]
    elif 'models' in cfg and cfg.models:
        try:
            mkey = _get_model_key(cfg, None)
            first_model = cfg.models[mkey]
            if hasattr(first_model, 'components') and first_model.components:
                return [OmegaConf.to_container(c, resolve=True) for c in first_model.components]
        except Exception:
            pass

    if hasattr(cfg, 'components') and cfg.components:
        return [OmegaConf.to_container(c, resolve=True) for c in cfg.components]

    from qpym2.fit.fit import comp_cfg
    return [dict(c) for c in comp_cfg]


def get_input_fname(cfg: DictConfig, fitdir: str, skip_roi: bool) -> str:
    """Build the fit input HDF5 path from config naming conventions."""
    naming = cfg.infrastructure.naming_conventions
    base = naming.get('input_noroi_name', 'input_noroi') if skip_roi else naming.get('input_name', 'input')
    return f"{fitdir}/{base}.h5"


def get_corrmat_path(cfg: DictConfig) -> Optional[str]:
    """Extract background systematic correlation matrix path from Hydra config."""
    return resolve_data_path(cfg, getattr(cfg.infrastructure.data_files, 'bkg_syst_corrmat', None))


def get_fit_dir(cfg: DictConfig) -> str:
    """Extract fit directory from Hydra config."""
    return resolve_data_path(cfg, cfg.infrastructure.output_paths.fit_dir)


def get_data_dir(cfg: DictConfig) -> Optional[str]:
    """Return the base results data_dir, or None if not set."""
    val = getattr(cfg.infrastructure.output_paths, 'data_dir', None)
    return resolve_data_path(cfg, val) if val else None


def get_plots_dir(cfg: DictConfig) -> str:
    """Return plots directory."""
    data_dir = get_data_dir(cfg)
    if data_dir:
        return str(Path(data_dir) / 'plots')
    plots = getattr(cfg.infrastructure.output_paths, 'plots_dir', None)
    if plots:
        return resolve_data_path(cfg, plots)
    return str(Path(get_fit_dir(cfg)) / 'plots')


def get_syst_dir(cfg: DictConfig) -> str:
    """Return systematic results directory."""
    data_dir = get_data_dir(cfg)
    if data_dir:
        return str(Path(data_dir) / 'syst_results')
    syst = getattr(cfg.infrastructure.output_paths, 'systematic_dir', None)
    if syst:
        return resolve_data_path(cfg, syst)
    raise ValueError(
        "No systematic results directory configured. "
        "Set data_dir or systematic_dir in output_paths."
    )


def get_toy_results_dir(cfg: DictConfig) -> str:
    """Return toy results directory."""
    data_dir = get_data_dir(cfg)
    if data_dir:
        return str(Path(data_dir) / 'toy_results')
    toy = getattr(cfg.infrastructure.output_paths, 'toy_results_dir', None)
    if toy:
        return resolve_data_path(cfg, toy)
    raise ValueError(
        "No toy results directory configured. "
        "Set data_dir or toy_results_dir in output_paths."
    )


def get_artifact_store_dir(cfg: DictConfig) -> str:
    """Extract artifact store directory path from Hydra config."""
    raw_path = getattr(cfg.infrastructure.output_paths, 'artifact_store', 'artifact_store')
    return resolve_data_path(cfg, raw_path)
