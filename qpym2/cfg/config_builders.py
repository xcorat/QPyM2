"""
Configuration builder functions for QPyM2 (v0.6).

Converts Hydra DictConfig objects into dict structures expected by
runners and analysis operations.

These builders provide the bridge between the Hydra YAML config system
and function signatures that expect specific dict structures.
"""

import math
from typing import Dict, Any, List, Tuple, Optional
from omegaconf import DictConfig, OmegaConf
from pathlib import Path

from qpym2.cfg.hydra_adapter import resolve_data_path


def _get_model_key(cfg: DictConfig, model_key: Optional[str] = "default") -> str:
    """Resolve model key, supporting 'default' via cfg.defaults.model_key.

    Args:
        cfg: Composed Hydra DictConfig object.
        model_key: Named model identifier, or 'default' / None.

    Returns:
        Canonical model key present in cfg.models.

    Raises:
        ValueError: If model_key is unresolvable or no models exist.
        KeyError: If model_key is not found in cfg.models.
    """
    if 'models' not in cfg or not cfg.models:
        raise ValueError("No models defined in configuration.")

    if model_key in (None, "default"):
        defaults_sec = getattr(cfg, 'defaults', None)
        default_from_cfg = getattr(defaults_sec, 'model_key', None) if defaults_sec else None
        if default_from_cfg and default_from_cfg in cfg.models:
            return str(default_from_cfg)

        default_model = getattr(cfg, 'default_model', None)
        if not default_model and hasattr(cfg, 'models'):
            default_model = getattr(cfg.models, 'default_model', None)
        if default_model and default_model in cfg.models:
            return str(default_model)

        if "default" in cfg.models:
            target = cfg.models["default"]
            if isinstance(target, str) and target in cfg.models:
                return target
            return "default"

        raise ValueError(
            "model_key 'default' cannot be resolved: no 'default' model, 'default_model', or 'defaults.model_key' "
            f"defined in configuration. Available models: {list(cfg.models.keys())}"
        )

    if model_key not in cfg.models:
        raise KeyError(
            f"model_key '{model_key}' not found in cfg.models. "
            f"Available models: {list(cfg.models.keys())}"
        )
    return str(model_key)


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
        model_key: Kept for signature compatibility; unused for staging.

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

    qroot = getattr(cfg, 'qroot', None) or getattr(cfg.staging, 'qroot', None)
    if qroot is None:
        raise ValueError("No qroot definition found in cfg or cfg.staging.")

    raw_outpath = cfg.infrastructure.output_paths.staging_dir

    return {
        'outpath': resolve_data_path(cfg, raw_outpath),
        'treename': qroot.output_tree_name,
        'out_cols': list(qroot.output_columns),
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

    input_root_key = getattr(getattr(cfg, 'staging', None), 'staging_input_root', 'cuore_data_root')
    raw_inpath = cfg.infrastructure.directories.bkg_mc_dir
    qroot = getattr(cfg, 'qroot', None) or getattr(cfg.staging, 'qroot', None)

    return {
        'inpath': resolve_data_path(cfg, raw_inpath, root_key=input_root_key),
        'treename': qroot.tree_name,
        'evar': qroot.energy_var,
        'esum_var': qroot.esum_var,
        'ch1': qroot.channel_vars[0],
        'ch2': qroot.channel_vars[1],
        'ds': 'Dataset',
        't': 'Time',
        'extra_cuts': _merge_cuts(multiplicity_cut, base_cuts, mc_cuts)
    }


def build_ndbd_config(cfg: DictConfig, model_key: Optional[str] = None) -> Dict[str, Any]:
    """
    Build NDBD (0νββ signal) configuration dict from Hydra config.

    The NDBD staging job applies only the common staging cuts. The
    model-specific ``ndbd_cutoff`` selection is applied downstream in input
    construction.

    Args:
        cfg: Hydra DictConfig object
        model_key: Kept for signature compatibility; unused for staging.

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
    input_root_key = getattr(getattr(cfg, 'staging', None), 'staging_input_root', 'cuore_data_root')
    raw_inpath = cfg.infrastructure.directories.ndbd_mc_dir
    qroot = getattr(cfg, 'qroot', None) or getattr(cfg.staging, 'qroot', None)

    return {
        'inpath': resolve_data_path(cfg, raw_inpath, root_key=input_root_key),
        'fname': ndbd.pattern,
        'out_fname': ndbd.out_fname,
        'treename': qroot.tree_name,
        'evar': qroot.energy_var,
        'esum_var': qroot.esum_var,
        'ch1': qroot.channel_vars[0],
        'ch2': qroot.channel_vars[1],
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
        model_key: Kept for signature compatibility; unused for staging.

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
    input_root_key = getattr(getattr(cfg, 'staging', None), 'staging_input_root', 'cuore_data_root')
    raw_inpath = cfg.infrastructure.directories.bm_dir
    reduced_fname = cfg.infrastructure.data_files.reduced_data
    qroot = getattr(cfg, 'qroot', None) or getattr(cfg.staging, 'qroot', None)

    out_fname = s_data.out_blinded_fname if blinded else s_data.out_fname

    return {
        'inpath': resolve_data_path(cfg, raw_inpath, root_key=input_root_key),
        'fname': reduced_fname,
        'out_fname': out_fname,
        'treename': qroot.tree_name,
        'evar': qroot.energy_var,
        'esum_var': qroot.esum_var,
        'ch1': qroot.channel_vars[0],
        'ch2': qroot.channel_vars[1],
        'ds': 'Dataset',
        't': 'Time',
        'extra_cuts': _merge_cuts(multiplicity_cut, base_cuts, data_cuts)
    }


def get_jags_path(cfg: DictConfig) -> str:
    """Extract JAGS HDF5 file path from Hydra config."""
    df = getattr(cfg.infrastructure, 'data_files', None)
    raw = getattr(df, 'jags_output', None) or getattr(df, 'bkg_model_jags', None)
    return resolve_data_path(cfg, raw)


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


def build_input_config(cfg: DictConfig, model_key: Optional[str] = "default") -> Dict[str, Any]:
    """Build fit input creation configuration dict from Hydra config."""
    staging_dir = resolve_data_path(cfg, cfg.infrastructure.output_paths.staging_dir)
    resolved_key = _get_model_key(cfg, model_key)
    model = cfg.models[resolved_key]

    qroot = getattr(cfg, 'qroot', None) or getattr(cfg.staging, 'qroot', None)
    treename = qroot.output_tree_name if qroot is not None else "uvtree"

    return {
        'outpath': staging_dir,
        'treename': treename,
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


def get_h2uv_model(
    cfg: DictConfig,
    model_key: Optional[str] = "default",
    smooth: Optional[Any] = None,
    skip_roi: Optional[Any] = None,
):
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

    # Smoothing resolution: explicit argument > model definition > default (1, 'k5b')
    if smooth is not None:
        effective_smooth = smooth
    elif hasattr(model, 'smoothing'):
        m_smoothing = model.smoothing
        if isinstance(m_smoothing, bool):
            effective_smooth = (1, 'k5b') if m_smoothing else False
        else:
            enabled = getattr(m_smoothing, 'enabled', True)
            if enabled:
                kernel = getattr(m_smoothing, 'kernel', 'k5b')
                iterations = getattr(m_smoothing, 'iterations', 1)
                effective_smooth = (iterations, kernel)
            else:
                effective_smooth = False
    elif hasattr(model, 'smooth'):
        effective_smooth = model.smooth
    else:
        effective_smooth = (1, 'k5b')

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
        smooth=effective_smooth,
        skip_roi=skip_roi,
        hist_type=HistType.h2uv,
        name=name,
    )


def get_comp_cfg(cfg: DictConfig, model_key: Optional[str] = "default") -> List[Dict[str, str]]:
    """Read background component grouping configuration from Hydra config."""
    resolved_key = _get_model_key(cfg, model_key)
    model = cfg.models[resolved_key]

    if hasattr(model, 'components') and model.components:
        return [OmegaConf.to_container(c, resolve=True) for c in model.components]

    if hasattr(cfg, 'components') and cfg.components:
        return [OmegaConf.to_container(c, resolve=True) for c in cfg.components]

    if hasattr(cfg, 'analysis') and hasattr(cfg.analysis, 'components') and cfg.analysis.components:
        return [OmegaConf.to_container(c, resolve=True) for c in cfg.analysis.components]

    raise ValueError(f"No components configured for model '{resolved_key}' or in analysis config.")


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


def get_fit_data_dir(cfg: DictConfig) -> Path:
    """Return fit_data directory path from Hydra config."""
    raw_path = getattr(cfg.infrastructure.output_paths, 'fit_data_dir', 'fit_data')
    return Path(resolve_data_path(cfg, raw_path))


def get_model_storage_coords(cfg: DictConfig, model_key: Optional[str] = "default") -> Dict[str, Any]:
    """Extract model storage coordinates (ROI, threshold, coord, bin_size) from config."""
    resolved_key = _get_model_key(cfg, model_key)
    model = cfg.models[resolved_key]
    return {
        'roi': (float(model.roi.esum_min), float(model.roi.esum_max)),
        'threshold': float(model.roi.threshold),
        'coord': getattr(model, 'coord', 'h2uv'),
        'bin_size': float(model.binning.bin_size),
    }
