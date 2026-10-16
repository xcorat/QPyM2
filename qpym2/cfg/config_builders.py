"""
Configuration builder functions for QPyM2.

Converts Hydra DictConfig objects into dict structures expected by
runner scripts (run_create_staging.py, run_create_input.py, run_fits.py).

These builders provide the bridge between the new Hydra YAML config system
and the existing function signatures that expect specific dict structures.
"""

import math
from typing import Dict, Any, List, Tuple, Optional
from omegaconf import DictConfig, OmegaConf
from pathlib import Path

from qpym2.utils.hist_cut_generator import for_qroot


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
    
    Args:
        cfg: Hydra DictConfig object
        model_key: Model configuration key (default: first model in cfg.models)
        
    Returns:
        Dict with keys: outpath, treename, out_cols, cuts
        
    Example:
        >>> cfg = get_config(environment='tech_note')
        >>> staging_cfg = build_staging_config(cfg)
        >>> staging_cfg['outpath']
        '/global/cfs/.../staging_1ton_v2'
    """
    # Get model configuration
    model = cfg.models[_get_model_key(cfg, model_key)]
    
    # Generate ROI cuts using hist_cut_generator
    roi_params = OmegaConf.to_container(model.roi, resolve=True)
    roi_cuts = for_qroot(
        esum_var=cfg.qroot.esum_var,
        evar=cfg.qroot.energy_var,
        **roi_params
    )
    
    return {
        'outpath': cfg.infrastructure.output_paths.staging_dir,
        'treename': cfg.qroot.output_tree_name,
        'out_cols': list(cfg.qroot.output_columns),
        'cuts': roi_cuts
    }


def build_mc_config(cfg: DictConfig) -> Dict[str, Any]:
    """
    Build background MC configuration dict from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        
    Returns:
        Dict with keys: inpath, treename, evar, esum_var, ch1, ch2, ds, t, extra_cuts
        
    Example:
        >>> cfg = get_config(environment='tech_note')
        >>> mc_cfg = build_mc_config(cfg)
        >>> mc_cfg['inpath']
        '/global/cfs/.../out_ares'
    """
    multiplicity_cut = (
        cfg.analysis_cuts.get('multiplicity', 'Multiplicity == 2')
        if 'analysis_cuts' in cfg else 'Multiplicity == 2'
    )
    base_cuts = cfg.cuts.get('base_cuts', 'true') if 'cuts' in cfg else 'true'
    mc_cuts = cfg.cuts.get('mc_cuts', None) if 'cuts' in cfg else None
    
    return {
        'inpath': cfg.infrastructure.directories.bkg_mc_dir,
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
    
    Args:
        cfg: Hydra DictConfig object
        model_key: Model configuration key (default: first model in cfg.models)
        
    Returns:
        Dict with keys: inpath, fname, out_fname, treename, evar, esum_var,
                       ch1, ch2, ds, t, extra_cuts
        
    Example:
        >>> cfg = get_config(environment='tech_note')
        >>> ndbd_cfg = build_ndbd_config(cfg)
        >>> ndbd_cfg['out_fname']
        'ndbd.root'
    """
    # Get model configuration
    model = cfg.models[_get_model_key(cfg, model_key)]
    
    # Get base cuts and model-specific cuts
    multiplicity_cut = (
        cfg.analysis_cuts.get('multiplicity', 'Multiplicity == 2')
        if 'analysis_cuts' in cfg else 'Multiplicity == 2'
    )
    base_cuts = cfg.cuts.get('base_cuts', 'true') if 'cuts' in cfg else 'true'
    mc_cuts = cfg.cuts.get('mc_cuts', None) if 'cuts' in cfg else None
    ndbd_cuts = model.cuts.get('ndbd_cuts', None)
    
    # Get NDBD file pattern from config
    ndbd_pattern = cfg.infrastructure.get('file_patterns', {}).get('ndbd_pattern', '0vbb_1M_*.root')
    
    # Get naming convention
    ndbd_staging_name = cfg.infrastructure.naming_conventions.ndbd_staging
    
    return {
        'inpath': cfg.infrastructure.directories.ndbd_mc_dir,
        'fname': ndbd_pattern,
        'out_fname': f'{ndbd_staging_name}.root',
        'treename': cfg.qroot.tree_name,
        'evar': cfg.qroot.energy_var,
        'esum_var': cfg.qroot.esum_var,
        'ch1': cfg.qroot.channel_vars[0],
        'ch2': cfg.qroot.channel_vars[1],
        'ds': 'Dataset',
        't': 'Time',
        'extra_cuts': _merge_cuts(multiplicity_cut, base_cuts, mc_cuts, ndbd_cuts)
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
        blinded: If True, use blinded data name
        model_key: Model configuration key (unused; reserved for future use)
        
    Returns:
        Dict with keys: inpath, fname, out_fname, treename, evar, esum_var,
                       ch1, ch2, ds, t, extra_cuts
        
    Example:
        >>> cfg = get_config(environment='tech_note')
        >>> data_cfg = build_data_config(cfg)
        >>> data_cfg['out_fname']
        'data.root'
    """
    # Get base/data cuts with legacy-compatible defaults
    multiplicity_cut = (
        cfg.analysis_cuts.get('multiplicity', 'Multiplicity == 2')
        if 'analysis_cuts' in cfg else 'Multiplicity == 2'
    )
    base_cuts = cfg.cuts.get('base_cuts', 'true') if 'cuts' in cfg else 'true'
    data_cuts = cfg.cuts.get('data_cuts', None) if 'cuts' in cfg else None
    
    # Get data file name
    reduced_fname = cfg.infrastructure.data_files.reduced_data
    
    # Get naming convention
    if blinded:
        out_name = cfg.infrastructure.naming_conventions.data_blinded_name
    else:
        out_name = cfg.infrastructure.naming_conventions.data_name
    
    return {
        'inpath': cfg.infrastructure.directories.bm_dir,
        'fname': reduced_fname,
        'out_fname': f'{out_name}.root',
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
    """
    Extract JAGS HDF5 file path from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        
    Returns:
        Path to JAGS output file
        
    Example:
        >>> cfg = get_config(environment='tech_note')
        >>> jags_path = get_jags_path(cfg)
        >>> jags_path
        '/global/cfs/.../Out.h5'
    """
    return cfg.infrastructure.data_files.jags_output


def get_data_names(cfg: DictConfig) -> Tuple[str, str]:
    """
    DEPRECATED: use the `get_data_name instead.
    
    Get data naming conventions from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        
    Returns:
        Tuple of (data_name, data_blinded_name)
        
    Example:
        >>> cfg = get_config(environment='tech_note')
        >>> data_name, blinded_name = get_data_names(cfg)
        >>> data_name
        'data'
    """
    naming = cfg.infrastructure.naming_conventions
    return (naming.data_name, naming.data_blinded_name)


def get_staging_dir(cfg: DictConfig) -> str:
    """
    Extract staging directory path from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        
    Returns:
        Path to staging directory
    """
    return cfg.infrastructure.output_paths.staging_dir


def get_bkg_mc_dir(cfg: DictConfig) -> str:
    """
    Extract background MC directory from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        
    Returns:
        Path to background MC directory
    """
    return cfg.infrastructure.directories.bkg_mc_dir


def get_ndbd_mc_dir(cfg: DictConfig) -> str:
    """
    Extract NDBD MC directory from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        
    Returns:
        Path to NDBD MC directory
    """
    return cfg.infrastructure.directories.ndbd_mc_dir


def build_input_config(cfg: DictConfig, model_key: Optional[str] = None) -> Dict[str, Any]:
    """
    Build fit input creation configuration dict from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        model_key: Model configuration key (default: first model in cfg.models)
        
    Returns:
        Dict with input creation parameters
    """
    # Get staging directory
    staging_dir = cfg.infrastructure.output_paths.staging_dir
    
    # Get model configuration
    model = cfg.models[_get_model_key(cfg, model_key)]

    return {
        # outpath of the input config section refers to the staging directory
        'outpath': staging_dir,
        'treename': cfg.qroot.output_tree_name,
        'data_fname': f"{cfg.infrastructure.naming_conventions.data_name}.root",
        'roi': OmegaConf.to_container(model.roi, resolve=True),
        'binning': OmegaConf.to_container(model.binning, resolve=True),
        'data_blinded_name': cfg.infrastructure.naming_conventions.data_blinded_name
    }


def build_write_config(cfg: DictConfig, model_key: Optional[str] = None) -> Dict[str, Any]:
    """
    Build fit output write configuration dict from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        model_key: Model configuration key (unused; reserved for future use)
        
    Returns:
        Dict with write configuration parameters
    """
    # Get fit directory
    fit_dir = cfg.infrastructure.output_paths.fit_dir
    
    return {
        'outdir': fit_dir,
        'overwrite': True,
    }


def get_h2uv_model(cfg: DictConfig, model_key: Optional[str] = None):
    """
    Build a TempConfig histogram model from Hydra config.
    
    Creates the TempConfig object with all derived attributes needed by
    downstream code (hists.get_empty_hist, fit.setup_fit, create_m2uv, etc.):
    umin/umax/vmin/vmax (√2 transformation), binsize, hist_type, smooth, etc.
    
    Args:
        cfg: Hydra DictConfig object
        model_key: Model configuration key (default: first model in cfg.models)
        
    Returns:
        TempConfig with flat attributes: esum_min, esum_max, threshold,
        binsize, umin, umax, vmin, vmax, smooth, skip_roi, hist_type, name
    """
    from qpym2.utils import TempConfig
    from qpym2.utils.hists import HistType

    resolved_key = _get_model_key(cfg, model_key)
    model = cfg.models[resolved_key]
    s2 = math.sqrt(2)

    # ROI parameters
    esum_min = model.roi.esum_min
    esum_max = model.roi.esum_max
    threshold = model.roi.threshold
    if 'u_range' in model.binning or 'v_range' in model.binning:
        raise ValueError(
            f"Model '{resolved_key}' uses deprecated binning.u_range/v_range. "
            "Active h2uv configs derive rotated-axis ranges from roi/esum and threshold. "
            "Please remove u_range/v_range and keep only binning.bin_size (defined in e1/e2 space)."
        )

    # bin_size in YAML is defined in e1/e2-space (keV); scale to rotated u/v-space (keV/√2)
    # so that get_hist_settings computes correct nbins:
    #   ceil((umax - umin) / binsize) = ceil(ΔESum / bin_size) = same # of bins
    binsize = model.binning.bin_size / s2

    # Derived u/v ranges (rotated coordinates)
    umin = esum_min / s2
    umax = esum_max / s2
    vmin = 0
    vmax = (esum_max - 2 * threshold) / s2

    # Smoothing: tuple (iterations, kernel) or False
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

    # skip_roi
    skip_roi = getattr(cfg.hyperparams, 'skip_roi', None) if hasattr(cfg, 'hyperparams') else None

    # Model name (for ROOT histogram naming)
    name = getattr(model, 'name', f'hmu_{model_key}')

    return TempConfig(
        esum_min=esum_min,
        esum_max=esum_max,
        threshold=threshold,
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


def get_comp_cfg(cfg: DictConfig) -> List[Dict[str, str]]:
    """
    Read background component grouping configuration from Hydra config.
    
    Returns the list-of-dicts format expected by make_bkgcomps_table,
    make_priors, and fit.setup_fit.
    
    Args:
        cfg: Hydra DictConfig with a ``components`` key
        
    Returns:
        List of dicts, each with 'name' and either 'filter_regex' or
        'filter_regex_inverse'
    """
    if not hasattr(cfg, 'components'):
        # Fallback to hardcoded default (backward compat)
        from qpym2.fit.fit import comp_cfg
        return [dict(c) for c in comp_cfg]

    return [OmegaConf.to_container(c, resolve=True) for c in cfg.components]


def get_input_fname(cfg: DictConfig, fitdir: str, skip_roi: bool) -> str:
    """Build the fit input HDF5 path from config naming conventions.

    Args:
        cfg: Hydra DictConfig object.
        fitdir: Directory containing the fit input file.
        skip_roi: If True, return the no-ROI variant of the input file.

    Returns:
        Full path to the input HDF5 file.
    """
    naming = cfg.infrastructure.naming_conventions
    base = naming.get('input_noroi_name', 'input_noroi') if skip_roi else naming.get('input_name', 'input')
    return f"{fitdir}/{base}.h5"


def get_corrmat_path(cfg: DictConfig) -> Optional[str]:
    """Extract background systematic correlation matrix path from Hydra config.

    Returns None if not set, allowing callers to skip systematic sampling.
    """
    return getattr(cfg.infrastructure.data_files, 'bkg_syst_corrmat', None)


def get_fit_dir(cfg: DictConfig) -> str:
    """
    Extract fit directory from Hydra config.
    
    Args:
        cfg: Hydra DictConfig object
        
    Returns:
        Path to fit directory
    """
    return cfg.infrastructure.output_paths.fit_dir


def get_data_dir(cfg: DictConfig) -> Optional[str]:
    """Return the base results data_dir, or None if not set."""
    val = getattr(cfg.infrastructure.output_paths, 'data_dir', None)
    return str(val) if val else None


def get_plots_dir(cfg: DictConfig) -> str:
    """Return plots directory.

    Priority:
    1. ``data_dir/plots`` if ``data_dir`` is set in config
    2. Explicit ``plots_dir`` config key
    3. ``fit_dir/plots`` as fallback
    """
    data_dir = get_data_dir(cfg)
    if data_dir:
        return str(Path(data_dir) / 'plots')
    plots = getattr(cfg.infrastructure.output_paths, 'plots_dir', None)
    if plots:
        return str(plots)
    return str(Path(get_fit_dir(cfg)) / 'plots')


def get_syst_dir(cfg: DictConfig) -> str:
    """Return systematic results directory.

    Priority:
    1. ``data_dir/syst_results`` if ``data_dir`` is set
    2. Explicit ``systematic_dir`` config key
    """
    data_dir = get_data_dir(cfg)
    if data_dir:
        return str(Path(data_dir) / 'syst_results')
    syst = getattr(cfg.infrastructure.output_paths, 'systematic_dir', None)
    if syst:
        return str(syst)
    raise ValueError(
        "No systematic results directory configured. "
        "Set data_dir or systematic_dir in output_paths."
    )


def get_toy_results_dir(cfg: DictConfig) -> str:
    """Return toy results directory.

    Priority:
    1. ``data_dir/toy_results`` if ``data_dir`` is set
    2. Explicit ``toy_results_dir`` config key
    """
    data_dir = get_data_dir(cfg)
    if data_dir:
        return str(Path(data_dir) / 'toy_results')
    toy = getattr(cfg.infrastructure.output_paths, 'toy_results_dir', None)
    if toy:
        return str(toy)
    raise ValueError(
        "No toy results directory configured. "
        "Set data_dir or toy_results_dir in output_paths."
    )


# NOTE: build_input_config and build_write_config are defined above (single definitions).
