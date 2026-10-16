"""
Systematics configuration for QPyM2.

Provides LSS systematic sample definitions and path retrieval
for systematic analysis outputs.
"""

# ---------------------------------------------------------------------------
# Default LSS parameters (fallback when not present in Hydra config)
# ---------------------------------------------------------------------------

DEFAULT_LSS_DATASETS = range(3601, 3616, 1)

DEFAULT_LSS_SYST_SAMPLES = {
    'sscale_1r': {'type': 'sscl', 'spt':  1},
    'dqres_1r':  {'type': 'dq',   'ptn': '1r'},
    'dqres_1l':  {'type': 'dq',   'ptn': '1l'},
}

# LSS file patterns (used with lss_data_path)
# {ds} is replaced with dataset number, {lss_data_path} with the base path
DEFAULT_LSS_SCALING_PATTERN = 'residual_and_width_vs_energy_ds{ds}.root'
DEFAULT_LSS_PEAKSHAPE_PATTERN = 'PeakShape_ds*.root'


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def get_systematic_dir(cfg):
    """Get systematic output directory from config.

    Resolution order:
    1. ``data_dir/syst_results`` if ``data_dir`` is set in output_paths
    2. Explicit ``systematic_dir`` key in output_paths

    Args:
        cfg: Hydra configuration object.

    Returns:
        str: Path to systematic output directory.

    Raises:
        ValueError: If neither data_dir nor systematic_dir is configured.
    """
    from qpym2.cfg.config_builders import get_syst_dir
    return get_syst_dir(cfg)


def get_lss_data_path(cfg):
    """Get LSS data path from config.

    Args:
        cfg: Hydra configuration object.

    Returns:
        str: Path to LSS ROOT data files.

    Raises:
        ValueError: If cfg.infrastructure.data_files.lss_data_path is not set.
    """
    try:
        return cfg.infrastructure.data_files.lss_data_path
    except (AttributeError, Exception) as e:
        if 'not in struct' in str(e) or isinstance(e, AttributeError):
            raise ValueError(
                "cfg.infrastructure.data_files.lss_data_path is not set. "
                "Please configure it in your environment YAML or use --lss-data-path."
            ) from e
        raise


def get_lss_file_patterns(cfg):
    """Get LSS file patterns from config.

    Args:
        cfg: Hydra configuration object.

    Returns:
        tuple: (scaling_pattern, peakshape_pattern)
    """
    try:
        patterns = cfg.infrastructure.data_files.lss_file_patterns
        scaling = getattr(patterns, 'scaling', DEFAULT_LSS_SCALING_PATTERN)
        peakshape = getattr(patterns, 'peakshape', DEFAULT_LSS_PEAKSHAPE_PATTERN)
        return scaling, peakshape
    except AttributeError:
        return DEFAULT_LSS_SCALING_PATTERN, DEFAULT_LSS_PEAKSHAPE_PATTERN


def build_lss_file_paths(lss_data_path, datasets, cfg=None):
    """Build LSS file paths from base path and datasets.

    Args:
        lss_data_path: Base path to LSS ROOT files.
        datasets: List/range of dataset numbers.
        cfg: Optional Hydra config for custom patterns.

    Returns:
        tuple: (scaling_fnames, peakshape_glob)
            - scaling_fnames: list of paths per dataset
            - peakshape_glob: glob pattern for peakshape files
    """
    import os

    if cfg is not None:
        scaling_pat, peakshape_pat = get_lss_file_patterns(cfg)
    else:
        scaling_pat = DEFAULT_LSS_SCALING_PATTERN
        peakshape_pat = DEFAULT_LSS_PEAKSHAPE_PATTERN

    scaling_fnames = [
        os.path.join(lss_data_path, scaling_pat.format(ds=ds))
        for ds in datasets
    ]
    peakshape_glob = os.path.join(lss_data_path, peakshape_pat)

    return scaling_fnames, peakshape_glob
