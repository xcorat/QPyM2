"""
Schema definitions and parameter normalization for qpym2 storage hashing.

Defines standardized, versioned dictionary structures for:
- Template inputs (<component_name>_<input_UUID>.h5)
- Fit outputs (fit_<output_UUID>.h5 / fit_<output_UUID>.nc)
"""

from typing import Dict, Any, Optional, Union


def normalize_float(val: Optional[Union[float, int]], precision: int = 6) -> Optional[float]:
    """Round float to fixed precision to avoid IEEE-754 string serialization discrepancies."""
    if val is None:
        return None
    return round(float(val), precision)


def normalize_cut_string(cut: Optional[str]) -> Optional[str]:
    """Normalize a cut expression by stripping whitespace and standardizing boolean tokens."""
    if cut is None:
        return None
    c = str(cut).strip()
    if not c or c.lower() == 'true':
        return None
    tokens = [part.strip() for part in c.split('&&') if part.strip() and part.strip().lower() != 'true']
    if not tokens:
        return None
    return ' && '.join(tokens)


def build_template_hash_dict(
    component: Optional[str] = None,
    staging_id: str = "v1",
    cuts: Optional[Dict[str, Optional[str]]] = None,
    transforms: Optional[Dict[str, Any]] = None,
    grouping: Optional[Dict[str, Any]] = None,
    norm: Optional[str] = None,
    schema_version: int = 1,
) -> Dict[str, Any]:
    """
    Construct canonical dictionary representing template parameter variations.

    The component name is excluded from the hash payload so that all sources
    sharing identical staging, cuts, and coordinate transforms produce the
    exact same parameter-set <input_UUID>.

    Args:
        component: Optional component name (ignored for hashing; retained for backward compatibility).
        staging_id: Identifier or contract version of the staging dataset.
        cuts: Dictionary of cuts affecting this variation (e.g. {'common': 'v>150', 'specific': None}).
        transforms: Dictionary of coordinate/energy transforms (e.g. {'shift_type': 'ushiftp', 'shift_val': 0.8}).
        grouping: Grouping rules (e.g. {'filter_regex': '60Co'}).
        norm: JAGS norm mode used for component weighting ('mode', 'mean', etc.).
        schema_version: Schema version for forward/backward compatibility.

    Returns:
        Canonical parameter dictionary ready for deterministic hashing.
    """
    normalized_cuts: Dict[str, Optional[str]] = {}
    if cuts:
        for k, v in sorted(cuts.items()):
            normalized_cuts[k] = normalize_cut_string(v)

    normalized_transforms: Optional[Dict[str, Any]] = None
    if transforms:
        normalized_transforms = {}
        for k, v in sorted(transforms.items()):
            if isinstance(v, (int, float)):
                normalized_transforms[k] = normalize_float(v)
            else:
                normalized_transforms[k] = v

    normalized_grouping: Optional[Dict[str, Any]] = None
    if grouping:
        normalized_grouping = {k: grouping[k] for k in sorted(grouping.keys())}

    return {
        'schema_version': int(schema_version),
        'staging_id': str(staging_id),
        'cuts': normalized_cuts,
        'transforms': normalized_transforms,
        'grouping': normalized_grouping,
        'norm': str(norm) if norm is not None else None,
    }


def build_fit_hash_dict(
    model_name: str,
    data_source: str,
    templates: Dict[str, str],
    smoothing: Optional[Dict[str, Any]] = None,
    priors: Optional[Dict[str, Any]] = None,
    roi_settings: Optional[Dict[str, Any]] = None,
    mcmc: Optional[Dict[str, Any]] = None,
    study_context: Optional[Dict[str, Any]] = None,
    schema_version: int = 1,
) -> Dict[str, Any]:
    """
    Construct canonical dictionary representing a fit configuration.

    Args:
        model_name: Fit formulation ('m2', 'm2_binned', etc.).
        data_source: Observed data source ('data', 'data_blinded').
        templates: Mapping of component names to their input UUIDs: {'co60': 'a3f8...', ...}.
        smoothing: Smoothing configuration (enabled, method, kernel).
        priors: Prior parameters (sigma_scale, use_sample_std, sample_mkchain).
        roi_settings: Likelihood / ROI exclusion settings (e.g. {'skip_roi': [2515, 2540]}).
        mcmc: MCMC sampling parameters (chains, draws, tune, target_accept).
        study_context: Optional study-specific settings (e.g. {'injected_n': 10, 'seed': 47}).
        schema_version: Schema version.

    Returns:
        Canonical parameter dictionary ready for deterministic hashing.
    """
    sorted_templates = {k: str(templates[k]) for k in sorted(templates.keys())}

    norm_smoothing = {k: smoothing[k] for k in sorted(smoothing.keys())} if smoothing else None

    norm_priors = None
    if priors:
        norm_priors = {}
        for k, v in sorted(priors.items()):
            if isinstance(v, (int, float)):
                norm_priors[k] = normalize_float(v)
            else:
                norm_priors[k] = v

    norm_roi = None
    if roi_settings:
        norm_roi = {}
        for k, v in sorted(roi_settings.items()):
            if isinstance(v, (list, tuple)):
                norm_roi[k] = [normalize_float(x) if isinstance(x, (int, float)) else x for x in v]
            elif isinstance(v, (int, float)):
                norm_roi[k] = normalize_float(v)
            else:
                norm_roi[k] = v

    norm_mcmc = {k: mcmc[k] for k in sorted(mcmc.keys())} if mcmc else None

    norm_study = None
    if study_context:
        norm_study = {}
        for k, v in sorted(study_context.items()):
            if isinstance(v, (int, float)):
                norm_study[k] = normalize_float(v)
            else:
                norm_study[k] = v

    return {
        'schema_version': int(schema_version),
        'model_name': str(model_name),
        'data_source': str(data_source),
        'templates': sorted_templates,
        'smoothing': norm_smoothing,
        'priors': norm_priors,
        'roi_settings': norm_roi,
        'mcmc': norm_mcmc,
        'study_context': norm_study,
    }
