"""
Fit assembly and orchestration routines for QPyM2 storage architecture (v0.6).

Bridges source-level unsmoothed HDF5 templates stored in:
fit_data/roi_<min>_<max>/thr<val>/<coord>_b<bin_size>/input/
into likelihood components, applies downstream model smoothing and JAGS normalization,
computes canonical <output_UUID>s, and configures PyMC Bayesian fits.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
from omegaconf import DictConfig

from qpym2.cfg.config_builders_new import (
    get_comp_cfg,
    get_corrmat_path,
    get_jags_path,
)
from qpym2.fit import fit as fit_module
from qpym2.fit.params import make_priors
from qpym2.io.bkg_model import get_mcsyst_samples, read_jagsh5, read_mkchain
from qpym2.models import m2
from qpym2.storage.backend import (
    fit_exists,
    read_observed_data,
    read_template,
)
from qpym2.storage.hasher import compute_input_uuid, compute_output_uuid
from qpym2.storage.inputs import resolve_storage_coords
from qpym2.storage.paths import (
    get_data_path,
    get_output_dir,
    get_template_path,
)
from qpym2.storage.schema import build_fit_hash_dict, build_template_hash_dict

_CONAME = '60Co'
_NDBDNAME = 'ndbd'


def resolve_fit_template_uuids(
    cfg: DictConfig,
    model_key: Optional[str] = None,
    co_cut: Optional[str] = None,
    ndbd_cutoff: Optional[str] = None,
    co_shift: float = 0.0,
    ndbd_shift: float = 0.0,
) -> Tuple[Dict[str, str], str, str, str]:
    """
    Resolve parameter-set UUIDs (base, co_shift, ndbd_shift) and map
    analysis likelihood components to their parameter-set UUID.

    Args:
        cfg: Hydra configuration DictConfig.
        model_key: Named model key (e.g. 'tech_note_200').
        co_cut: Co-60 calibration cut key or expression.
        ndbd_cutoff: Signal lower-edge cutoff key or expression.
        co_shift: Co-60 energy scale shift in keV.
        ndbd_shift: Signal energy scale shift in keV.

    Returns:
        Tuple of (component_uuid_map, uuid_base, uuid_co, uuid_ndbd).
    """
    staging_id = str(getattr(getattr(cfg, 'staging', None), 'version', 'v1'))
    co_cuts_dict = dict(getattr(cfg.cuts, 'co_cuts', {}))
    ndbd_cuts_dict = dict(getattr(cfg.cuts, 'ndbd_cutoff', {}))

    common_cuts: Dict[str, Optional[str]] = {}
    if co_cut and co_cut in co_cuts_dict:
        common_cuts['co_cut'] = co_cuts_dict[co_cut]
    elif co_cut and co_cut != 'None':
        common_cuts['co_cut'] = co_cut

    if ndbd_cutoff and ndbd_cutoff in ndbd_cuts_dict:
        common_cuts['ndbd_cutoff'] = ndbd_cuts_dict[ndbd_cutoff]
    elif ndbd_cutoff and ndbd_cutoff != 'None':
        common_cuts['ndbd_cutoff'] = ndbd_cutoff

    # 1. Base UUID (invariant parameters)
    p_base = build_template_hash_dict(staging_id=staging_id, cuts=common_cuts, transforms=None)
    uuid_base = compute_input_uuid(p_base)

    # 2. Co-60 shifted UUID
    if co_shift != 0.0:
        transforms_co = {'shift_type': 'ushiftp', 'shift_val': co_shift}
        p_co = build_template_hash_dict(staging_id=staging_id, cuts=common_cuts, transforms=transforms_co)
        uuid_co = compute_input_uuid(p_co)
    else:
        uuid_co = uuid_base

    # 3. NDBD shifted UUID
    if ndbd_shift != 0.0:
        transforms_ndbd = {'shift_type': 'ushiftp', 'shift_val': ndbd_shift}
        p_ndbd = build_template_hash_dict(staging_id=staging_id, cuts=common_cuts, transforms=transforms_ndbd)
        uuid_ndbd = compute_input_uuid(p_ndbd)
    else:
        uuid_ndbd = uuid_base

    # Map analysis components to their parameter-set UUID
    comp_cfg = get_comp_cfg(cfg)
    component_uuid_map: Dict[str, str] = {}
    for c in comp_cfg:
        cname = c['name']
        f_regex = c.get('filter_regex', '')
        if '60Co' in f_regex or 'co60' in cname.lower():
            component_uuid_map[cname] = uuid_co
        else:
            component_uuid_map[cname] = uuid_base

    # Signal
    component_uuid_map['ndbd'] = uuid_ndbd

    return component_uuid_map, uuid_base, uuid_co, uuid_ndbd


def assemble_fit_components(
    input_dir: Union[str, Path],
    jags_fpath: Union[str, Path],
    comp_cfg: List[Dict[str, Any]],
    hm: Any,
    uuid_base: str,
    uuid_co: str,
    uuid_ndbd: str,
    smooth: Optional[Any] = None,
    read_jagsh5_fn: Optional[Any] = None,
    read_template_fn: Optional[Any] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[List[np.ndarray]], pd.DataFrame]:
    """
    Assemble fit components table from source-level unsmoothed HDF5 templates.

    Reads source templates, applies JAGS normalization and exposure weighting,
    groups into likelihood components, applies downstream smoothing, and formats
    the final component DataFrame for PyMC.

    Args:
        input_dir: Path to fit_data/<coords>/input/ directory.
        jags_fpath: Path to JAGS h5 file.
        comp_cfg: Component configuration list from analysis tier.
        hm: Histogram model (TempConfig).
        uuid_base: Baseline parameter-set UUID.
        uuid_co: Co-60 parameter-set UUID.
        uuid_ndbd: NDBD parameter-set UUID.
        smooth: Optional smoothing override (defaults to hm.smooth).
        read_jagsh5_fn: Custom JAGS reader for testing.
        read_template_fn: Custom template reader for testing.

    Returns:
        Tuple of (mctable_with_mchist, signal_df, edges, comps_df).
    """
    inp_dir = Path(input_dir)
    jreader = read_jagsh5_fn or read_jagsh5
    treader = read_template_fn or read_template

    mctable, _ = jreader(str(jags_fpath))

    # Read all source histograms
    edges_captured = None
    mchists = []
    for source_name in mctable.index:
        is_co = '60Co' in source_name or 'co60' in source_name.lower()
        uuid = uuid_co if is_co else uuid_base
        tpl_path = get_template_path(inp_dir, source_name, uuid)
        if not tpl_path.exists():
            raise FileNotFoundError(f"Missing template for source '{source_name}': {tpl_path}")
        hist, edges, _ = treader(tpl_path)
        mchists.append(hist)
        if edges_captured is None and edges is not None:
            edges_captured = edges

    mctable['mchist'] = mchists

    # Signal
    ndbd_path = get_template_path(inp_dir, _NDBDNAME, uuid_ndbd)
    if not ndbd_path.exists():
        raise FileNotFoundError(f"Missing signal template: {ndbd_path}")
    sig_hist, sig_edges, _ = treader(ndbd_path)
    if edges_captured is None and sig_edges is not None:
        edges_captured = sig_edges

    sig_integral = float(np.sum(sig_hist))
    signal_df = pd.DataFrame(
        {
            'mchist': [sig_hist],
            'integral': [sig_integral],
            'prior': [None],
        },
        index=['ndbd'],
    )

    # Downstream smoothing applied during fit setup
    smoothing_param = smooth if smooth is not None else getattr(hm, 'smooth', None)

    # Ensure comp_cfg has 'std' if get_comps is called
    safe_comp_cfg = []
    for c in comp_cfg:
        c_copy = dict(c)
        if 'std' not in c_copy:
            c_copy['std'] = 0.0
        safe_comp_cfg.append(c_copy)

    # Build likelihood components
    comps_df = m2.get_comps(
        mctable=mctable,
        signal_df=signal_df,
        comp_cfg=safe_comp_cfg,
        smooth=smoothing_param,
    )

    return mctable, signal_df, edges_captured, comps_df


def build_fit_hash(
    model_name: str,
    data_source: str,
    template_uuids: Dict[str, str],
    smoothing: Optional[Dict[str, Any]] = None,
    priors: Optional[Dict[str, Any]] = None,
    roi_settings: Optional[Dict[str, Any]] = None,
    mcmc: Optional[Dict[str, Any]] = None,
    study_context: Optional[Dict[str, Any]] = None,
) -> Tuple[Dict[str, Any], str]:
    """
    Construct canonical fit hash dictionary and compute deterministic 12-char output UUID.

    Args:
        model_name: Model formulation ('m2', 'm2_binned').
        data_source: Data source name ('data', 'data_blinded').
        template_uuids: Mapping of component names to parameter-set UUIDs.
        smoothing: Smoothing settings dictionary.
        priors: Priors settings dictionary.
        roi_settings: Likelihood / ROI exclusion settings.
        mcmc: Sampler configuration (chains, draws, tune).
        study_context: Optional study parameter coordinates.

    Returns:
        Tuple of (fit_params_dict, output_uuid).
    """
    fit_params = build_fit_hash_dict(
        model_name=model_name,
        data_source=data_source,
        templates=template_uuids,
        smoothing=smoothing,
        priors=priors,
        roi_settings=roi_settings,
        mcmc=mcmc,
        study_context=study_context,
    )
    output_uuid = compute_output_uuid(fit_params)
    return fit_params, output_uuid


def setup_storage_fit(
    cfg: DictConfig,
    model_key: Optional[str] = None,
    co_cut: Optional[str] = None,
    ndbd_cutoff: Optional[str] = None,
    co_shift: float = 0.0,
    ndbd_shift: float = 0.0,
    blinded: Union[bool, str] = False,
    skip_roi: Optional[Union[List[float], Tuple[float, float]]] = None,
    prior_width: str = '5sig',
    mkchain_samples: int = 100,
    chains: int = 8,
    draws: int = 20000,
    tune: int = 1000,
    target_accept: float = 0.9,
    study_context: Optional[Dict[str, Any]] = None,
    model_name: str = 'm2',
    mksamples: Optional[pd.DataFrame] = None,
    read_jagsh5_fn: Optional[Any] = None,
    read_chain_fn: Optional[Any] = None,
    read_template_fn: Optional[Any] = None,
    read_data_fn: Optional[Any] = None,
    dry_run: bool = False,
) -> Dict[str, Any]:
    """
    Orchestrate full setup for a Bayesian fit in the QPyM2 storage architecture.

    Resolves coordinates, template UUIDs, checks file readiness, loads inputs,
    computes priors, constructs PyMC model, and calculates output UUID and paths.

    Returns:
        Dictionary containing:
        - coords, param_dir, input_dir, output_dir
        - output_uuid, fit_prefix, h5_path, nc_path
        - fit_params (metadata dictionary)
        - is_cached (bool: whether fit_<output_UUID>.h5 already exists)
        - pymc_model, fit_comps, data, edges
    """
    coords, hm, param_dir, input_dir = resolve_storage_coords(cfg, model_key=model_key)
    out_dir = get_output_dir(param_dir)
    jags_fpath = get_jags_path(cfg)
    comp_cfg = [dict(c) for c in get_comp_cfg(cfg)]

    # 1. Resolve parameter-set UUIDs
    tpl_uuids, uuid_base, uuid_co, uuid_ndbd = resolve_fit_template_uuids(
        cfg=cfg,
        model_key=model_key,
        co_cut=co_cut,
        ndbd_cutoff=ndbd_cutoff,
        co_shift=co_shift,
        ndbd_shift=ndbd_shift,
    )

    # 2. Build canonical fit parameter hash and determine output UUID
    smoothing_cfg = {
        'enabled': bool(getattr(hm, 'smooth', None)),
        'method': 'smooth_nph2',
        'smooth': list(hm.smooth) if getattr(hm, 'smooth', None) else None,
    }

    if prior_width in ('syst', '2syst'):
        sigma_scale = 1.0 if prior_width == 'syst' else 2.0
        use_sample_std = True
    else:
        sigma_scale = 5.0 if prior_width == '5sig' else 10.0
        use_sample_std = False

    prior_cfg = {
        'prior_width': prior_width,
        'sigma_scale': sigma_scale,
        'use_sample_std': use_sample_std,
        'mkchain_samples': mkchain_samples,
    }

    roi_cfg = {'skip_roi': list(skip_roi) if skip_roi else None}
    mcmc_cfg = {
        'chains': chains,
        'draws': draws,
        'tune': tune,
        'target_accept': target_accept,
    }
    data_source_name = 'data_blinded' if blinded else 'data'

    fit_params, output_uuid = build_fit_hash(
        model_name=model_name,
        data_source=data_source_name,
        template_uuids=tpl_uuids,
        smoothing=smoothing_cfg,
        priors=prior_cfg,
        roi_settings=roi_cfg,
        mcmc=mcmc_cfg,
        study_context=study_context,
    )

    fit_prefix = out_dir / f"fit_{output_uuid}"
    h5_path = fit_prefix.with_suffix('.h5')
    nc_path = fit_prefix.with_suffix('.nc')
    cached = fit_exists(fit_prefix)

    if dry_run:
        return {
            'coords': coords,
            'param_dir': param_dir,
            'input_dir': input_dir,
            'output_dir': out_dir,
            'template_uuids': tpl_uuids,
            'uuid_base': uuid_base,
            'uuid_co': uuid_co,
            'uuid_ndbd': uuid_ndbd,
            'output_uuid': output_uuid,
            'fit_prefix': fit_prefix,
            'h5_path': h5_path,
            'nc_path': nc_path,
            'is_cached': cached,
            'fit_params': fit_params,
            'pymc_model': None,
            'fit_comps': None,
            'data': None,
            'edges': None,
            'mcmc_cfg': mcmc_cfg,
        }

    # 3. Read observed data
    d_reader = read_data_fn or read_observed_data
    data_path = get_data_path(input_dir, blinded=bool(blinded))
    if not data_path.exists():
        raise FileNotFoundError(f"Observed data file not found: {data_path}")
    data_events, edges_data, _ = d_reader(data_path)

    # 4. Assemble components and compute priors
    mctable, signal_df, edges_tpl, _ = assemble_fit_components(
        input_dir=input_dir,
        jags_fpath=jags_fpath,
        comp_cfg=comp_cfg,
        hm=hm,
        uuid_base=uuid_base,
        uuid_co=uuid_co,
        uuid_ndbd=uuid_ndbd,
        read_jagsh5_fn=read_jagsh5_fn,
        read_template_fn=read_template_fn,
    )

    edges = edges_tpl if edges_tpl is not None else edges_data

    # Priors generation
    if mksamples is None:
        chain_reader = read_chain_fn or read_mkchain
        chain = chain_reader(str(jags_fpath))
        if use_sample_std:
            corrmat_fpath = str(get_corrmat_path(cfg))
            mksamples = get_mcsyst_samples(chain, corrmat_fpath=corrmat_fpath, n=mkchain_samples).T
        else:
            mksamples = chain.sample(n=mkchain_samples).T

    comp_cfg_priors = make_priors(
        comp_cfg=comp_cfg,
        mctable=mctable,
        mksamples=mksamples,
        sigma_scale=sigma_scale,
        use_sample_std=use_sample_std,
    )

    # PyMC model construction
    pymc_model, fit_comps = fit_module.setup_fit(
        mctable=mctable,
        data=data_events,
        signal_df=signal_df,
        hm=hm,
        smooth=hm.smooth,
        comp_cfg=comp_cfg_priors,
        edges=edges,
    )

    return {
        'coords': coords,
        'param_dir': param_dir,
        'input_dir': input_dir,
        'output_dir': out_dir,
        'template_uuids': tpl_uuids,
        'uuid_base': uuid_base,
        'uuid_co': uuid_co,
        'uuid_ndbd': uuid_ndbd,
        'output_uuid': output_uuid,
        'fit_prefix': fit_prefix,
        'h5_path': h5_path,
        'nc_path': nc_path,
        'is_cached': cached,
        'fit_params': fit_params,
        'pymc_model': pymc_model,
        'fit_comps': fit_comps,
        'data': data_events,
        'edges': edges,
        'mcmc_cfg': mcmc_cfg,
    }
