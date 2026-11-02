"""
Input generation routines for QPyM2 storage architecture (v0.6).

Transforms staged ROOT files into identity-addressed HDF5 inputs:
- Source-level background MC templates: <component_name>_<input_UUID>.h5
- Signal template: ndbd_<input_UUID>.h5
- Binned / event observed data: data.h5 and data_blinded.h5

Stored under:
fit_data/roi_<min>_<max>/thr<val>/<coord>_b<bin_size>/input/

Per docs/ai/reviews/2026-10-29-data-storage-review.md:
- Inputs are persisted at source granularity (1-to-1 with staging ROOT files).
- Templates store original, unsmoothed histograms.
- Component grouping and smoothing take place downstream during fit setup.
- Shifts (e.g. co_shift on 60Co) only modify the UUID of the affected source;
  invariant sources are cached and reused across shifts.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
from omegaconf import DictConfig

from qpym2.cfg.config_builders_new import (
    get_fit_data_dir,
    get_model_storage_coords,
    get_h2uv_model,
    get_jags_path,
    get_staging_dir,
)
from qpym2.io import staging
from qpym2.io.bkg_model import read_jagsh5
from qpym2.storage.paths import (
    get_parameter_dir,
    get_input_dir,
    get_template_path,
    get_data_path,
)
from qpym2.storage.schema import build_template_hash_dict
from qpym2.storage.hasher import compute_input_uuid
from qpym2.storage.backend import (
    write_template,
    template_exists,
    write_observed_data,
    data_exists,
)
from qpym2.utils.hist_model import get_hm_cuts
from qpym2.utils.logging import debug

_CONAME = '60Co'
_NDBDNAME = 'ndbd'


def resolve_storage_coords(
    cfg: DictConfig,
    model_key: Optional[str] = None,
) -> Tuple[Dict[str, Any], Any, Path, Path]:
    """
    Resolve model storage coordinates, histogram model, and destination directories.

    Args:
        cfg: Hydra configuration DictConfig.
        model_key: Named model key (e.g. 'tech_note_200').

    Returns:
        Tuple of (coords_dict, h2uv_model, param_dir, input_dir).
    """
    coords = get_model_storage_coords(cfg, model_key=model_key)
    hm = get_h2uv_model(cfg, model_key=model_key)
    fit_data_root = get_fit_data_dir(cfg)

    param_dir = get_parameter_dir(
        fit_data_root,
        roi=coords['roi'],
        threshold=coords['threshold'],
        coord=coords['coord'],
        bin_size=coords['bin_size'],
    )
    input_dir = get_input_dir(param_dir)

    return coords, hm, param_dir, input_dir


def build_shift_defs(
    shift_val: float,
    shift_type: str = 'ushiftp',
    invars: Tuple[str, str] = ('u0', 'v0'),
    outvars: Tuple[str, str] = ('u', 'v'),
) -> List[Tuple[str, str]]:
    """
    Build RDF branch definition tuples for coordinate transforms.

    Args:
        shift_val: Energy shift in keV.
        shift_type: Transform type ('ushiftp' or 'none').
        invars: Input column names ('u0', 'v0').
        outvars: Output column names ('u', 'v').

    Returns:
        List of (alias, definition_string) tuples for RDataFrame.
    """
    u_out, v_out = outvars
    u_in, v_in = invars

    if shift_val == 0.0 or shift_type == 'none':
        return [(u_out, u_in), (v_out, v_in)]

    if shift_type == 'ushiftp':
        de1 = f'({u_in}+{v_in})/2*{shift_val / 1000.0}/sqrt(2)'
        return [
            (u_out, f'{u_in} + {de1}'),
            (v_out, f'{v_in} + {de1}'),
        ]

    raise NotImplementedError(f"Unsupported shift type: {shift_type}")


def build_source_hash_params(
    component: Optional[str] = None,
    staging_id: str = "v1",
    cuts: Optional[Dict[str, Optional[str]]] = None,
    transforms: Optional[Dict[str, Any]] = None,
    grouping: Optional[Dict[str, Any]] = None,
    norm: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Construct canonical dictionary for template parameter hashing.

    The component name is excluded from the hash payload so that all sources
    sharing identical staging, cuts, and coordinate transforms produce the
    exact same parameter-set <input_UUID>.

    Args:
        component: Optional source name (ignored for hashing; retained for backward compatibility).
        staging_id: Staging dataset identifier.
        cuts: Active cuts dictionary affecting this source.
        transforms: Transforms affecting this source.
        grouping: Optional source grouping metadata.
        norm: JAGS norm mode ('mode', 'mean').

    Returns:
        Canonical parameter dictionary.
    """
    return build_template_hash_dict(
        component=component,
        staging_id=staging_id,
        cuts=cuts,
        transforms=transforms,
        grouping=grouping,
        norm=norm,
    )


def generate_source_template(
    source_name: str,
    staging_fpath: Union[str, Path],
    treename: str,
    hm: Any,
    cuts_dict: Dict[str, Optional[str]],
    transforms_dict: Optional[Dict[str, Any]] = None,
    staging_id: str = "v1",
    input_dir: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    hist_cache: Optional[Dict[Any, Any]] = None,
    read_hist_fn: Optional[Any] = None,
) -> Tuple[Path, str, bool]:
    """
    Generate or reuse an unsmoothed source-level template histogram from staging.

    Args:
        source_name: Staged source name (e.g. '60Co', '208Tl', 'ndbd').
        staging_fpath: Path to staged ROOT file.
        treename: Tree name in ROOT file ('uvtree').
        hm: Histogram model (TempConfig).
        cuts_dict: Dictionary of cuts affecting this source.
        transforms_dict: Dictionary of transforms (e.g. {'shift_type': 'ushiftp', 'shift_val': 0.8}).
        staging_id: Staging version string.
        input_dir: Destination input directory.
        overwrite: If True, regenerate and overwrite existing file.
        hist_cache: In-memory histogram cache across shifts.
        read_hist_fn: Optional custom reader (for dependency injection / testing).

    Returns:
        Tuple of (target_path, input_uuid, is_cached).
    """
    tpl_params = build_source_hash_params(
        component=source_name,
        staging_id=staging_id,
        cuts=cuts_dict,
        transforms=transforms_dict,
    )
    input_uuid = compute_input_uuid(tpl_params)

    if input_dir is None:
        raise ValueError("input_dir must be provided")
    target_path = get_template_path(input_dir, source_name, input_uuid)

    if not overwrite and template_exists(target_path):
        debug(f"[CACHE HIT] Template exists: {target_path.name}")
        return target_path, input_uuid, True

    # Build branch definitions
    shift_val = float(transforms_dict.get('shift_val', 0.0)) if transforms_dict else 0.0
    shift_type = str(transforms_dict.get('shift_type', 'ushiftp')) if transforms_dict else 'none'
    defs = build_shift_defs(shift_val=shift_val, shift_type=shift_type)

    # Build filters
    hm_cut = get_hm_cuts(hm, inpars=('u', 'v'))
    filters = [hm_cut]
    for cut_val in cuts_dict.values():
        if cut_val and str(cut_val).strip() and str(cut_val).strip().lower() != 'true':
            filters.append(str(cut_val).strip())

    # In-memory caching check
    cache_key = (str(staging_fpath), tuple(defs), tuple(filters))
    if hist_cache is not None and cache_key in hist_cache:
        debug(f"[MEM CACHE] Hit for {source_name}")
        hist, edges = hist_cache[cache_key]
    else:
        reader = read_hist_fn if read_hist_fn is not None else staging.read_hist
        read_result = reader(
            str(staging_fpath),
            hm,
            treename=treename,
            defs=defs,
            filters=filters,
            rtype='numpy',
        )
        hist = np.asarray(read_result[0], dtype=float)
        edges = list(read_result[1:])
        # Ensure no negative bins
        if np.any(hist < 0):
            debug(f"Clamping negative bins to 0 for {source_name}")
            hist[hist < 0] = 0.0

        if hist_cache is not None:
            hist_cache[cache_key] = (hist, edges)

    meta_to_write = dict(tpl_params)
    meta_to_write['component'] = source_name
    write_template(
        target_path,
        hist=hist,
        edges=edges,
        metadata=meta_to_write,
        overwrite=overwrite,
    )
    debug(f"[WRITE] Persisted template: {target_path.name}")
    return target_path, input_uuid, False


def generate_observed_data(
    data_fpath: Union[str, Path],
    treename: str,
    hm: Any,
    cuts_dict: Dict[str, Optional[str]],
    blinded: bool = False,
    staging_id: str = "v1",
    input_dir: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
    hist_cache: Optional[Dict[Any, Any]] = None,
    read_evlist_fn: Optional[Any] = None,
) -> Tuple[Path, bool]:
    """
    Generate or reuse observed data (data.h5 or data_blinded.h5).

    Args:
        data_fpath: Path to staged ROOT data file.
        treename: Tree name in ROOT file ('uvtree').
        hm: Histogram model (TempConfig).
        cuts_dict: Active cuts dictionary for data.
        blinded: If True, destination is data_blinded.h5.
        staging_id: Staging dataset version.
        input_dir: Destination input directory.
        overwrite: If True, overwrite existing file.
        hist_cache: In-memory cache.
        read_evlist_fn: Optional custom reader (for dependency injection / testing).

    Returns:
        Tuple of (target_path, is_cached).
    """
    if input_dir is None:
        raise ValueError("input_dir must be provided")
    target_path = get_data_path(input_dir, blinded=blinded)

    if not overwrite and data_exists(target_path):
        debug(f"[CACHE HIT] Observed data exists: {target_path.name}")
        return target_path, True

    defs = [('u', 'u0'), ('v', 'v0')]
    hm_cut = get_hm_cuts(hm, inpars=('u', 'v'))
    filters = [hm_cut]
    for cut_val in cuts_dict.values():
        if cut_val and str(cut_val).strip() and str(cut_val).strip().lower() != 'true':
            filters.append(str(cut_val).strip())

    cache_key = (str(data_fpath), tuple(defs), tuple(filters))
    if hist_cache is not None and cache_key in hist_cache:
        data = hist_cache[cache_key]
    else:
        reader = read_evlist_fn if read_evlist_fn is not None else staging.read_evlist
        data = reader(
            str(data_fpath),
            treename,
            defs=defs,
            filters=filters,
            rtype='numpy',
        )
        data = np.asarray(data)
        if hist_cache is not None:
            hist_cache[cache_key] = data

    metadata = {
        'staging_id': str(staging_id),
        'blinded': bool(blinded),
        'cuts': {k: str(v) for k, v in cuts_dict.items() if v},
    }
    write_observed_data(
        target_path,
        data=data,
        metadata=metadata,
        overwrite=overwrite,
    )
    debug(f"[WRITE] Persisted observed data: {target_path.name}")
    return target_path, False


def create_storage_inputs_for_point(
    cfg: DictConfig,
    model_key: Optional[str] = None,
    co_cut: Optional[str] = 'gamma10_v150',
    ndbd_cutoff: Optional[str] = 'u1740',
    co_shift: float = 0.0,
    ndbd_shift: float = 0.0,
    blinded: bool = False,
    overwrite: bool = False,
    dry_run: bool = False,
    hist_cache: Optional[Dict[Any, Any]] = None,
    read_hist_fn: Optional[Any] = None,
    read_evlist_fn: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Orchestrate storage input creation for a single parameter point.

    Discovers all background MC sources from JAGS output and staging directory,
    processes NDBD signal, and processes observed data.

    Args:
        cfg: Hydra configuration DictConfig.
        model_key: Named model key (e.g. 'tech_note_200').
        co_cut: Co-60 cut key (e.g. 'gamma10_v150').
        ndbd_cutoff: NDBD cutoff key (e.g. 'u1740').
        co_shift: Co-60 energy shift in keV.
        ndbd_shift: NDBD energy shift in keV.
        blinded: If True, also process data_blinded.root.
        overwrite: If True, overwrite existing files.
        dry_run: If True, preview paths and UUIDs without writing files.
        hist_cache: In-memory cache shared across shifts.
        read_hist_fn: Custom reader for testing.
        read_evlist_fn: Custom reader for testing.

    Returns:
        Manifest dictionary with point metadata, template paths, UUIDs, and cache statuses.
    """
    coords, hm, param_dir, input_dir = resolve_storage_coords(cfg, model_key=model_key)
    staging_id = str(getattr(getattr(cfg, 'staging', None), 'version', 'v1'))
    staging_dir = Path(get_staging_dir(cfg))
    treename = str(cfg.qroot.output_tree_name)
    jags_fpath = get_jags_path(cfg)

    # Cut definitions
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

    # Load background source summary
    mctable, _ = read_jagsh5(jags_fpath)
    sources = list(mctable.index)

    point_manifest: Dict[str, Any] = {
        'coords': coords,
        'param_dir': str(param_dir),
        'input_dir': str(input_dir),
        'co_cut': co_cut,
        'ndbd_cutoff': ndbd_cutoff,
        'co_shift': co_shift,
        'ndbd_shift': ndbd_shift,
        'blinded': blinded,
        'templates': {},
        'data': {},
    }

    # Precompute parameter-level UUIDs
    p_base = build_source_hash_params(staging_id=staging_id, cuts=common_cuts, transforms=None)
    uuid_base = compute_input_uuid(p_base)

    if co_shift != 0.0:
        transforms_co = {'shift_type': 'ushiftp', 'shift_val': co_shift}
        p_co = build_source_hash_params(staging_id=staging_id, cuts=common_cuts, transforms=transforms_co)
        uuid_co = compute_input_uuid(p_co)
    else:
        transforms_co = None
        uuid_co = uuid_base

    if ndbd_shift != 0.0:
        transforms_ndbd = {'shift_type': 'ushiftp', 'shift_val': ndbd_shift}
        p_ndbd = build_source_hash_params(staging_id=staging_id, cuts=common_cuts, transforms=transforms_ndbd)
        uuid_ndbd = compute_input_uuid(p_ndbd)
    else:
        transforms_ndbd = None
        uuid_ndbd = uuid_base

    if dry_run:
        # Calculate expected UUIDs and paths without I/O
        for source_name in sources:
            is_co = '60Co' in source_name or 'co60' in source_name.lower()
            uuid = uuid_co if is_co else uuid_base
            path = get_template_path(input_dir, source_name, uuid)
            point_manifest['templates'][source_name] = {
                'path': str(path),
                'uuid': uuid,
                'cached': template_exists(path),
            }

        # NDBD
        path_ndbd = get_template_path(input_dir, _NDBDNAME, uuid_ndbd)
        point_manifest['templates'][_NDBDNAME] = {
            'path': str(path_ndbd),
            'uuid': uuid_ndbd,
            'cached': template_exists(path_ndbd),
        }

        # Data
        data_path = get_data_path(input_dir, blinded=blinded)
        point_manifest['data'] = {
            'path': str(data_path),
            'cached': data_exists(data_path),
        }
        return point_manifest

    # 1. Process observed data
    data_fname = 'data_blinded.root' if blinded else 'data.root'
    data_fpath = staging_dir / data_fname
    d_path, d_cached = generate_observed_data(
        data_fpath=data_fpath,
        treename=treename,
        hm=hm,
        cuts_dict=common_cuts,
        blinded=blinded,
        staging_id=staging_id,
        input_dir=input_dir,
        overwrite=overwrite,
        hist_cache=hist_cache,
        read_evlist_fn=read_evlist_fn,
    )
    point_manifest['data'] = {'path': str(d_path), 'cached': d_cached}

    # 2. Process background MC source templates
    for source_name, row in mctable.iterrows():
        fname = row['fname']
        staging_fpath = staging_dir / fname
        is_co = '60Co' in str(source_name) or 'co60' in str(source_name).lower()
        transforms = {'shift_type': 'ushiftp', 'shift_val': co_shift} if is_co and co_shift != 0.0 else None

        t_path, t_uuid, t_cached = generate_source_template(
            source_name=str(source_name),
            staging_fpath=staging_fpath,
            treename=treename,
            hm=hm,
            cuts_dict=common_cuts,
            transforms_dict=transforms,
            staging_id=staging_id,
            input_dir=input_dir,
            overwrite=overwrite,
            hist_cache=hist_cache,
            read_hist_fn=read_hist_fn,
        )
        point_manifest['templates'][str(source_name)] = {
            'path': str(t_path),
            'uuid': t_uuid,
            'cached': t_cached,
        }

    # 3. Process NDBD signal template
    ndbd_fpath = staging_dir / 'ndbd.root'
    ndbd_transforms = {'shift_type': 'ushiftp', 'shift_val': ndbd_shift} if ndbd_shift != 0.0 else None
    t_ndbd_path, t_ndbd_uuid, t_ndbd_cached = generate_source_template(
        source_name=_NDBDNAME,
        staging_fpath=ndbd_fpath,
        treename=treename,
        hm=hm,
        cuts_dict=common_cuts,
        transforms_dict=ndbd_transforms,
        staging_id=staging_id,
        input_dir=input_dir,
        overwrite=overwrite,
        hist_cache=hist_cache,
        read_hist_fn=read_hist_fn,
    )
    point_manifest['templates'][_NDBDNAME] = {
        'path': str(t_ndbd_path),
        'uuid': t_ndbd_uuid,
        'cached': t_ndbd_cached,
    }

    return point_manifest
