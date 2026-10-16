"""
Utilities for constructing model directory names from Hydra configurations.

Provides backward compatibility with the old model_cfg naming scheme.
"""

from typing import Optional, Dict, Any
from omegaconf import DictConfig


def get_shiftname(shift_value: float, shift_type: str = 'ushiftp') -> str:
    """
    Generate shift name from value.
    
    Args:
        shift_value: Shift value in keV (e.g., 1.0)
        shift_type: Type of shift (default: 'ushiftp')
    
    Returns:
        Shift name string (e.g., 'ushiftp10' for shift_value=1.0)
    """
    return f"{shift_type}{shift_value*10:02.0f}"


def build_model_dir(
    base_name: str = 'model',
    co_cut: Optional[str] = None,
    ndbd_cutoff: Optional[str] = None,
    smooth: bool = True
) -> str:
    """
    Build model directory name (without shift parameters).

    Directory groups all co_shift / ndbd_shift / noroi variants together.

    Args:
        base_name: Base model name (default: 'model')
        co_cut: Co-60 calibration cut key (e.g., 'gamma10_v150', 'v150')
        ndbd_cutoff: 0νββ signal cutoff key (e.g., 'u1740', 'u1775')
        smooth: Whether smoothing is enabled (default: True)

    Returns:
        Model directory name (e.g., 'model_gamma10_v150_u1740')

    Examples:
        >>> build_model_dir()
        'model'
        >>> build_model_dir(co_cut='gamma10_v150', ndbd_cutoff='u1740')
        'model_gamma10_v150_u1740'
    """
    mname = base_name

    if co_cut:
        mname = mname + '_' + co_cut
    if ndbd_cutoff:
        mname = mname + '_' + ndbd_cutoff
    if not smooth:
        mname = mname + '_nosmooth'

    return mname


def get_shift_suffix(co_shift: float = 0.0, ndbd_shift: float = 0.0) -> str:
    """
    Build filename suffix from shift parameters.

    Zero/None values produce no suffix (backward-compatible).
    Ordering: co_shift first, then ndbd_shift.

    Args:
        co_shift: Co-60 energy shift in keV (default: 0.0)
        ndbd_shift: 0νββ energy shift in keV (default: 0.0)

    Returns:
        Suffix string (e.g., '_ushiftp10', '_ushiftp10_ndbdushiftp05', or '')

    Examples:
        >>> get_shift_suffix(0.0)
        ''
        >>> get_shift_suffix(0.5)
        '_ushiftp05'
        >>> get_shift_suffix(1.0, 0.5)
        '_ushiftp10_ndbdushiftp05'
    """
    suffix = ''
    if co_shift:
        suffix += '_' + get_shiftname(co_shift)
    if ndbd_shift:
        suffix += '_ndbd' + get_shiftname(ndbd_shift)
    return suffix


def build_model_name(
    base_name: str = 'model',
    co_cut: Optional[str] = None,
    ndbd_cutoff: Optional[str] = None, 
    co_shift: float = 0.0,
    ndbd_shift: float = 0.0,
    smooth: bool = True
) -> str:
    """
    Build full model name from parameters (directory + shift suffix).

    Kept for backward compatibility and logging.
    
    Examples:
        >>> build_model_name()
        'model'
        >>> build_model_name(co_cut='v150')
        'model_v150'
        >>> build_model_name(co_cut='gamma10_v150', ndbd_cutoff='u1740', co_shift=1.0)
        'model_gamma10_v150_u1740_ushiftp10'
    """
    return build_model_dir(base_name, co_cut, ndbd_cutoff, smooth) + \
           get_shift_suffix(co_shift, ndbd_shift)


def get_model_dir_from_config(cfg: DictConfig, model_key: Optional[str] = None) -> str:
    """
    Extract model directory name from Hydra configuration.
    
    Args:
        cfg: Hydra DictConfig object
        model_key: Model configuration key (default: first model in cfg.models)
    
    Returns:
        Full path to model directory
    
    Example:
        >>> cfg = HydraConfigManager.get_config(environment='tech_note')
        >>> model_dir = get_model_dir_from_config(cfg)
        '/path/to/fits/model_gamma10_v150'
    """
    from pathlib import Path
    
    # Get hyperparameters
    hyperparams = cfg.get('hyperparams', {})
    co_cut = hyperparams.get('co_cuts', None)
    ndbd_cutoff = hyperparams.get('ndbd_cutoff', None)
    co_shift = hyperparams.get('co_shift', 0.0)
    ndbd_shift = hyperparams.get('ndbd_shift', 0.0)
    smooth = hyperparams.get('smooth', True)
    
    # Build model directory (without shifts)
    model_dir = build_model_dir(
        base_name='model',
        co_cut=co_cut,
        ndbd_cutoff=ndbd_cutoff,
        smooth=smooth
    )
    
    # Construct full path
    fit_dir = Path(cfg.infrastructure.output_paths.fit_dir)
    return fit_dir / model_dir


def parse_model_name(model_name: str) -> Dict[str, Any]:
    """
    Parse model directory name back into components.
    
    Args:
        model_name: Model directory name (e.g., 'model_gamma10_v150_u1740')
    
    Returns:
        Dictionary with parsed components
    
    Example:
        >>> parse_model_name('model_gamma10_v150_u1740_ushiftp10')
        {'base': 'model', 'co_cut': 'gamma10_v150', 'ndbd_cutoff': 'u1740', 
         'co_shift': 1.0, 'smooth': True}
    """
    parts = model_name.split('_')
    result = {
        'base': parts[0] if parts else 'model',
        'co_cut': None,
        'ndbd_cutoff': None,
        'co_shift': 0.0,
        'ndbd_shift': 0.0,
        'smooth': True
    }
    
    # Parse co_cut (can be multi-part like gamma10_v150)
    i = 1
    co_parts = []
    while i < len(parts):
        part = parts[i]
        # Check if this is a cutoff (starts with 'u'), shift, or other flag
        if part.startswith('u') and len(part) > 1 and part[1:].isdigit():
            result['ndbd_cutoff'] = part
            i += 1
            break
        elif part.startswith('ushiftp') or part.startswith('ndbd'):
            # This is a shift parameter
            break
        elif part == 'nosmooth':
            break
        else:
            co_parts.append(part)
            i += 1
    
    if co_parts:
        result['co_cut'] = '_'.join(co_parts)
    
    # Parse remaining parts (shifts, flags)
    while i < len(parts):
        part = parts[i]
        if part.startswith('ushiftp'):
            # Extract shift value (e.g., ushiftp10 → 1.0)
            shift_str = part.replace('ushiftp', '')
            result['co_shift'] = float(shift_str) / 10.0 if shift_str else 0.0
        elif part.startswith('ndbdushiftp'):
            shift_str = part.replace('ndbdushiftp', '')
            result['ndbd_shift'] = float(shift_str) / 10.0 if shift_str else 0.0
        elif part == 'nosmooth':
            result['smooth'] = False
        i += 1
    
    return result
