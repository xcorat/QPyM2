"""
Utilities for constructing model directory names from Hydra configurations.

Provides backward compatibility with the old model_cfg naming scheme.
"""

from typing import Optional, Dict, Any
from pathlib import Path
from omegaconf import DictConfig

from qpym2.cfg.hydra_adapter_new import resolve_data_path


def get_shiftname(shift_value: float, shift_type: str = 'ushiftp') -> str:
    """Generate shift name from value."""
    return f"{shift_type}{shift_value*10:02.0f}"


def build_model_dir(
    base_name: str = 'model',
    co_cut: Optional[str] = None,
    ndbd_cutoff: Optional[str] = None,
    smooth: bool = True
) -> str:
    """Build model directory name (without shift parameters)."""
    mname = base_name
    if co_cut:
        mname = mname + '_' + co_cut
    if ndbd_cutoff:
        mname = mname + '_' + ndbd_cutoff
    if not smooth:
        mname = mname + '_nosmooth'
    return mname


def get_shift_suffix(co_shift: float = 0.0, ndbd_shift: float = 0.0) -> str:
    """Build filename suffix from shift parameters."""
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
    """Build full model name from parameters (directory + shift suffix)."""
    return build_model_dir(base_name, co_cut, ndbd_cutoff, smooth) + \
           get_shift_suffix(co_shift, ndbd_shift)


def get_model_dir_from_config(cfg: DictConfig, model_key: Optional[str] = None) -> Path:
    """Extract model directory name from Hydra configuration."""
    hyperparams = cfg.get('hyperparams', {})
    co_cut = hyperparams.get('co_cuts', None)
    ndbd_cutoff = hyperparams.get('ndbd_cutoff', None)
    smooth = hyperparams.get('smooth', True)
    
    model_dir = build_model_dir(
        base_name='model',
        co_cut=co_cut,
        ndbd_cutoff=ndbd_cutoff,
        smooth=smooth
    )
    
    raw_fit_dir = cfg.infrastructure.output_paths.fit_dir
    fit_dir = Path(resolve_data_path(cfg, raw_fit_dir))
    return fit_dir / model_dir


def parse_model_name(model_name: str) -> Dict[str, Any]:
    """Parse model directory name back into components."""
    parts = model_name.split('_')
    result = {
        'base': parts[0] if parts else 'model',
        'co_cut': None,
        'ndbd_cutoff': None,
        'co_shift': 0.0,
        'ndbd_shift': 0.0,
        'smooth': True
    }
    
    i = 1
    co_parts = []
    while i < len(parts):
        part = parts[i]
        if part.startswith('u') and len(part) > 1 and part[1:].isdigit():
            result['ndbd_cutoff'] = part
            i += 1
            break
        elif part.startswith('ushiftp') or part.startswith('ndbd'):
            break
        elif part == 'nosmooth':
            break
        else:
            co_parts.append(part)
            i += 1
    
    if co_parts:
        result['co_cut'] = '_'.join(co_parts)
    
    while i < len(parts):
        part = parts[i]
        if part.startswith('ushiftp'):
            shift_str = part.replace('ushiftp', '')
            result['co_shift'] = float(shift_str) / 10.0 if shift_str else 0.0
        elif part.startswith('ndbdushiftp'):
            shift_str = part.replace('ndbdushiftp', '')
            result['ndbd_shift'] = float(shift_str) / 10.0 if shift_str else 0.0
        elif part == 'nosmooth':
            result['smooth'] = False
        i += 1
    
    return result
