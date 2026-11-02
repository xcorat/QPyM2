"""
Path generation and directory structure routing for qpym2 storage backend.

Implements the path contract:
  fit_data/
    roi_<min>_<max>/thr<val>/<coord>_b<bin_size>/
      input/
        data.h5
        data_blinded.h5
        <component_name>_<input_UUID>.h5
      output/
        fit_<output_UUID>.h5
        fit_<output_UUID>.nc
    studies/
      <study_name>/
        <study_name>.manifest
        <study_name>.sqlite
        results.csv
        plots/
"""

from pathlib import Path
from typing import Dict, Tuple, Union


def format_roi_dir(esum_min: Union[int, float], esum_max: Union[int, float]) -> str:
    """Format ROI directory string, e.g. roi_2400_2570."""
    return f"roi_{int(esum_min)}_{int(esum_max)}"


def format_thr_dir(threshold: Union[int, float]) -> str:
    """Format threshold directory string, e.g. thr200."""
    thr_val = int(threshold) if float(threshold).is_integer() else threshold
    return f"thr{thr_val}"


def format_coord_bin_dir(coord: str, bin_size: Union[int, float]) -> str:
    """Format coordinate and binning directory string, e.g. h2uv_b2.5 or h2uv_b1."""
    return f"{coord}_b{bin_size:g}"


def get_parameter_dir(
    fit_data_root: Union[str, Path],
    roi: Tuple[Union[int, float], Union[int, float]],
    threshold: Union[int, float],
    coord: str,
    bin_size: Union[int, float],
) -> Path:
    """Return the base parameter coordinate directory Path."""
    root = Path(fit_data_root)
    roi_dir = format_roi_dir(roi[0], roi[1])
    thr_dir = format_thr_dir(threshold)
    cb_dir = format_coord_bin_dir(coord, bin_size)
    return root / roi_dir / thr_dir / cb_dir


def get_input_dir(param_dir: Union[str, Path]) -> Path:
    """Return the input directory under a parameter coordinate directory."""
    return Path(param_dir) / "input"


def get_output_dir(param_dir: Union[str, Path]) -> Path:
    """Return the output directory under a parameter coordinate directory."""
    return Path(param_dir) / "output"


def get_template_path(
    input_dir: Union[str, Path],
    component_name: str,
    input_uuid: str,
) -> Path:
    """Return path for a composite template HDF5 file (<component_name>_<input_UUID>.h5)."""
    return Path(input_dir) / f"{component_name}_{input_uuid}.h5"


def get_data_path(input_dir: Union[str, Path], blinded: bool = False) -> Path:
    """Return path for binned observed data (data.h5 or data_blinded.h5)."""
    fname = "data_blinded.h5" if blinded else "data.h5"
    return Path(input_dir) / fname


def get_fit_path(
    output_dir: Union[str, Path],
    output_uuid: str,
    ext: str = ".h5",
) -> Path:
    """Return path for fit result file (fit_<output_UUID>.h5 or fit_<output_UUID>.nc)."""
    clean_ext = ext if ext.startswith('.') else f".{ext}"
    return Path(output_dir) / f"fit_{output_uuid}{clean_ext}"


def get_study_dir(fit_data_root: Union[str, Path], study_name: str) -> Path:
    """Return base directory for a study."""
    return Path(fit_data_root) / "studies" / study_name


def get_study_paths(study_dir: Union[str, Path], study_name: str) -> Dict[str, Path]:
    """Return dict of standard file paths for a study."""
    sdir = Path(study_dir)
    return {
        'dir': sdir,
        'manifest': sdir / f"{study_name}.manifest",
        'sqlite': sdir / f"{study_name}.sqlite",
        'results_csv': sdir / "results.csv",
        'plots_dir': sdir / "plots",
    }
