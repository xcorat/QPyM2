"""
Core storage backend routines for reading and writing templates, data, and fits.

Delegates storage layout and file access without modifying legacy qpym2/io modules.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import h5py
import numpy as np
import pandas as pd

from qpym2.storage.metadata import write_hdf5_metadata, read_hdf5_metadata

HIST_DATASET_KEY = "hist"
DATA_DATASET_KEY = "data"
EDGES_GROUP_KEY = "hist_edges"
COMPONENTS_KEY = "components"


def _ensure_parent_dir(path: Path) -> None:
    """Ensure the parent directory of a path exists."""
    path.parent.mkdir(parents=True, exist_ok=True)


def write_template(
    path: Union[str, Path],
    hist: np.ndarray,
    edges: Optional[List[np.ndarray]] = None,
    metadata: Optional[Dict[str, Any]] = None,
    overwrite: bool = False,
) -> Path:
    """
    Write a composite template histogram and edges to HDF5.

    Args:
        path: Destination path (<component_name>_<input_UUID>.h5).
        hist: Numpy array representing the 2D template histogram.
        edges: Optional list of bin edge arrays.
        metadata: Optional metadata dictionary to persist in HDF5 attributes.
        overwrite: If True, overwrite existing file.

    Returns:
        Path to written file.
    """
    p = Path(path)
    if p.exists() and not overwrite:
        raise FileExistsError(f"Template already exists: {p}. Pass overwrite=True to replace.")

    _ensure_parent_dir(p)

    with h5py.File(p, "w") as h5f:
        h5f.create_dataset(HIST_DATASET_KEY, data=hist)
        if edges is not None:
            egrp = h5f.create_group(EDGES_GROUP_KEY)
            for i, e in enumerate(edges):
                egrp.create_dataset(str(i), data=e)

    if metadata:
        write_hdf5_metadata(p, metadata)

    return p


def read_template(
    path: Union[str, Path],
) -> Tuple[np.ndarray, Optional[List[np.ndarray]], Dict[str, Any]]:
    """
    Read a composite template histogram, edges, and metadata from HDF5.

    Args:
        path: Path to template HDF5 file.

    Returns:
        Tuple of (hist, edges, metadata).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Template file not found: {p}")

    with h5py.File(p, "r") as h5f:
        hist = h5f[HIST_DATASET_KEY][:]
        edges = None
        if EDGES_GROUP_KEY in h5f:
            egrp = h5f[EDGES_GROUP_KEY]
            edges = [egrp[str(i)][:] for i in sorted(int(k) for k in egrp.keys())]

    metadata = read_hdf5_metadata(p)
    return hist, edges, metadata


def write_observed_data(
    path: Union[str, Path],
    data: np.ndarray,
    edges: Optional[List[np.ndarray]] = None,
    metadata: Optional[Dict[str, Any]] = None,
    overwrite: bool = False,
) -> Path:
    """
    Write observed data (or blinded data) array and optional edges to HDF5.

    Args:
        path: Destination path (data.h5 or data_blinded.h5).
        data: Numpy array of observed data.
        edges: Optional list of bin edge arrays.
        metadata: Optional metadata dictionary.
        overwrite: If True, overwrite existing file.

    Returns:
        Path to written file.
    """
    p = Path(path)
    if p.exists() and not overwrite:
        raise FileExistsError(f"Data file already exists: {p}. Pass overwrite=True to replace.")

    _ensure_parent_dir(p)

    with h5py.File(p, "w") as h5f:
        h5f.create_dataset(DATA_DATASET_KEY, data=data)
        if edges is not None:
            egrp = h5f.create_group(EDGES_GROUP_KEY)
            for i, e in enumerate(edges):
                egrp.create_dataset(str(i), data=e)

    if metadata:
        write_hdf5_metadata(p, metadata)

    return p


def read_observed_data(
    path: Union[str, Path],
) -> Tuple[np.ndarray, Optional[List[np.ndarray]], Dict[str, Any]]:
    """
    Read observed data array, edges, and metadata from HDF5.

    Args:
        path: Path to data.h5 or data_blinded.h5.

    Returns:
        Tuple of (data, edges, metadata).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Data file not found: {p}")

    with h5py.File(p, "r") as h5f:
        data = h5f[DATA_DATASET_KEY][:]
        edges = None
        if EDGES_GROUP_KEY in h5f:
            egrp = h5f[EDGES_GROUP_KEY]
            edges = [egrp[str(i)][:] for i in sorted(int(k) for k in egrp.keys())]

    metadata = read_hdf5_metadata(p)
    return data, edges, metadata


def write_fit(
    path_prefix: Union[str, Path],
    table: pd.DataFrame,
    trace: Optional[Any] = None,
    data: Optional[np.ndarray] = None,
    metadata: Optional[Dict[str, Any]] = None,
    save_trace: bool = False,
    overwrite: bool = False,
) -> Tuple[Path, Optional[Path]]:
    """
    Write fit results to <path_prefix>.h5 and optionally <path_prefix>.nc.

    Args:
        path_prefix: Destination file prefix without extension (e.g. output_dir / "fit_a3f81c9d04e2").
        table: Fitted components DataFrame.
        trace: ArviZ InferenceData posterior trace (optional).
        data: Observed data events array (optional).
        metadata: Metadata dictionary to persist.
        save_trace: If True, persist trace to <path_prefix>.nc via ArviZ.
        overwrite: Overwrite existing files if True.

    Returns:
        Tuple of (h5_path, nc_path_or_None).
    """
    prefix = Path(path_prefix)
    _ensure_parent_dir(prefix)

    h5_path = prefix.with_suffix(".h5")
    nc_path = prefix.with_suffix(".nc") if (save_trace and trace is not None) else None

    if h5_path.exists() and not overwrite:
        raise FileExistsError(f"Fit file already exists: {h5_path}")
    if nc_path and nc_path.exists() and not overwrite:
        raise FileExistsError(f"Trace file already exists: {nc_path}")

    # Write HDF5 components table and optional data
    table.to_hdf(h5_path, key=COMPONENTS_KEY, mode="w")

    if data is not None:
        s = pd.Series([p for p in data])
        s.to_hdf(h5_path, key=DATA_DATASET_KEY, mode="a")

    if metadata:
        write_hdf5_metadata(h5_path, metadata)

    # Optionally persist ArviZ NetCDF trace
    if nc_path is not None:
        import arviz as az
        az.to_netcdf(trace, str(nc_path))

    return h5_path, nc_path


def read_fit(
    path_prefix: Union[str, Path],
) -> Tuple[pd.DataFrame, Optional[Any], Optional[np.ndarray], Dict[str, Any]]:
    """
    Read fit results from <path_prefix>.h5 and <path_prefix>.nc (if trace exists).

    Args:
        path_prefix: Prefix or path to .h5 file.

    Returns:
        Tuple of (table, trace_or_None, data_or_None, metadata).
    """
    prefix = Path(path_prefix)
    if prefix.suffix in (".h5", ".nc"):
        prefix = prefix.with_suffix("")

    h5_path = prefix.with_suffix(".h5")
    nc_path = prefix.with_suffix(".nc")

    if not h5_path.exists():
        raise FileNotFoundError(f"Fit HDF5 file not found: {h5_path}")

    table = pd.read_hdf(h5_path, key=COMPONENTS_KEY)

    data = None
    try:
        data_s = pd.read_hdf(h5_path, key=DATA_DATASET_KEY)
        data = data_s.values
    except (KeyError, ValueError):
        pass

    metadata = read_hdf5_metadata(h5_path)

    trace = None
    if nc_path.exists():
        import arviz as az
        trace = az.from_netcdf(str(nc_path))

    return table, trace, data, metadata


def template_exists(path: Union[str, Path]) -> bool:
    """Check if template HDF5 file exists."""
    return Path(path).is_file()


def data_exists(path: Union[str, Path]) -> bool:
    """Check if observed data HDF5 file exists."""
    return Path(path).is_file()


def fit_exists(path_prefix: Union[str, Path], check_trace: bool = False) -> bool:
    """Check if fit results exist on disk."""
    prefix = Path(path_prefix)
    if prefix.suffix in (".h5", ".nc"):
        prefix = prefix.with_suffix("")

    h5_ok = prefix.with_suffix(".h5").is_file()
    if not h5_ok:
        return False
    if check_trace:
        return prefix.with_suffix(".nc").is_file()
    return True
