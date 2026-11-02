"""
HDF5 file metadata serialization and retrieval routines for qpym2 storage.

Saves full resolved parameters, timestamps, and provenance into HDF5 file-level attributes.
"""

import json
from pathlib import Path
from typing import Any, Dict, Union
import h5py

from qpym2.storage.hasher import canonical_json


METADATA_JSON_KEY = "metadata_json"


def write_hdf5_metadata(fpath: Union[str, Path], metadata: Dict[str, Any]) -> None:
    """
    Write metadata dictionary to HDF5 file-level attributes.

    Stores the full dictionary as a canonical JSON string under 'metadata_json'
    and copies top-level primitive values directly to attrs for quick inspection.

    Args:
        fpath: Path to existing HDF5 file.
        metadata: Metadata dictionary to persist.
    """
    p = Path(fpath)
    if not p.exists():
        raise FileNotFoundError(f"Cannot write metadata: file {p} does not exist")

    with h5py.File(p, "a") as h5f:
        # Write canonical JSON representation
        json_str = canonical_json(metadata)
        h5f.attrs[METADATA_JSON_KEY] = json_str

        # Also store top-level primitive attributes for convenient h5dump/h5ls inspection
        for k, v in metadata.items():
            if isinstance(v, (str, int, float, bool)):
                h5f.attrs[k] = v


def read_hdf5_metadata(fpath: Union[str, Path]) -> Dict[str, Any]:
    """
    Read metadata dictionary from HDF5 file-level attributes.

    Args:
        fpath: Path to HDF5 file.

    Returns:
        Metadata dictionary.
    """
    p = Path(fpath)
    if not p.exists():
        raise FileNotFoundError(f"Cannot read metadata: file {p} does not exist")

    with h5py.File(p, "r") as h5f:
        if METADATA_JSON_KEY in h5f.attrs:
            raw_json = h5f.attrs[METADATA_JSON_KEY]
            if isinstance(raw_json, bytes):
                raw_json = raw_json.decode("utf-8")
            return json.loads(raw_json)

        # Fallback: read all available attributes
        result = {}
        for k, v in h5f.attrs.items():
            if isinstance(v, bytes):
                result[k] = v.decode("utf-8")
            else:
                result[k] = v
        return result
