"""
Storage backend module for qpym2.

Implements the path and deterministic hashing storage architecture described in
docs/ai/reviews/2026-10-29-data-storage-review.md.
"""

from qpym2.storage.schema import (
    normalize_float,
    normalize_cut_string,
    build_template_hash_dict,
    build_fit_hash_dict,
)
from qpym2.storage.hasher import (
    canonical_json,
    hash_dict,
    compute_input_uuid,
    compute_output_uuid,
)
from qpym2.storage.paths import (
    format_roi_dir,
    format_thr_dir,
    format_coord_bin_dir,
    get_parameter_dir,
    get_input_dir,
    get_output_dir,
    get_template_path,
    get_data_path,
    get_fit_path,
    get_study_dir,
    get_study_paths,
)
from qpym2.storage.metadata import (
    write_hdf5_metadata,
    read_hdf5_metadata,
)
from qpym2.storage.backend import (
    write_template,
    read_template,
    write_observed_data,
    read_observed_data,
    write_fit,
    read_fit,
    template_exists,
    data_exists,
    fit_exists,
)
from qpym2.storage.study import (
    init_study_db,
    record_study_point,
    write_study_manifest,
    read_study_manifest,
    collect_study_results,
)
from qpym2.storage.fits import (
    resolve_fit_template_uuids,
    assemble_fit_components,
    build_fit_hash,
    setup_storage_fit,
)

__all__ = [
    # Schema
    "normalize_float",
    "normalize_cut_string",
    "build_template_hash_dict",
    "build_fit_hash_dict",
    # Hasher
    "canonical_json",
    "hash_dict",
    "compute_input_uuid",
    "compute_output_uuid",
    # Paths
    "format_roi_dir",
    "format_thr_dir",
    "format_coord_bin_dir",
    "get_parameter_dir",
    "get_input_dir",
    "get_output_dir",
    "get_template_path",
    "get_data_path",
    "get_fit_path",
    "get_study_dir",
    "get_study_paths",
    # Metadata
    "write_hdf5_metadata",
    "read_hdf5_metadata",
    # Backend
    "write_template",
    "read_template",
    "write_observed_data",
    "read_observed_data",
    "write_fit",
    "read_fit",
    "template_exists",
    "data_exists",
    "fit_exists",
    # Study
    "init_study_db",
    "record_study_point",
    "write_study_manifest",
    "read_study_manifest",
    "collect_study_results",
    # Fits
    "resolve_fit_template_uuids",
    "assemble_fit_components",
    "build_fit_hash",
    "setup_storage_fit",
]
