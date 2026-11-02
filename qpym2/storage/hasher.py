"""
Deterministic hashing routines for qpym2 storage backend.

Produces reproducible, collision-resistant UUIDs from canonical JSON representations
of resolved parameters.
"""

import hashlib
import json
from typing import Any, Dict


def canonical_json(obj: Any) -> str:
    """
    Serialize Python object into canonical, whitespace-free, key-sorted JSON string.

    Args:
        obj: Dictionary or primitive to serialize.

    Returns:
        Canonical JSON string.
    """
    return json.dumps(
        obj,
        sort_keys=True,
        separators=(',', ':'),
        ensure_ascii=True,
        default=lambda o: str(o)
    )


def hash_dict(data: Dict[str, Any], length: int = 12) -> str:
    """
    Generate deterministic hex hash from dictionary data.

    Args:
        data: Parameter dictionary.
        length: Desired character length of hash digest (default 12).

    Returns:
        Truncated hex digest string.
    """
    canonical_str = canonical_json(data)
    digest = hashlib.sha256(canonical_str.encode('utf-8')).hexdigest()
    return digest[:length]


def compute_input_uuid(template_params: Dict[str, Any], length: int = 12) -> str:
    """
    Compute deterministic <input_UUID> for a template parameter configuration.

    The UUID identifies the parameter set (staging contract, cuts, coordinate transforms).
    Sources sharing the same parameters share the same <input_UUID>.

    Args:
        template_params: Template parameters (from schema.build_template_hash_dict).
        length: Hash length (default 12).

    Returns:
        12-character hex UUID.
    """
    return hash_dict(template_params, length=length)


def compute_output_uuid(fit_params: Dict[str, Any], length: int = 12) -> str:
    """
    Compute deterministic <output_UUID> for a fit result.

    Args:
        fit_params: Fit parameters (from schema.build_fit_hash_dict).
        length: Hash length (default 12).

    Returns:
        12-character hex UUID.
    """
    return hash_dict(fit_params, length=length)
