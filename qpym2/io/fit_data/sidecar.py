"""JSON sidecars for fit-data payloads."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Dict, Mapping, Sequence, Union

from qpym2.storage.hasher import hash_file


def build_sidecar(
    *,
    object_id: str,
    kind: str,
    identity: Mapping[str, Any],
    files: Sequence[Union[str, Path]],
    partition: Mapping[str, Any],
    metadata: Mapping[str, Any] | None = None,
    upstream: Mapping[str, str] | None = None,
    profile_version: int = 1,
) -> Dict[str, Any]:
    """Build sidecar metadata and checksums for a completed payload."""
    file_records = []
    for value in files:
        path = Path(value)
        file_records.append(
            {
                "name": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": hash_file(path),
            }
        )
    return {
        "schema_version": 1,
        "profile_version": int(profile_version),
        "object_id": str(object_id),
        "kind": str(kind),
        "identity": dict(identity),
        "partition": dict(partition),
        "upstream": dict(upstream or {}),
        "files": file_records,
        "metadata": dict(metadata or {}),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }


def write_sidecar(path: Union[str, Path], sidecar: Mapping[str, Any]) -> None:
    """Write formatted, deterministic JSON to a path."""
    Path(path).write_text(
        json.dumps(dict(sidecar), sort_keys=True, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def read_sidecar(path: Union[str, Path]) -> Dict[str, Any]:
    """Read and validate the basic structure of a JSON sidecar."""
    source = Path(path)
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read fit-data sidecar {source}: {exc}") from exc
    required = {"schema_version", "object_id", "kind", "identity", "files"}
    missing = required - set(value)
    if missing:
        raise ValueError(f"Fit-data sidecar {source} is missing fields: {sorted(missing)}")
    return value
