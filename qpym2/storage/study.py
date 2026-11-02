"""
Study tracking, indexing, and result collection routines for qpym2 storage backend.

Manages:
- <study_name>.sqlite: Central collected index for parameter points and artifact references.
- <study_name>.manifest: Metadata and parameter sweep definition.
- results.csv: Aggregated summary metrics across study points.
"""

from datetime import datetime
import json
from pathlib import Path
import sqlite3
from typing import Any, Dict, List, Optional, Union
import pandas as pd


def _ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def init_study_db(sqlite_path: Union[str, Path]) -> Path:
    """
    Initialize SQLite index database for a study.

    Args:
        sqlite_path: Path to <study_name>.sqlite.

    Returns:
        Path to initialized database.
    """
    p = Path(sqlite_path)
    _ensure_parent_dir(p)

    with sqlite3.connect(p) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS study_points (
                point_id TEXT PRIMARY KEY,
                study_name TEXT NOT NULL,
                task_id INTEGER,
                roi_dir TEXT,
                thr_dir TEXT,
                coord_bin_dir TEXT,
                output_uuid TEXT,
                status TEXT NOT NULL,
                summary_json TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
            """
        )
        cursor.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_study_status
            ON study_points (study_name, status)
            """
        )
        cursor.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_output_uuid
            ON study_points (output_uuid)
            """
        )
        conn.commit()

    return p


def record_study_point(
    sqlite_path: Union[str, Path],
    point_dict: Dict[str, Any],
) -> None:
    """
    Insert or update a study point record in SQLite index.

    Args:
        sqlite_path: Path to <study_name>.sqlite.
        point_dict: Dictionary containing point fields (point_id, study_name, output_uuid, status, etc.).
    """
    p = Path(sqlite_path)
    init_study_db(p)

    now_iso = datetime.utcnow().isoformat()
    point_id = str(point_dict.get("point_id", point_dict.get("output_uuid", "")))
    if not point_id:
        raise ValueError("point_dict must contain 'point_id' or 'output_uuid'")

    summary_json = point_dict.get("summary_json")
    if summary_json is not None and not isinstance(summary_json, str):
        summary_json = json.dumps(summary_json, sort_keys=True)

    with sqlite3.connect(p) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT INTO study_points (
                point_id, study_name, task_id, roi_dir, thr_dir, coord_bin_dir,
                output_uuid, status, summary_json, created_at, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(point_id) DO UPDATE SET
                task_id=excluded.task_id,
                roi_dir=excluded.roi_dir,
                thr_dir=excluded.thr_dir,
                coord_bin_dir=excluded.coord_bin_dir,
                output_uuid=excluded.output_uuid,
                status=excluded.status,
                summary_json=excluded.summary_json,
                updated_at=excluded.updated_at
            """,
            (
                point_id,
                str(point_dict.get("study_name", "")),
                point_dict.get("task_id"),
                point_dict.get("roi_dir"),
                point_dict.get("thr_dir"),
                point_dict.get("coord_bin_dir"),
                point_dict.get("output_uuid"),
                str(point_dict.get("status", "pending")),
                summary_json,
                point_dict.get("created_at", now_iso),
                now_iso,
            ),
        )
        conn.commit()


def write_study_manifest(
    manifest_path: Union[str, Path],
    manifest_data: Dict[str, Any],
) -> Path:
    """
    Write study specification and sweep definition manifest.

    Args:
        manifest_path: Path to <study_name>.manifest (JSON).
        manifest_data: Dictionary describing the study.

    Returns:
        Path to manifest file.
    """
    p = Path(manifest_path)
    _ensure_parent_dir(p)

    if "created_at" not in manifest_data:
        manifest_data["created_at"] = datetime.utcnow().isoformat()

    with open(p, "w", encoding="utf-8") as f:
        json.dump(manifest_data, f, indent=2, sort_keys=True)

    return p


def read_study_manifest(manifest_path: Union[str, Path]) -> Dict[str, Any]:
    """Read study manifest JSON file."""
    p = Path(manifest_path)
    if not p.exists():
        raise FileNotFoundError(f"Study manifest not found: {p}")
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


def collect_study_results(
    study_dir: Union[str, Path],
    study_name: str,
    result_csv_paths: Optional[List[Union[str, Path]]] = None,
) -> pd.DataFrame:
    """
    Collect, validate, and merge worker task result CSVs into results.csv.

    Also synchronizes the study SQLite index database.

    Args:
        study_dir: Base study directory (fit_data/studies/<study_name>/).
        study_name: Name of the study.
        result_csv_paths: Optional explicit list of task CSV files to merge.
            If None, scans study_dir / "tasks" for *.csv files.

    Returns:
        Merged results pandas DataFrame.
    """
    sdir = Path(study_dir)
    sqlite_path = sdir / f"{study_name}.sqlite"
    out_csv = sdir / "results.csv"

    csv_files = result_csv_paths
    if csv_files is None:
        tasks_dir = sdir / "tasks"
        csv_files = list(tasks_dir.glob("*.csv")) if tasks_dir.exists() else []

    if not csv_files:
        empty_df = pd.DataFrame()
        empty_df.to_csv(out_csv, index=False)
        return empty_df

    dfs = []
    for f in csv_files:
        fp = Path(f)
        if fp.exists() and fp.stat().st_size > 0:
            try:
                df = pd.read_csv(fp)
                dfs.append(df)
            except Exception as e:
                print(f"Warning: could not read task result CSV {fp}: {e}")

    if not dfs:
        empty_df = pd.DataFrame()
        empty_df.to_csv(out_csv, index=False)
        return empty_df

    merged_df = pd.concat(dfs, ignore_index=True)
    _ensure_parent_dir(out_csv)
    merged_df.to_csv(out_csv, index=False)

    # Sync into SQLite index if sqlite exists
    if sqlite_path.exists():
        for _, row in merged_df.iterrows():
            output_uuid = str(row.get("output_uuid", row.get("task_id", "")))
            if output_uuid:
                record_study_point(
                    sqlite_path,
                    {
                        "point_id": output_uuid,
                        "study_name": study_name,
                        "task_id": row.get("task_id"),
                        "output_uuid": row.get("output_uuid"),
                        "status": "completed",
                        "summary_json": row.to_dict(),
                    },
                )

    return merged_df
