"""CSV I/O for toy and systematic fit results.

For concurrent SLURM array tasks, each task writes its own per-task CSV
via :func:`task_csv_path` (scatter), then a single collect step merges
them with :func:`collect_task_csvs` (gather).  No filesystem locking is
required or used.
"""

import glob
import os
import pandas as pd


def flatten_columns(df):
    """Flatten tuple columns (e.g. from bias fits) into 'var_stat' strings.

    Columns that are already plain strings are left unchanged.

    Returns:
        pd.DataFrame with string-only columns (copy).
    """
    df = df.copy()
    new_cols = []
    for c in df.columns:
        if isinstance(c, tuple):
            new_cols.append(f"{c[0]}_{c[1]}")
        else:
            new_cols.append(str(c))
    df.columns = new_cols
    return df


def task_csv_path(base_path, task_id):
    """Return the per-task CSV path derived from *base_path*.

    E.g. ``toys_batch_toys.csv`` → ``toys_batch_toys_task42.csv``.

    Args:
        base_path (str): Configured base CSV path.
        task_id (int): SLURM array task ID (or any unique integer).

    Returns:
        str: Per-task path, written exclusively by one process.
    """
    root, ext = os.path.splitext(base_path)
    return f"{root}_task{task_id}{ext}"


def append_csv(df, path):
    """Append *df* rows to a CSV file (serial use only).

    Safe for serial appends (interactive runs).  Do NOT
    call this from concurrent SLURM tasks — use :func:`task_csv_path`
    and :func:`collect_task_csvs` instead.

    Args:
        df (pd.DataFrame): Rows to append.
        path (str): Path to the output CSV file.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    write_header = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a") as f:
        df.to_csv(f, header=write_header, index=False)


def collect_task_csvs(base_path, output_path=None):
    """Merge all per-task CSVs from a scatter run into one file.

    Globs ``{stem}_task*.csv`` next to *base_path*, concatenates them,
    and writes the result to *output_path* (defaults to *base_path*).

    Args:
        base_path (str): Configured base CSV path (e.g. ``toys_batch_toys.csv``).
        output_path (str | None): Destination for merged output.  Defaults to
            *base_path*.

    Returns:
        pd.DataFrame: Merged data.

    Raises:
        FileNotFoundError: If no per-task CSV files are found.
    """
    root, ext = os.path.splitext(base_path)
    pattern = f"{root}_task*{ext}"
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No per-task CSVs found matching: {pattern}")
    df = pd.concat((pd.read_csv(p) for p in paths), ignore_index=True)
    dest = output_path or base_path
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    df.to_csv(dest, index=False)
    return df


def read_csv(path):
    """Read a CSV written by :func:`append_csv`.

    Returns:
        pd.DataFrame
    """
    return pd.read_csv(path)
