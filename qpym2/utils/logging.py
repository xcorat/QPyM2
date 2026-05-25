"""Logging setup for QPyM2 entry points.

Usage in any ``main()``::

    from qpym2.utils.logging import setup_logging
    cfg = get_config(...)
    setup_logging(cfg)

Reads ``cfg.logging.{level,format,log_dirpath}``.  When ``log_dirpath`` is
non-empty the root logger is tee'd to
``<log_dirpath>/<caller_name>_<YYYYMMDD_HHMMSS>.log``.
"""

from __future__ import annotations

import datetime
import inspect
import logging
import pathlib
from typing import Optional


def setup_logging(cfg, script_name: Optional[str] = None) -> None:
    """Configure the root logger from *cfg.logging*.

    Args:
        cfg: Hydra DictConfig with a ``logging`` sub-key containing
             ``level``, ``format``, and ``log_dirpath``.
        script_name: Stem used in the log filename.  Defaults to the
                     ``__name__`` of the calling module.
    """
    log_cfg = cfg.logging
    fmt = logging.Formatter(log_cfg.format)

    root = logging.getLogger()
    root.setLevel(log_cfg.level)

    # Console handler (always)
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    root.addHandler(sh)

    # File handler (optional)
    if log_cfg.log_dirpath:
        if script_name is None:
            frame = inspect.stack()[1]
            script_name = pathlib.Path(frame.filename).stem

        logdir = pathlib.Path(log_cfg.log_dirpath)
        logdir.mkdir(parents=True, exist_ok=True)
        ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        logfile = logdir / f'{script_name}_{ts}.log'

        fh = logging.FileHandler(logfile)
        fh.setFormatter(fmt)
        root.addHandler(fh)
        logging.debug('Logging to %s', logfile)
