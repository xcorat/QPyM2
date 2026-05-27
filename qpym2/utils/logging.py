"""Logging setup for QPyM2 entry points.

Usage in any ``main()``::

    from qpym2.utils.logging import setup_logging
    cfg = get_config(...)
    setup_logging(cfg)

Reads ``cfg.logging.{level,format,log_dirpath,quiet_others}``.

When ``log_dirpath`` is non-empty:
  - Python log records are tee'd to ``<log_dirpath>/<caller_name>_<YYYYMMDD_HHMMSS>.log``.
  - C++ stderr (ROOT ``Error in <...>`` messages) is also tee'd to the
    same file via a low-level fd-2 pipe redirect + background thread.

When ``quiet_others`` is true (default), the root logger stays at WARNING so
third-party libraries (numba, matplotlib, arviz, …) stay silent; only the
``qpym2`` package logger is set to ``level``.
"""

from __future__ import annotations

import datetime
import inspect
import logging
import os
import pathlib
import threading
from typing import Optional


def _start_stderr_tee(logfile_path: pathlib.Path) -> None:
    """Redirect C++ fd-2 stderr to both console and *logfile_path*.

    Replaces fd 2 with the write end of a pipe; a daemon thread reads from
    the read end and copies bytes to both the original stderr fd and the log
    file.  On process exit the write end closes automatically, signalling EOF
    and letting the thread drain cleanly.
    """
    orig_fd = os.dup(2)          # save real stderr fd
    r_fd, w_fd = os.pipe()
    os.dup2(w_fd, 2)             # fd 2 now writes into the pipe
    os.close(w_fd)

    def _tee() -> None:
        orig = os.fdopen(orig_fd, 'wb', buffering=0)
        log = open(logfile_path, 'ab', buffering=0)
        try:
            with os.fdopen(r_fd, 'rb', buffering=0) as pipe_r:
                while True:
                    data = pipe_r.read(4096)
                    if not data:
                        break
                    orig.write(data)
                    orig.flush()
                    log.write(data)
                    log.flush()
        finally:
            orig.close()
            log.close()

    threading.Thread(target=_tee, daemon=True, name='stderr-tee').start()


def setup_logging(cfg, script_name: Optional[str] = None) -> None:
    """Configure logging from *cfg.logging*.

    Args:
        cfg: Hydra DictConfig with a ``logging`` sub-key containing
             ``level``, ``format``, ``log_dirpath``, and ``quiet_others``.
        script_name: Stem used in the log filename.  Defaults to the
                     filename stem of the calling module.
    """
    log_cfg = cfg.logging
    fmt = logging.Formatter(log_cfg.format)

    # Root logger: WARNING when quiet_others=true (suppresses third-party noise)
    root = logging.getLogger()
    root_level = logging.WARNING if log_cfg.get('quiet_others', True) else log_cfg.level
    root.setLevel(root_level)

    # qpym2 package logger honours the configured level
    logging.getLogger('qpym2').setLevel(log_cfg.level)

    # Console handler attached to root (catches root_level+ and qpym2 propagation)
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    root.addHandler(sh)

    # File handler + C++ stderr tee (optional)
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

        # Tee C++ stderr (ROOT Error/Warning messages) into the same file
        _start_stderr_tee(logfile)

        logging.getLogger('qpym2').debug('Logging to %s', logfile)
