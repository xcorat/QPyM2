"""Generic multiprocessing pool runner.

Two entry points:

* :func:`run_pool` — homogeneous batch: apply one function to a list of args.
* :class:`JobPool` — heterogeneous queue: ``add`` / ``extend`` arbitrary
  ``(fn, args, kwargs)`` tuples, then ``run`` them under one pool.

Pickling constraint: every ``fn`` (and every value in ``args`` / ``kwargs``)
must be importable / picklable, since workers are spawned in separate
processes.

``ncores`` resolution: ``None`` → :func:`os.cpu_count`; ``1`` (or a single
job) → run sequentially in-process without spawning a pool.
"""

from __future__ import annotations

import os
from multiprocessing import Pool
from typing import Any, Callable, Iterable, List, Optional, Tuple

__all__ = ['run_pool', 'JobPool']


def _invoke(task: Tuple[Callable[..., Any], tuple, dict]) -> Any:
    fn, args, kwargs = task
    return fn(*args, **kwargs)


def _resolve_ncores(ncores: Optional[int], njobs: int) -> int:
    if ncores is None or ncores < 1:
        return os.cpu_count() or 1
    return min(ncores, njobs) if njobs > 0 else 1


def _run_tasks(tasks: List[Tuple[Callable[..., Any], tuple, dict]],
               ncores: Optional[int]) -> list:
    if not tasks:
        return []
    n = _resolve_ncores(ncores, len(tasks))
    if n == 1 or len(tasks) == 1:
        return [_invoke(t) for t in tasks]
    with Pool(processes=n) as pool:
        # imap_unordered drains as workers complete; results come back in
        # completion order, not submission order.
        return list(pool.imap_unordered(_invoke, tasks))


def run_pool(fn: Callable[[Any], Any],
             args_list: Iterable[Any],
             ncores: Optional[int] = None) -> list:
    """Apply ``fn`` to each item in ``args_list`` in a pool.

    ``fn`` is called with a *single* positional argument (the item). To pass
    multiple arguments per call, package them into a tuple and have ``fn``
    unpack — or use :class:`JobPool` instead.

    Returns results in completion order (or submission order when running
    sequentially).
    """
    tasks = [(fn, (item,), {}) for item in args_list]
    return _run_tasks(tasks, ncores)


class JobPool:
    """Accumulate heterogeneous tasks and run them under one pool.

    Example::

        pool = JobPool()
        pool.add(stage_one_file, 'a.root', cfg)
        pool.extend(stage_one_file, [(f, cfg) for f in many_files])
        pool.add(stage_data, data_fname, cfg, blinded=True)
        results = pool.run(ncores=8)
    """

    def __init__(self) -> None:
        self._tasks: List[Tuple[Callable[..., Any], tuple, dict]] = []

    def add(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        """Queue a single task: ``fn(*args, **kwargs)``."""
        self._tasks.append((fn, args, kwargs))

    def extend(self, fn: Callable[..., Any],
               args_iter: Iterable[Any]) -> None:
        """Queue many tasks all using the same ``fn``.

        Each item in ``args_iter`` is either a tuple (passed as positional
        args) or a single non-tuple value (passed as one positional arg).
        """
        for item in args_iter:
            if isinstance(item, tuple):
                self._tasks.append((fn, item, {}))
            else:
                self._tasks.append((fn, (item,), {}))

    def __len__(self) -> int:
        return len(self._tasks)

    def run(self, ncores: Optional[int] = None) -> list:
        """Execute all queued tasks; the queue is *not* cleared on return."""
        return _run_tasks(self._tasks, ncores)
