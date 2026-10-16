"""
Implement any utility functions related to computing hardware.

TODO: rename this to jobs/cluster or something
"""
import os
import time


def get_ncores(print_state=False):
    """ Check the number of cores available on the node. Copied from NERSC doc.

    Returns:
        int: number of cores available on the node.
    """
    import multiprocessing as mp
    import os
    import psutil

    nthreads = psutil.cpu_count(logical=True)
    ncores = psutil.cpu_count(logical=False)
    nthreads_per_core = nthreads // ncores
    nthreads_available = len(os.sched_getaffinity(0))
    ncores_available = nthreads_available // nthreads_per_core

    assert nthreads == os.cpu_count()
    assert nthreads == mp.cpu_count()

    if print_state:
        print(f'{nthreads=}')
        print(f'{ncores=}')
        print(f'{nthreads_per_core=}')
        print(f'{nthreads_available=}')
        print(f'{ncores_available=}')

    return ncores_available


def get_runtime_usage(include_children=False):
    """Return runtime resource usage snapshot for debugging.

    Args:
        include_children (bool): include child process RSS in the total.

    Returns:
        dict: usage fields suitable for logging.
    """
    try:
        import psutil
    except Exception:
        return {
            'timestamp': time.time(),
            'pid': os.getpid(),
            'psutil_available': False,
        }

    proc = psutil.Process(os.getpid())
    mem_info = proc.memory_info()
    rss_bytes = mem_info.rss

    if include_children:
        for child in proc.children(recursive=True):
            try:
                rss_bytes += child.memory_info().rss
            except Exception:
                continue

    load = os.getloadavg() if hasattr(os, "getloadavg") else (None, None, None)
    cpu_pct = psutil.cpu_percent(interval=None)

    return {
        'timestamp': time.time(),
        'pid': os.getpid(),
        'rss_mb': rss_bytes / (1024 ** 2),
        'vms_mb': mem_info.vms / (1024 ** 2),
        'cpu_percent': cpu_pct,
        'threads': proc.num_threads(),
        'load1': load[0],
        'load5': load[1],
        'load15': load[2],
        'psutil_available': True,
    }
