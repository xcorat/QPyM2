"""Utilities for monitored multiprocessing pool execution."""

import time

from qpym2.utils import hw
from qpym2.utils import logging as qlog


def emit_runtime_telemetry(label, elapsed_seconds):
    """Print periodic runtime telemetry for long-running pool jobs."""
    if not qlog.RUNTIME_TELEMETRY_ENABLED:
        return

    usage = hw.get_runtime_usage(include_children=True)
    elapsed = qlog.format_elapsed(elapsed_seconds)

    if not usage.get('psutil_available', True):
        print(f"telemetry: {label}: elapsed={elapsed} psutil=unavailable pid={usage.get('pid')}")
        return

    print(
        f"telemetry: {label}: "
        f"elapsed={elapsed} "
        f"rss_mb={usage['rss_mb']:.1f} "
        f"vms_mb={usage['vms_mb']:.1f} "
        f"cpu_pct={usage['cpu_percent']:.1f} "
        f"threads={usage['threads']} "
        f"load1={usage['load1']:.2f}"
    )


def run_pool_map_with_monitor(pool, label, fn, items, timeout_s=60 * 80):
    """Run pool map asynchronously with timeout and optional telemetry."""
    start_time = qlog.start_timer()
    async_result = pool.map_async(fn, items)

    telemetry_interval = max(float(qlog.RUNTIME_TELEMETRY_INTERVAL_S), 1.0)
    next_telemetry = 0.0

    while not async_result.ready():
        elapsed = time.perf_counter() - start_time

        if elapsed >= next_telemetry:
            emit_runtime_telemetry(label, elapsed)
            next_telemetry += telemetry_interval

        if elapsed >= timeout_s:
            print(f"Timeout in run_sample_post ({label}) after {qlog.format_elapsed(elapsed)}")
            pool.terminate()
            pool.join()
            return []

        time.sleep(1)

    results = async_result.get()
    qlog.log_elapsed(f"{label} pool", start_time)
    return results
