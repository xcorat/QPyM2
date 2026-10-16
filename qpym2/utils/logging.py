""" This file is intended to provide easy capabilities

Debug levels are named:

    5: Debug
    0: No log output

"""
import time

_DEBUG_ = 5

# TODO: Import debug level from main configurtion file
DEBUG = _DEBUG_

# Global diagnostics switches
TIMING_ENABLED = True
RUNTIME_TELEMETRY_ENABLED = True
RUNTIME_TELEMETRY_INTERVAL_S = 10

def log(level:int, *args, **kwargs):
    """ log output (print) if log `level` is less than the `DEBUG` level 
    
    level < 0 is skipped. Use this for specific debugging while keeping print statements

    Parameters
    ----------
    level : int
        The log level.
    *args and **kwargs : passed to print
    """ 
    if level < 0:
        # negative levels are skipped.
        return
    if DEBUG >= level:
        print(*args, **kwargs)

def debug(*args, **kwargs):
    log(5, *args, **kwargs)


def start_timer():
    """Return a monotonic timer start value."""
    return time.perf_counter()


def format_elapsed(seconds):
    """Format elapsed seconds for log readability."""
    return f"{seconds:.2f}s"


def log_elapsed(label, start_time, level=0):
    """Log elapsed time from start_time, if timing is enabled."""
    elapsed = time.perf_counter() - start_time
    if TIMING_ENABLED:
        log(level, f"timing: {label}: {format_elapsed(elapsed)}")
    return elapsed