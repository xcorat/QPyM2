"""Run signal-injection bias tests for existing standard and mu bin matrices.

Each matrix-point directory containing ``input.h5`` and ``fit.nc`` receives a
``toys_bias.csv`` with per-toy posterior summaries and the injected signal.
Run with ``python -m qpym2.run.run_bin_size_bias_fits --dry-run`` to inspect
which matrix points will be processed. Packed SLURM jobs run one matrix point
per GNU Parallel task, sequentially fitting toys with PyMC parallelizing chains.
"""

import argparse
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from qpym2.utils.hists import get_sum_hist, sample_hist
from qpym2.utils.mathy import get_mode

_WORKER_SHARED = {}
_MATRIX_DIR_RE = re.compile(r"^(?P<family>mu_)?bin_(?P<bin>.+)_cs_(?P<shift>.+)$")


def _worker_pool_init(shared):
    """Install per-matrix-point data for sequential toy execution."""
    global _WORKER_SHARED
    _WORKER_SHARED = shared


def _disable_pymc_chain_threadpool_limits():
    """Replace PyMC's per-chain ``threadpool_limits`` with a no-op.

    ``pymc.sampling.ParallelSampler`` wraps every sampling chain in
    ``threadpoolctl.threadpool_limits(limits=blas_cores)``. Constructing that
    context manager scans all loaded shared libraries
    (``dl_iterate_phdr`` + ``dlopen``) even when ``limits`` is ``None``. Run
    inside a child forked from a process that already initialized BLAS/OpenMP,
    that scan intermittently aborts the child with
    ``Fatal Python error: Aborted`` in ``threadpoolctl._load_libraries``; the
    parent then dies with ``ConnectionResetError: [Errno 104]`` in
    ``recv_draw``.

    This runner already pins BLAS/OpenMP to one thread per chain through the
    ``*_NUM_THREADS`` environment variables, so PyMC's limiting is redundant.
    Disabling it removes the aborting library scan while keeping fork-based
    chains (``spawn`` avoids the abort too but is ~3x slower).
    """
    import contextlib

    import pymc.sampling.parallel as _parallel

    def _noop(**kwargs):
        return contextlib.nullcontext()

    _parallel.threadpool_limits = _noop


def _bias_worker(job):
    """Fit one injected toy dataset; PyMC parallelizes its chains."""
    pid = os.getpid()
    worker_dir = os.path.join(os.environ.get("TMPDIR", "/tmp"), f".pytensor_worker_{pid}")
    os.makedirs(worker_dir, exist_ok=True)
    os.environ["PYTENSOR_FLAGS"] = f"compiledir={worker_dir}"
    import pytensor  # noqa: F401 -- initialize compiledir before importing PyMC

    index, data, signal_norm, random_seed = job
    # Keep this in the main process rather than a Pool worker: PyMC must be
    # able to create its own chain processes. GNU Parallel controls how many
    # matrix-point commands run concurrently on the node.
    import pymc as pm
    from qpym2.models import m2

    # PyMC forks one process per chain; its default per-chain threadpool_limits
    # scan aborts in forked children in this environment (see helper docstring).
    _disable_pymc_chain_threadpool_limits()

    comps = _WORKER_SHARED["comps"]
    fik = m2.make_fik(
        comps, data, _WORKER_SHARED["xedges"], _WORKER_SHARED["yedges"]
    )
    model = m2.m2_model(comps, fik)
    with model:
        trace = pm.sample(
            draws=_WORKER_SHARED["draws"],
            tune=_WORKER_SHARED["tune"],
            chains=_WORKER_SHARED["chains"],
            cores=_WORKER_SHARED["chains"],
            return_inferencedata=True,
            target_accept=0.9,
            progressbar=False,
            idata_kwargs={"log_likelihood": True},
            random_seed=random_seed,
        )
    return index, trace, signal_norm


def _default_root(environment, analysis, hyperparams):
    from qpym2.cfg.hydra_adapter import get_config

    cfg = get_config(environment=environment, analysis=analysis, hyperparams=hyperparams)
    root = cfg.infrastructure.output_paths.matrix_fit_dir
    return Path(root).expanduser()


def _matrix_points(root, family="all"):
    """Find complete matrix fit directories, optionally filtering by family."""
    points = []
    for fit_dir in sorted(root.iterdir()):
        if not fit_dir.is_dir():
            continue
        match = _MATRIX_DIR_RE.fullmatch(fit_dir.name)
        if match is None:
            continue
        label = "mu" if match.group("family") else "standard"
        if family != "all" and label != family:
            continue
        points.append((fit_dir, label, match.group("bin"), match.group("shift")))
    return points


def _parse_tag(tag):
    """Convert the numeric filename tag to a float."""
    return float(tag.replace("m", "-").replace("p", "."))


def _resolve_element(points, element):
    """Return exactly one matrix point by its directory name."""
    matches = [point for point in points if point[0].name == element]
    if not matches:
        available = ", ".join(point[0].name for point in points)
        raise ValueError(f"Unknown matrix element {element!r}. Available: {available}")
    return matches[0]


def _run_point(fit_dir, family, bin_size, co_shift, n, bias_min, bias_max,
               chains, draws, tune, seed, force):
    csv_path = fit_dir / "toys_bias.csv"
    lock_path = fit_dir / ".toys_bias.lock"
    try:
        lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError as exc:
        raise RuntimeError(f"Bias job already active for {fit_dir} (lock: {lock_path})") from exc
    os.close(lock_fd)

    try:
        if csv_path.exists() and not force:
            existing_count = len(pd.read_csv(csv_path))
            if existing_count == n:
                print(f"Skip {fit_dir}: found the requested {n} rows in {csv_path}")
                return
            raise FileExistsError(
                f"{csv_path} has {existing_count} rows, requested {n}; use --force to replace it"
            )

        missing = [name for name in ("input.h5", "fit.nc", "fit.h5")
                   if not (fit_dir / name).is_file()]
        if missing:
            raise FileNotFoundError(f"{fit_dir} is missing required fit artifacts: {', '.join(missing)}")
        if bias_max <= bias_min:
            raise ValueError("--bias-max must be greater than --bias-min")

        # Keep arviz/PyMC imports out of the module top level: spawn workers
        # import this module before setting their isolated PyTensor cache.
        from qpym2.io.fit import read_fit
        from qpym2.io.hdf import read_hist_edges

        comps, trace, _ = read_fit(str(fit_dir), fname="fit")
        # Only the bin edges are needed from input.h5: the per-component
        # pdf_hist lives in fit.h5 and the toy data is regenerated per fit, so
        # reading the mctable/signal/data would be wasted I/O.
        xedges, yedges = read_hist_edges(str(fit_dir / "input.h5"))
        component_names = list(comps.index)
        if "ndbd" not in component_names:
            raise ValueError(f"No 'ndbd' component in {fit_dir}; found {component_names}")

        stacked = trace.posterior[component_names].stack(draws=("chain", "draw"))
        mode_point = np.asarray([get_mode(stacked[name].values) for name in component_names])
        rng = np.random.default_rng(seed)
        injections = rng.uniform(bias_min, bias_max, size=n)
        points = np.repeat(mode_point[None, :], n, axis=0)
        points[:, component_names.index("ndbd")] = injections
        np.random.seed(seed)

        shared = {
            "comps": comps,
            "xedges": xedges,
            "yedges": yedges,
            "chains": chains,
            "draws": draws,
            "tune": tune,
        }
        print(
            f"{fit_dir.name} [{family}, bin={bin_size}, shift={co_shift}]: "
            f"running {n} bias fits, {chains} chains per fit, "
            "one active toy fit in this matrix-point task",
            flush=True,
        )

        _worker_pool_init(shared)
        # The command executes under GNU Parallel, so permit PyMC's chain
        # subprocesses while limiting BLAS/OpenMP to one thread per chain.
        for variable in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            os.environ[variable] = "1"
        from qpym2.fit.results import _get_stats

        rows = []
        for completed, (norms, injection) in enumerate(zip(points, injections), start=1):
            dataset = sample_hist(get_sum_hist(comps, norms), xedges, yedges)
            job = (completed - 1, dataset, float(injection), int(seed + completed))
            index, fitted_trace, signal_norm = _bias_worker(job)
            row = {
                "toy_id": index,
                "logp_mode": float(get_mode(fitted_trace.log_likelihood.pik.values.flatten())),
                "ndbd_inj_": float(signal_norm),
                "model_": fit_dir.name,
                "type_": "bias",
                "family": family,
                "bin_size_keV": _parse_tag(bin_size),
                "co_shift_keV": _parse_tag(co_shift),
            }
            for variable in component_names:
                for stat, value in _get_stats(fitted_trace.posterior[variable].values.flatten()).items():
                    row[f"{variable}_{stat}"] = float(value)
            rows.append(row)
            del fitted_trace
            if completed % 10 == 0 or completed == n:
                print(f"Completed {completed}/{n} bias toys", flush=True)

        if len(rows) != n:
            raise RuntimeError(f"Expected {n} bias results, received {len(rows)}")

        result_df = pd.DataFrame(rows).sort_values("toy_id", ignore_index=True)
        temp_csv = fit_dir / f".toys_bias.{os.getpid()}.csv.tmp"
        result_df.to_csv(temp_csv, index=False)
        temp_csv.replace(csv_path)
        print(f"Saved {len(result_df)} bias results: {csv_path}")
    finally:
        lock_path.unlink(missing_ok=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run signal-injection bias fits across the standard/mu bin matrices."
    )
    parser.add_argument("--root", type=Path, default=None,
                        help="Matrix root containing bin_* and mu_bin_* directories")
    parser.add_argument("--environment", default="m2v5_mu",
                        help="Hydra environment for the default matrix root")
    parser.add_argument("--analysis", default="tech_note",
                        help="Hydra analysis for the default matrix root")
    parser.add_argument("--hyperparams", default="tech_note",
                        help="Hydra hyperparameters for the default matrix root")
    parser.add_argument("--family", choices=("all", "standard", "mu"), default="all",
                        help="Which matrix family to process (default: all)")
    parser.add_argument("--element", default=None,
                        help="Run exactly one matrix directory, e.g. bin_1_cs_0p3 or mu_bin_5_cs_1")
    parser.add_argument("--n", type=int, default=200,
                        help="Bias fits per matrix point (default: 200)")
    parser.add_argument("--bias-min", type=float, default=5.0,
                        help="Minimum injected signal count (default: 5)")
    parser.add_argument("--bias-max", type=float, default=200.0,
                        help="Maximum injected signal count (default: 200)")
    parser.add_argument("--chains", type=int, default=4,
                        help="Chains per toy fit (each uses one core; default: 4)")
    parser.add_argument("--draws", type=int, default=4000,
                        help="MCMC draws per chain (default: 4000)")
    parser.add_argument("--tune", type=int, default=1000,
                        help="MCMC tune steps per chain (default: 1000)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Base random seed; each matrix point gets a distinct seed")
    parser.add_argument("--force", action="store_true",
                        help="Replace existing toys_bias.csv files")
    parser.add_argument("--dry-run", action="store_true",
                        help="List matching fit points without launching fits")
    parser.add_argument("--list-elements", action="store_true",
                        help="Print complete element directory names only (for scheduler array indexing)")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.n < 1 or args.chains < 1 or args.draws < 1 or args.tune < 0:
        raise SystemExit("--n, --chains, --draws must be positive and --tune non-negative")
    if args.bias_max <= args.bias_min:
        raise SystemExit("--bias-max must be greater than --bias-min")

    root = (args.root or _default_root(
        args.environment, args.analysis, args.hyperparams
    )).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Matrix root does not exist: {root}")
    points = _matrix_points(root, args.family)
    if not points:
        raise FileNotFoundError(f"No matching bin-matrix directories found under {root}")
    if args.element is not None:
        points = [_resolve_element(points, args.element)]
    if args.list_elements:
        for fit_dir, _, _, _ in points:
            if all((fit_dir / name).is_file() for name in ("input.h5", "fit.h5", "fit.nc")):
                print(fit_dir.name)
        return

    print(f"Matrix root: {root}")
    print(f"Selected fit points: {len(points)}; bias fits per point: {args.n}")
    for fit_dir, family, bin_size, co_shift in points:
        artifacts = all((fit_dir / name).is_file() for name in ("input.h5", "fit.h5", "fit.nc"))
        print(f"  {fit_dir.name} ({family}; artifacts={'ready' if artifacts else 'incomplete'})")
    if args.dry_run:
        return

    start = time.perf_counter()
    for point_index, (fit_dir, family, bin_size, co_shift) in enumerate(points):
        _run_point(
            fit_dir, family, bin_size, co_shift, args.n,
            args.bias_min, args.bias_max, args.chains, args.draws, args.tune,
            args.seed + point_index, args.force,
        )
    print(f"Completed bias matrix in {time.perf_counter() - start:.1f}s")


if __name__ == "__main__":
    main()
