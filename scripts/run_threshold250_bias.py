"""
Signal-injection (bias) fits for a BM threshold systematic variant.

All settings are read from settings.yaml inside the syst directory.
Output toys_bias.csv is written into the same directory.

Usage:
    python scripts/run_threshold250_bias.py
    python scripts/run_threshold250_bias.py --syst-dir path/to/Threshold250
    python scripts/run_threshold250_bias.py --n 200 --bias-min 5 --bias-max 200 --force
"""

import argparse
import multiprocessing
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))

from qpym2.cfg import get_config
from qpym2.cfg.config_builders import get_h2uv_model
from qpym2.fit.params import model_cfg
from qpym2.utils.hists import get_empty_hist, get_sum_hist, sample_hist
from qpym2.utils.mathy import get_mode
from qpym2.utils.pool_monitor import run_pool_map_with_monitor
from omegaconf import OmegaConf
# NOTE: qpym2.io.fit (arviz) and qpym2.fit.posterior_checks (pymc) are imported
# lazily inside functions to keep the parent process free of threading locks
# before the spawn Pool is created.

# ---------------------------------------------------------------------------
# Spawn-safe worker infrastructure (mirrors run_syst_lss.py pattern)
# ---------------------------------------------------------------------------
_WORKER_SHARED = {}


def _worker_pool_init(shared):
    """Spawn worker initializer: receive shared data via pickle.

    spawn gives each worker a clean interpreter with no inherited threading
    locks, so arviz/scipy/matplotlib never cause fork-deadlocks in pm.sample.
    """
    global _WORKER_SHARED
    _WORKER_SHARED = shared


def _bias_worker(job):
    """Run one bias fit in a spawn worker."""
    import os, time
    _pid = os.getpid()
    _t0 = time.perf_counter()
    print(f"  [worker pid={_pid}] START", flush=True)

    # Set per-PID compiledir BEFORE first pytensor/pymc import.
    _worker_dir = os.path.join(
        os.environ.get('TMPDIR', '/tmp'), f'.pytensor_worker_{_pid}'
    )
    os.makedirs(_worker_dir, exist_ok=True)
    os.environ['PYTENSOR_FLAGS'] = f'compiledir={_worker_dir}'

    import pytensor  # noqa: F401 — triggers compiledir setup in fresh process
    print(f"  [worker pid={_pid}] pytensor t={time.perf_counter()-_t0:.1f}s", flush=True)

    import pathlib
    (pathlib.Path.home() / '.cache' / 'arviz').mkdir(parents=True, exist_ok=True)

    i, data_i, sig_norm = job
    comps  = _WORKER_SHARED['comps']
    xedges = _WORKER_SHARED['xedges']
    yedges = _WORKER_SHARED['yedges']

    from qpym2.fit.posterior_checks import fit_m2  # lazy — pymc first loaded here
    trace = fit_m2(data_i, comps, xedges, yedges)
    print(f"  [worker pid={_pid}] DONE t={time.perf_counter()-_t0:.1f}s", flush=True)
    return (i, trace, float(sig_norm))

DEFAULT_SYST_DIR = (
    ROOT_DIR / "data_2syst" / "syst_results" / "bm" / "thresholds" / "Threshold250"
)


def main():
    parser = argparse.ArgumentParser(description="Bias fits for a BM threshold systematic")
    parser.add_argument("--syst-dir",  default=str(DEFAULT_SYST_DIR),
                        help="Path to syst variant dir containing settings.yaml and fit files")
    parser.add_argument("--n",         type=int,   default=500)
    parser.add_argument("--bias-min",  type=float, default=5)
    parser.add_argument("--bias-max",  type=float, default=200)
    parser.add_argument("--ncores",    type=int,   default=25)
    parser.add_argument("--force",     action="store_true",
                        help="Re-run even if toys_bias.csv already exists")
    args = parser.parse_args()

    syst_dir = Path(args.syst_dir)
    bias_csv = syst_dir / "toys_bias.csv"

    if bias_csv.exists() and not args.force:
        print(f"toys_bias.csv already exists in {syst_dir} — skipping (use --force to re-run)")
        sys.exit(0)

    # --- Read all settings from the yaml in the syst dir ---
    settings_path = syst_dir / "settings.yaml"
    if not settings_path.exists():
        print(f"ERROR: settings.yaml not found in {syst_dir}")
        sys.exit(1)
    with open(settings_path) as f:
        s = yaml.safe_load(f)

    environment = s["environment"]
    analysis    = s["analysis"]
    model_key   = s["model_key"]
    co_cut      = s["co_cut"]
    ndbd_cutoff = s["ndbd_cutoff"]
    co_shift    = float(s["co_shift"])

    print(f"syst_dir    : {syst_dir}")
    print(f"environment : {environment}  analysis : {analysis}  model_key : {model_key}")
    print(f"co_cut={co_cut}  ndbd_cutoff={ndbd_cutoff}  co_shift={co_shift}")
    print(f"bias range  : [{args.bias_min}, {args.bias_max}]  n={args.n}  ncores={args.ncores}")

    # --- Build mcfg from config + settings ---
    cfg  = get_config(environment=environment, analysis=analysis)
    hm   = get_h2uv_model(cfg, model_key)
    cc   = OmegaConf.to_container(cfg.cuts.co_cuts,     resolve=True)
    nc   = OmegaConf.to_container(cfg.cuts.ndbd_cutoff, resolve=True)
    mcfg = model_cfg(hm, co_cut=co_cut, ndbd_cutoff=ndbd_cutoff, co_shift=co_shift,
                     co_cuts_dict=cc, ndbd_cuts_dict=nc)

    fit_fname = f"fit{mcfg['shift_suffix']}"   # e.g. fit_ushiftp10
    print(f"fit files   : {fit_fname}.{{nc,h5}}")

    # --- Load the existing fit (lazy import of arviz-heavy qpym2.io.fit) ---
    from qpym2.io.fit import read_fit
    fit_comps, idata, _ = read_fit(str(syst_dir), fname=fit_fname)
    print(f"components  : {list(fit_comps.index)}")

    # --- Pre-sample toy datasets in parent (pure numpy — no pymc/arviz Pool) ---
    _, xedges, yedges = get_empty_hist(mcfg['h2uv_model'])
    vars_list = list(fit_comps.index)
    stacked = idata.posterior[vars_list].stack(draws=("chain", "draw"))
    point = [get_mode(stacked[v].values) for v in vars_list]
    points = np.array([point] * args.n)
    isig = vars_list.index('ndbd')
    points[:, isig] = np.random.uniform(args.bias_min, args.bias_max, size=args.n)
    signal_norms = points[:, isig].copy()

    print(f"Sampling {args.n} toy datasets ...", flush=True)
    data_list = [sample_hist(get_sum_hist(fit_comps, p), xedges, yedges) for p in points]
    print(f"Toy datasets ready. Running {args.n} bias fits (ncores={min(args.ncores, args.n)}, spawn) ...", flush=True)

    # --- Spawn Pool (no inherited threading locks) ---
    global _WORKER_SHARED
    _WORKER_SHARED = {'comps': fit_comps, 'xedges': xedges, 'yedges': yedges}
    jobs = [(i, data_list[i], float(signal_norms[i])) for i in range(args.n)]
    ncores = min(args.ncores, args.n)

    mp_ctx = multiprocessing.get_context('spawn')
    with mp_ctx.Pool(processes=ncores, maxtasksperchild=1,
                     initializer=_worker_pool_init,
                     initargs=(_WORKER_SHARED,)) as pool:
        results = run_pool_map_with_monitor(pool, 'bias_fits', _bias_worker, jobs)

    rows = []
    for i, t, sig_norm in results:
        ndbd_post = t.posterior["ndbd"].values.flatten()
        rows.append({
            "ndbd_mode":   float(get_mode(ndbd_post)),
            "ndbd_p90":    float(np.percentile(ndbd_post, 90)),
            "ndbd_inj_":   float(sig_norm),
            "co_cut":      co_cut,
            "ndbd_cutoff": ndbd_cutoff,
            "co_shift":    co_shift,
        })

    df_bias = pd.DataFrame(rows)
    df_bias.to_csv(str(bias_csv), index=False)
    print(f"Saved {len(df_bias)} rows → {bias_csv}")


if __name__ == "__main__":
    main()
