"""Run the bin-size × Co-60-shift fit matrix for the muon analysis.

Inputs are generated from the configured staging ROOT files and written beside
each fit under ``infrastructure.output_paths.matrix_fit_dir``. By default,
the matrix and priors follow ``hyperparams/tech_note.yaml``.

Example::

    python -m qpym2.run.run_bin_size_fits --dry-run
    python -m qpym2.run.run_bin_size_fits --chains 4 --draws 5000
"""

import argparse
import copy
import math
from pathlib import Path

from qpym2.cfg import get_comp_cfg, get_h2uv_model
from qpym2.cfg.config_builders import build_input_config, get_jags_path
from qpym2.cfg.hydra_adapter import get_config
from qpym2.fit import fit as fit_module
from qpym2.fit.create_m2uv import create_input
from qpym2.fit.params import make_priors, model_cfg
from qpym2.io.bkg_model import read_mkchain
from qpym2.io.fit import read_input, write_fit
from qpym2.utils.TempConfig import TempConfig


def _value_tag(value):
    """Format a numeric matrix coordinate safely for directory names."""
    return f"{value:g}".replace(".", "p").replace("-", "m")


def _make_job_name(bin_size, co_shift):
    return f"bin_{_value_tag(bin_size)}_cs_{_value_tag(co_shift)}"


def _get_matrix_values(cfg):
    """Read and validate matrix axes from the selected hyperparameter config."""
    hyperparams = cfg.hyperparams
    shifts = [float(value) for value in hyperparams.co_shift_values]
    bin_sizes = [float(value) for value in hyperparams.bin_size_values]
    if not shifts or not bin_sizes:
        raise ValueError("co_shift_values and bin_size_values must both be non-empty")
    if any(value <= 0 for value in bin_sizes):
        raise ValueError("All bin sizes must be positive")
    if 1.0 not in bin_sizes:
        raise ValueError("bin_size_values must include 1 keV for nominal priors")
    return shifts, bin_sizes


def run_matrix(chains=4, draws=5000, mkchain_samples=None, dry_run=False):
    """Build inputs and run the configured bin-size × Co-60-shift matrix."""
    cfg = get_config(
        environment="m2v5_mu",
        analysis="tech_note", # Use tech_note_mu for seperate mu compoenent fit, standard otherwise.
        hyperparams="tech_note",
    )
    shifts, bin_sizes = _get_matrix_values(cfg)
    output_root = Path(cfg.infrastructure.output_paths.matrix_fit_dir).expanduser()
    if not output_root.is_absolute():
        raise ValueError(f"matrix_fit_dir must be absolute: {output_root}")

    co_cut = cfg.hyperparams.co_cuts
    ndbd_cutoff = cfg.hyperparams.ndbd_cutoff
    # Process the nominal 1-keV point first for each shift so its templates
    # can supply the prior means/widths reused at other bin sizes.
    ordered_bin_sizes = sorted(bin_sizes, key=lambda value: value != 1.0)
    jobs = [(shift, bin_size, output_root / _make_job_name(bin_size, shift))
            for shift in shifts for bin_size in ordered_bin_sizes]

    print("Config: m2v5_mu / tech_note_mu / tech_note")
    print(f"Cuts: co_cut={co_cut}, ndbd_cutoff={ndbd_cutoff}")
    print(f"Output: {output_root}")
    print(f"Matrix: {len(shifts)} co-shifts × {len(bin_sizes)} bin sizes = {len(jobs)} fits")
    for shift, bin_size, job_dir in jobs:
        print(f"  co_shift={shift:g} keV, bin_size={bin_size:g} keV -> {job_dir}")
    if dry_run:
        return

    output_root.mkdir(parents=True, exist_ok=True)
    input_cfg = build_input_config(cfg)
    jags_path = get_jags_path(cfg)
    base_hm = get_h2uv_model(cfg)
    comp_cfg_base = get_comp_cfg(cfg)
    co_cuts_dict = dict(cfg.cuts.co_cuts)
    ndbd_cuts_dict = dict(cfg.cuts.ndbd_cutoff)

    mkchain = read_mkchain(jags_path)
    sample_count = mkchain_samples
    if sample_count is None:
        sample_count = int(cfg.hyperparams.sample_mkchain)
    sample_count = min(sample_count, len(mkchain))
    if sample_count < 1:
        raise ValueError("At least one JAGS chain sample is required for priors")
    mksamples = mkchain.sample(n=sample_count, random_state=42).T

    priors_by_shift = {}
    input_write_cfg = {"outdir": str(output_root), "overwrite": True}
    total = len(jobs)

    for index, (co_shift, bin_size, job_dir) in enumerate(jobs, start=1):
        print(f"\n[{index}/{total}] co_shift={co_shift:g}, bin_size={bin_size:g}")

        hm = TempConfig(**{key: getattr(base_hm, key) for key in vars(base_hm)})
        hm.binsize = bin_size / math.sqrt(2)
        fit_cfg = model_cfg(
            hm,
            co_cut=co_cut,
            ndbd_cutoff=ndbd_cutoff,
            co_shift=co_shift,
            co_cuts_dict=co_cuts_dict,
            ndbd_cuts_dict=ndbd_cuts_dict,
            common_cuts=[],
        )
        # Each matrix point has its own directory, so standard file names are
        # sufficient and cannot collide across points.
        fit_cfg["shift_suffix"] = ""

        create_input(
            fit_cfg,
            input_cfg,
            jags_path,
            fit_name=job_dir.name,
            write_cfg=input_write_cfg,
        )
        input_path = job_dir / "input.h5"
        mctable, signal_df, data, edges = read_input(str(input_path))

        if bin_size == 1.0:
            priors_by_shift[co_shift] = make_priors(
                copy.deepcopy(comp_cfg_base),
                mctable,
                mksamples,
                sigma_scale=10.0,
                use_sample_std=False,
            )
        if co_shift not in priors_by_shift:
            raise RuntimeError(
                f"The 1-keV bin size is required to initialize priors for shift {co_shift:g}"
            )

        model, fit_comps = fit_module.setup_fit(
            mctable,
            data,
            signal_df,
            hm,
            smooth=hm.smooth,
            comp_cfg=copy.deepcopy(priors_by_shift[co_shift]),
            edges=edges,
        )
        trace = fit_module.run_fit(
            model,
            single_core=True,
            chains=chains,
            draws=draws,
        )
        write_fit(fit_comps, trace, outdir=str(job_dir), fname="fit", data=data)
        print(f"Saved input and fit under {job_dir}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the configured muon-analysis bin-size × Co-60-shift fit matrix."
    )
    parser.add_argument("--chains", type=int, default=4,
                        help="MCMC chains per fit (default: 4)")
    parser.add_argument("--draws", type=int, default=5000,
                        help="MCMC draws per chain (default: 5000)")
    parser.add_argument("--mkchain-samples", type=int, default=None,
                        help="JAGS samples for priors (default: hyperparameter config)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print configured matrix and destinations without running")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.chains < 1 or args.draws < 1:
        raise SystemExit("--chains and --draws must be positive")
    if args.mkchain_samples is not None and args.mkchain_samples < 1:
        raise SystemExit("--mkchain-samples must be positive")
    run_matrix(
        chains=args.chains,
        draws=args.draws,
        mkchain_samples=args.mkchain_samples,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()