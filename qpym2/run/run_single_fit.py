"""
Atomic PyMC Bayesian fit runner for the M2 analysis.

Runs a single fit for one parameter point (co_cut, ndbd_cutoff, co_shift, etc.).
Designed to be called standalone, from a SLURM array job, or imported as a library
by other runners (run_fits.py, run_toy_fits.py).

Usage:
    # Interactive single fit
    python -m qpym2.run.run_single_fit \\
        --co-cut v150 --ndbd-cutoff u1740 --co-shift 0.8

    # SLURM array (via gen_params.py fit + run_job.sh)
    python -m qpym2.run.run_single_fit \\
        --environment tech_note --co-cut gamma10_v150 \\
        --ndbd-cutoff u1740 --co-shift 1.0 --blinded false

    # Test prior generation
    python -m qpym2.run.run_single_fit --test-priors
"""

import argparse

from qpym2.cfg.hydra_adapter import get_config
from qpym2.cfg.config_builders import (
    get_jags_path, get_data_names, get_fit_dir, get_h2uv_model, get_comp_cfg,
    get_corrmat_path,
)
from qpym2.io import fit as fit_io
from qpym2.io.fit import read_input
from qpym2.io.bkg_model import read_mkchain, get_mcsyst_samples
from qpym2.fit.params import model_cfg, make_priors
from qpym2.fit import fit as fit
from qpym2.run.run_toy_fits import _build_model_grid


# ---------------------------------------------------------------------------
# Fit orchestration
# ---------------------------------------------------------------------------

def run_fit(input_fname, hm, chains=4, draws=5000,
            dataname='data', mksamples=None, comp_cfg=None, single_core=False,
            sigma_scale=1.0, use_sample_std=False):
    """Run a single PyMC MCMC fit.

    Reads the input HDF5, optionally computes priors from MC chain samples,
    builds the PyMC model, runs MCMC sampling, and returns the results.

    Args:
        input_fname: Path to input HDF5 file.
        hm: TempConfig histogram model (carries ROI geometry and smoothing).
        chains: Number of MCMC chains.
        draws: Number of MCMC draws per chain.
        dataname: HDF5 dataset key for data.
        mksamples: Optional MC chain samples for prior generation.
        comp_cfg: Component configuration list. If None, falls back to
            the hardcoded default in fit.fit.

    Returns:
        Tuple of (m2model, trace, fit_comps, data).
    """
    mctable, signal_df, data, edges = read_input(
        input_fname, dataname=dataname, signal_name='ndbd')

    if comp_cfg is None:
        from qpym2.fit.fit import comp_cfg

    if mksamples is not None:
        comp_cfg_new = make_priors(comp_cfg, mctable, mksamples,
                                   sigma_scale=sigma_scale,
                                   use_sample_std=use_sample_std)
    else:
        comp_cfg_new = comp_cfg

    m2model, fit_comps = fit.setup_fit(
        mctable, data, signal_df, hm, smooth=hm.smooth,
        comp_cfg=comp_cfg_new, edges=edges)

    trace = fit.run_fit(m2model, single_core=single_core, chains=chains, draws=draws)

    print(fit_comps['prior'])

    return m2model, trace, fit_comps, data


def setup_and_run(fitdir, model_cfg, data_name, data_blinded_name,
                  mksamples=None, write=False, comp_cfg=None,
                  input_name='input', input_noroi_name='input_noroi',
                  chains=8, draws=20000, single_core=False,
                  sigma_scale=1.0, use_sample_std=False):
    """Setup and run a single fit for one model configuration point.

    Resolves the input file path (standard vs no-ROI) and data dataset
    name (blinded vs unblinded), then delegates to ``run_fit``.

    Args:
        fitdir: Directory containing the fit input HDF5 and where output
            (trace .nc + components .h5) will be written.
        model_cfg: Model configuration dict (from ``fit.params.model_cfg``).
        data_name: HDF5 dataset key for unblinded data.
        data_blinded_name: HDF5 dataset key for blinded data.
        mksamples: Optional MC chain samples for prior generation.
        write: If True, write fit results to disk.
        comp_cfg: Component configuration list.
        input_name: Base filename (without .h5) for the standard input file.
        input_noroi_name: Base filename (without .h5) for the no-ROI input file.
        chains: Number of MCMC chains (default: 8).
        draws: Number of MCMC draws per chain (default: 20000).

    Returns:
        Tuple of (m2model, trace, fit_comps, data).
    """
    h2uv_model = model_cfg['h2uv_model']

    shift_suffix = model_cfg.get('shift_suffix', '')
    fname = input_noroi_name if model_cfg['skip_roi'] else input_name
    infn = f"{fitdir}/{fname}{shift_suffix}.h5"

    if model_cfg['blinded']:
        if isinstance(model_cfg['blinded'], str):
            dataname = model_cfg['blinded']
        else:
            dataname = data_blinded_name
    else:
        dataname = data_name

    m2model, trace, fit_comps, data = run_fit(
        infn, hm=h2uv_model, chains=chains, draws=draws,
        dataname=dataname, mksamples=mksamples, comp_cfg=comp_cfg,
        single_core=single_core, sigma_scale=sigma_scale,
        use_sample_std=use_sample_std)

    if write:
        fname = 'fit'
        if model_cfg['skip_roi']:
            fname += '_noroi'
        if model_cfg['blinded']:
            fname += '_blinded'
        fname += shift_suffix

        fit_io.write_fit(fit_comps, trace, outdir=fitdir, fname=fname, data=data)

    return m2model, trace, fit_comps, data


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    """Parse command line arguments for a single fit point."""
    parser = argparse.ArgumentParser(
        description='Run a single PyMC Bayesian fit for one parameter point.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
            Examples:
            # Single interactive fit with default cuts
            python -m qpym2.run.run_single_fit --co-shift 0.8

            # Specific model point, unblinded
            python -m qpym2.run.run_single_fit \\
                --co-cut gamma10_v150 --ndbd-cutoff u1740 --co-shift 1.0

            # Blinded fit with ROI exclusion
            python -m qpym2.run.run_single_fit \\
                --blinded true --skip-roi 2515 2540

            # Test prior generation
            python -m qpym2.run.run_single_fit --test-priors
        """,
    )

    parser.add_argument(
        '--environment', '-e', default=None,
        help='Hydra environment config name (default from config: local)')
    parser.add_argument(
        '--analysis', '-a', default=None,
        help='Analysis profile name (default from config: base)')
    parser.add_argument(
        '--model-key', '-m', default='default',
        help='Model key in config/models/ (default: default)')

    # Model point parameters
    parser.add_argument(
        '--co-cut', default=None,
        help='Co-60 calibration cut key (e.g. v150, gamma10, gamma10_v150)')
    parser.add_argument(
        '--ndbd-cutoff', default=None,
        help='Signal-region lower-edge cutoff key (e.g. u1740, u1775)')
    parser.add_argument(
        '--co-shift', type=float, default=0.0,
        help='Co-60 energy scale shift in keV (default: 0.0)')

    # Blinding
    parser.add_argument(
        '--blinded', default='false',
        help='Run on blinded data. Accepts true/false or a custom HDF5 '
             'dataset name (default: false)')

    # ROI exclusion
    parser.add_argument(
        '--skip-roi', nargs=2, type=float, metavar=('MIN', 'MAX'),
        help='Exclude an energy range from fit inputs (e.g. --skip-roi 2515 2540)')

    # MCMC settings
    parser.add_argument(
        '--chains', type=int, default=8,
        help='Number of MCMC chains (default: 8)')
    parser.add_argument(
        '--draws', type=int, default=20000,
        help='Number of MCMC draws per chain (default: 20000)')
    parser.add_argument(
        '--mkchain-samples', type=int, default=100,
        help='Number of JAGS MC chain samples for priors (default: 100)')

    # Parallelism
    parser.add_argument(
        '--single-core', action='store_true', default=False,
        help='Run MCMC with cores=1 (no subprocess workers). Use in packed '
             'mode where GNU Parallel already provides node-level parallelism.')

    # Output
    parser.add_argument(
        '--no-write', action='store_true', default=False,
        help='Disable writing fit results to disk (dry-run)')

    # Diagnostics
    parser.add_argument(
        '--test-priors', action='store_true',
        help='Test prior generation and exit')
    parser.add_argument(
        '--prior-width', choices=['5sig', '10sig', 'syst', '2syst'], default='5sig',
        help='Prior sigma mode: 5sig = 5×raw chain std, 10sig = 10×raw chain std, '
             'syst = 1× correlation-matrix perturbation std, '
             '2syst = 2× correlation-matrix perturbation std (previous syst behavior). '
             'Default: 5sig')

    return parser.parse_args()


def _parse_blinded(value):
    """Parse the --blinded flag into a bool or custom dataset name.

    Returns:
        False, True, or a string (custom HDF5 dataset name).
    """
    if value.lower() in ('false', '0', 'no'):
        return False
    if value.lower() in ('true', '1', 'yes'):
        return True
    # Interpret as a custom HDF5 dataset name
    return value


def test_prior_generation(cfg, fit_dir, h2uv_model, comp_cfg, mksamples, sigma_scale=1.0,
                           co_cut=None, ndbd_cutoff=None, co_shift=0.0,
                           use_sample_std=False):
    """Run prior generation as a diagnostic and print results."""
    mcfg = model_cfg(h2uv_model, co_cut=co_cut,
                     ndbd_cutoff=ndbd_cutoff, co_shift=co_shift)
    fitdir = f"{fit_dir}/{mcfg['dir_name']}"
    hm = mcfg['h2uv_model']
    shift_suffix = mcfg.get('shift_suffix', '')

    data_name, _ = get_data_names(cfg)
    mctable, signal_df, data, edges = read_input(f"{fitdir}/input{shift_suffix}.h5", dataname=data_name)

    comp_cfg_new = make_priors(comp_cfg, mctable, mksamples,
                               sigma_scale=sigma_scale,
                               use_sample_std=use_sample_std)

    print("\n--- Prior diagnostics ---")
    for comp in comp_cfg_new:
        print(f"  {comp.get('name', '?'):30s}  "
              f"mean={comp['mean']:.2f}  std={comp['std']:.2f}")
    print("--- Done ---\n")


def main():
    """Run a single fit for one parameter point."""
    args = parse_args()

    # Load Hydra configuration
    cfg = get_config(environment=args.environment, analysis=args.analysis)

    # Derive paths and objects from config
    jags_fpath = get_jags_path(cfg)
    data_name, data_blinded_name = get_data_names(cfg)
    fit_dir = get_fit_dir(cfg)
    h2uv_model = get_h2uv_model(cfg, model_key=args.model_key)
    comp_cfg = get_comp_cfg(cfg, model_key=args.model_key)

    # Prepare JAGS MC chain samples for prior computation
    chain = read_mkchain(jags_fpath)
    if args.prior_width in ('syst', '2syst'):
        corrmat_fpath = get_corrmat_path(cfg)
        mksamples = get_mcsyst_samples(chain, corrmat_fpath=corrmat_fpath, n=args.mkchain_samples).T
        sigma_scale = 1.0 if args.prior_width == 'syst' else 2.0
        use_sample_std = True
    else:
        mksamples = chain.sample(n=args.mkchain_samples).T
        sigma_scale = 5.0 if args.prior_width == '5sig' else 10.0
        use_sample_std = False
    print(f"[run_single_fit] prior_width={args.prior_width}  sigma_scale={sigma_scale}")

    if args.test_priors:
        co_cut = None if args.co_cut in (None, 'None', 'null') else args.co_cut
        ndbd_cutoff = (None if args.ndbd_cutoff in (None, 'None', 'null')
                       else args.ndbd_cutoff)
        test_prior_generation(cfg, fit_dir, h2uv_model, comp_cfg, mksamples,
                               sigma_scale=sigma_scale,
                               co_cut=co_cut, ndbd_cutoff=ndbd_cutoff,
                               co_shift=args.co_shift,
                               use_sample_std=use_sample_std)
        return

    # Parse model point parameters
    blinded = _parse_blinded(args.blinded)
    skip_roi = tuple(args.skip_roi) if args.skip_roi else False

    # Normalise "None" strings from CLI/TSV to actual None
    co_cut = None if args.co_cut in (None, 'None', 'null') else args.co_cut
    ndbd_cutoff = (None if args.ndbd_cutoff in (None, 'None', 'null')
                   else args.ndbd_cutoff)

    # Build the single-point model grid (reuses the same function as toys)
    fit_models = _build_model_grid(
        cfg, h2uv_model,
        co_cuts=[co_cut],
        ndbd_cutoffs=[ndbd_cutoff],
        co_shifts=[args.co_shift],
        skip_roi=skip_roi,
        blinded=blinded,
    )
    mcfg = fit_models.cfg.values.flatten()[0]

    fitdir = f"{fit_dir}/{mcfg['dir_name']}"
    write = not args.no_write

    naming = cfg.infrastructure.naming_conventions
    print(f"[run_single_fit] model={mcfg['name']}  dir={mcfg['dir_name']}  "
          f"write={write}  chains={args.chains}  draws={args.draws}")

    _, trace, fit_comps, _ = setup_and_run(
        fitdir, mcfg, data_name, data_blinded_name,
        mksamples=mksamples, write=write, comp_cfg=comp_cfg,
        chains=args.chains, draws=args.draws,
        single_core=args.single_core, sigma_scale=sigma_scale,
        use_sample_std=use_sample_std,
        input_name=naming.get('input_name', 'input'),
        input_noroi_name=naming.get('input_noroi_name', 'input_noroi'))

    print("\n--- Posterior means ---")
    for var in fit_comps.index:
        print(f"  {var:30s}  posterior_mean={trace.posterior[var].values.mean():.2f}")
    print("--- Done ---\n")


if __name__ == '__main__':
    main()
