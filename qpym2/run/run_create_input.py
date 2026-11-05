"""
Create fit input HDF5 files from staging ROOT files for M2 analysis.

Fit inputs contain histogrammed data and MC components in 2D u,v space,
prepared for Bayesian fitting with PyMC.

Usage:
    python run_create_input.py [--environment ENV] [--analysis ANALYSIS]
    python run_create_input.py --co-cut gamma10_v150 --ndbd-cutoff u1740 --co-shifts 0.8 1.0 1.2
    python run_create_input.py --ndbd-syst  # Create NDBD systematic variations
"""

import argparse
import math
import time
import xarray as xr
import tracemalloc

from qpym2.cfg.hydra_adapter import get_config
from qpym2.cfg.config_builders import (
    get_jags_path, build_input_config, build_write_config, get_h2uv_model
)
from qpym2.cfg.params_grid import ALL_CO_CUTS, ALL_NDBD_CUTS, ALL_CO_SHIFTS
from qpym2.fit.create_m2uv import create_input
from qpym2.fit.params import model_cfg
from qpym2.utils.logging import debug


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Create fit input HDF5 files from staging ROOT files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default tech_note environment
  python run_create_input.py
  
  # Run with specific environment and analysis profile
  python run_create_input.py --environment production --analysis threshold_350
  
  # Run one (co_cut, ndbd_cutoff) point with multiple co_shifts (SLURM mode)
  python run_create_input.py --co-cut gamma10_v150 --ndbd-cutoff u1740 --co-shifts 0.8 1.0 1.2
  
  # Create NDBD systematic variations
  python run_create_input.py --ndbd-syst
        """
    )
    
    parser.add_argument(
        '--environment', '-e',
        default=None,
        help='Environment name (default from config: local)'
    )
    
    parser.add_argument(
        '--analysis', '-a',
        default=None,
        help='Analysis profile name (default from config: base)'
    )

    parser.add_argument(
        '--model-key', '-m',
        default='default',
        help='Model key in config/models/ (default: default)'
    )
    
    parser.add_argument(
        '--ndbd-syst',
        action='store_true',
        help='Create NDBD systematic variations'
    )
    
    parser.add_argument(
        '--blinded',
        action='store_true',
        help='Use blinded data (NOT FULLY IMPLEMENTED)'
    )
    
    # Single grid-point args (for SLURM array parallelism)
    parser.add_argument(
        '--co-cut',
        default=None,
        help='Single co_cut value (e.g. gamma10_v150).'
    )
    parser.add_argument(
        '--ndbd-cutoff',
        default=None,
        help='Single ndbd_cutoff value (e.g. u1740)'
    )
    parser.add_argument(
        '--co-shifts',
        nargs='+',
        type=float,
        default=None,
        help='One or more co_shift values (e.g. 0.8 1.0 1.2). '
             'Shared ROOT reads are cached across shifts within the same task.'
    )
    
    return parser.parse_args()





def main():
    """Create fit input files for main analysis."""
    args = parse_args()
    
    # Load Hydra configuration
    cfg = get_config(
        environment=args.environment,
        analysis=args.analysis
    )
    
    # Build configuration dicts using builders
    jags_fpath = get_jags_path(cfg)
    input_cfg = build_input_config(cfg, model_key=args.model_key)
    write_cfg = build_write_config(cfg)
    
    # Get h2uv_model from config
    h2uv_model = get_h2uv_model(cfg, model_key=args.model_key)
    
    # Extract cut dictionaries from config for model_cfg
    co_cuts_dict = dict(cfg.cuts.co_cuts)
    ndbd_cuts_dict = dict(cfg.cuts.ndbd_cutoff)
    
    if args.blinded:
        print('WARNING: Blinded data mode not fully tested!')
        # TODO: Enable blinded mode properly
        # Note: This requires proper implementation in create_input
    
    # Set DEBUG level
    from qpym2.utils.logging import DEBUG
    DEBUG = 0
    
    # Single grid-point mode (for SLURM array jobs)
    single_point = (args.co_cut is not None or args.ndbd_cutoff is not None
                    or args.co_shifts is not None)
    if single_point:
        co_cuts      = [args.co_cut]
        ndbd_cutoffs = [args.ndbd_cutoff]
        co_shifts    = args.co_shifts if args.co_shifts is not None else [0.0]
    else:
        # Default grid (full-grid / single-node mode)
        co_cuts      = ALL_CO_CUTS[3:4]    # ['gamma10_v150']
        ndbd_cutoffs = ALL_NDBD_CUTS[2:3]  # ['u1740']
        co_shifts    = ALL_CO_SHIFTS[3:]   # [0.8, 1.0, 1.2, 1.5]

    fit_models = xr.Dataset(
        coords={
            'ndbd_cutoff': ndbd_cutoffs,
            'co_cuts': co_cuts,
            'co_shift': co_shifts,
        },
    )

    fit_models['cfg'] = xr.apply_ufunc(
        lambda coc, ndc, cos: model_cfg(h2uv_model,
                co_cut=coc,
                ndbd_cutoff=ndc,
                co_shift=cos,
                co_cuts_dict=co_cuts_dict,
                ndbd_cuts_dict=ndbd_cuts_dict),
            fit_models.co_cuts, fit_models.ndbd_cutoff, fit_models.co_shift, vectorize=True)
    
    # Create fit inputs — shared cache avoids re-reading ROOT files
    # for components whose (fpath, defs, filters) are identical across co_shifts.
    hist_cache = {}
    t_total = time.perf_counter()
    for cfg in fit_models.cfg.values.flatten():
        debug('running:', cfg['name'])
        t0 = time.perf_counter()
        create_input(cfg, input_cfg, jags_fpath, fit_name=None,
                     write_cfg=write_cfg, hist_cache=hist_cache)
        print(f"  [{cfg['name']}] done in {time.perf_counter() - t0:.1f}s")
    print(f"Total time: {time.perf_counter() - t_total:.1f}s")


def main_blinded():
    """Create fit input files for blinded analysis.
    
    TODO: This function needs proper testing and blinded data implementation.
    """
    print('ERROR: Blinded mode requires proper implementation.')
    print('Remove this check after verifying blinded data is correctly configured.')
    exit(1)
    
    args = parse_args()
    
    # Load Hydra configuration
    analysis = args.analysis if args.analysis else args.environment
    cfg = get_config(
        environment=args.environment,
        analysis=analysis
    )
    
    # Build configuration dicts using builders
    jags_fpath = get_jags_path(cfg)
    input_cfg = build_input_config(cfg)
    write_cfg = build_write_config(cfg)
    
    # Get h2uv_model from config
    h2uv_model = get_h2uv_model(cfg)
    
    from qpym2.utils.logging import DEBUG
    DEBUG = 0
    
    co_cuts      = ALL_CO_CUTS[3:4]    # ['gamma10_v150']
    ndbd_cutoffs = ALL_NDBD_CUTS[2:3]  # ['u1740']
    co_shifts    = ALL_CO_SHIFTS[3:]   # [0.8, 1.0, 1.2, 1.5]
    
    fit_models = xr.Dataset(
        coords={
            'ndbd_cutoff': ndbd_cutoffs,
            'co_cuts': co_cuts,
            'co_shift': co_shifts,
        },
    )

    # For blinded fits- need to implement blinded=True parameter
    # fit_models['cfg'] = xr.apply_ufunc(
    #     lambda coc, ndc, cos: model_cfg(h2uv_model,
    #             co_cut=coc,
    #             ndbd_cutoff=ndc,
    #             co_shift=cos,
    #             blinded=True),
    #         fit_models.co_cuts, fit_models.ndbd_cutoff, fit_models.co_shift, vectorize=True)
    
    # For NO ROI inputs:
    # ROI_MIN, ROI_MAX = 2515, 2540
    # fit_models['cfg'] = xr.apply_ufunc(
    #     lambda coc, ndc, cos: model_cfg(h2uv_model,
    #             co_cut=coc,
    #             ndbd_cutoff=ndc,
    #             co_shift=cos,
    #             skip_roi=(ROI_MIN, ROI_MAX)),
    #         fit_models.co_cuts, fit_models.ndbd_cutoff, fit_models.co_shift, vectorize=True)
    
    fit_models['cfg'] = xr.apply_ufunc(
        lambda coc, ndc, cos: model_cfg(h2uv_model,
                co_cut=coc,
                ndbd_cutoff=ndc,
                co_shift=cos),
            fit_models.co_cuts, fit_models.ndbd_cutoff, fit_models.co_shift, vectorize=True)
    
    hist_cache = {}
    for cfg in fit_models.cfg.values.flatten():
        debug('running:', cfg['name'])
        create_input(cfg, input_cfg, jags_fpath, fit_name=None,
                     write_cfg=write_cfg, hist_cache=hist_cache)


def ndbd_syst_main():
    """Create NDBD systematic variation inputs."""
    args = parse_args()
    
    # Load Hydra configuration
    analysis = args.analysis if args.analysis else args.environment
    cfg = get_config(
        environment=args.environment,
        analysis=analysis
    )
    
    # Build configuration dicts using builders
    jags_fpath = get_jags_path(cfg)
    input_cfg = build_input_config(cfg)
    write_cfg = build_write_config(cfg)
    
    # Get h2uv_model from config
    h2uv_model = get_h2uv_model(cfg)
    
    from qpym2.utils.logging import DEBUG
    DEBUG = 0
    
    # Fixed model configuration with NDBD shift variations
    coc, ndc, cos = 'gamma10_v150', 'u1740', 1.0
    ndbd_shifts = [0, 0.3, 0.5, 0.8, 1, 1.2, 1.5][1:5]
    
    cfgs = [model_cfg(h2uv_model, co_cut=coc, ndbd_cutoff=ndc, co_shift=cos, ndbd_shift=nds)
            for nds in ndbd_shifts]
    
    hist_cache = {}
    for cfg in cfgs:
        debug('running:', cfg['name'])
        create_input(cfg, input_cfg, jags_fpath, fit_name=None,
                     write_cfg=write_cfg, hist_cache=hist_cache)


if __name__ == '__main__':
    args = parse_args()
    
    if args.ndbd_syst:
        ndbd_syst_main()
    elif args.blinded:
        main_blinded()
    else:
        main()
