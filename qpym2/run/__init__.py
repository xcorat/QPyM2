"""
QPyM2 Run Module — pipeline stage runners.

Modules:
    run_create_staging: Create staging ROOT files from raw data.
    run_create_input: Create fit input HDF5 files from staging.
    run_single_fit: Atomic runner for a single Bayesian fit point (CLI + library).
    run_fits: Local multiprocessing.Pool runner that executes every row of a
        fit_params.tsv via run_single_fit.setup_and_run.
    run_toy_fits: Run toy / sensitivity / bias fits for validation.
    run_toy_batch: SLURM array task wrapper for run_toy_fits.
    run_syst_lss: Systematic uncertainty fits from LSS variations.
    run_syst_mc: Systematic uncertainty fits from background MC variations.
    run_ndbd_shift_syst: Systematic fits for the NDBD energy-scale shift.
"""
