== Analysis Process Guide <sec:qpym2_guide>

=== Introduction

==== Project structure:

- *Code*: `./qpym2/` package contains all code for staging, input creation, fitting, and systematics.

- * Setup configuration*: `./config/`
    - Setup the environment, mainly input files, paths and output folders
    - Setup the main analysis config
    - Find more on `config/README.md`

- *Data folder*: `./data/`
    - Usually symlinked to a scratch folder (the subfolders)
    - `fits/`: Contains all input files (`input.h5`, component histograms) and fit outputs (trace `fit.nc`, components `fit.h5`), organized by model directory (e.g., `model_gamma10_v150_u1740/`)
    - `plots/`: Contains all plots (fits, bias, systematics, etc.)
    - `syst_results/`: Systematic fit results, organized by type (lss, mcsyst); includes summaries and traces.
    - `toy_results/`: Toy, sensitivity, and bias fit results (ag gregate CSV results and HDF5)

- *Scripts*: `./scripts/`
    - Job submission and management scripts. The run commands with parameters are in `params`, the `slurm` folder contain submission scripts that sometimes differ slightly (multiprocess vs gnu_parallel, etc.)
    - `gen_params.py`: Generates parameter grid TSV files for batch jobs
    - `submit.sh`: SLURM array job wrapper (use `--test` for debug queue)
    - `submit_syst.sh`: Systematic study submission wrapper
    - `run_job.sh`: Per-task runner (reads command from TSV, sets environment)

- *Notebooks*: `./notebooks/`
    - Analysis and validation notebooks (bias correction, posterior plots, systematic validation)

- *Tech Note*: `./tech_note/`
    - Typst source for this analysis note (`main.typ`)



=== Steps

1. *Create Staging Files*
    - Configure filtering for data and MC in `./qpym2/cfg/staging.py`.
    - Run `run_create_staging.py` either directly  or through the batch scripts.
    - Ensure background model and $0 nu beta beta$ MC inputs are consistent; adjust cuts in config as needed.

2. *Create Fit Inputs*
    - Configure fit input parameters.    
    - *Note:* The ROI differs from the staging selection to include necessary padding.
    - Run `run_create_input.py`. Run in batches (e.g., ~50 models) to manage memory.

3. *Run Fits*
    - Execute fits with `run_fits.py`. Split into batches if needed to prevent crashes. Each fit takes ~1 minute.

4. *Run Bias Fits*
    - Use `run_toy_fits.py`. Modify `main()` to select bias fits (vs. sensitivity).
    - Run only for selected models (specific Cobalt/NDBD cuts).

5. *Run Systematics Fits*
    - Use `run_syst.py`. Update `main()` to execute the desired study:
      - Line shape scaling: `run_lss_syst()`
      - Background model variations: `run_mcsyst()`
      - Signal shift systematic: `run_ndbdshift_syst()`