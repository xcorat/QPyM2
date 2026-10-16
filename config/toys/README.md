# Toy Fit Configuration

This directory contains configuration files for toy fit runs. Toy configs follow the same pattern as other QPyM2 config tiers (environments/, analysis/, hyperparams/): `base.yaml` provides defaults, with optional override variants.

## Configuration Files

### base.yaml
Default configuration for interactive toy fit runs.
- **Minimal model grid**: Single parameter point for quick testing
- **Interactive output paths**: Local CSV and log files
- **Usage**: Automatically loaded when no `--toys` argument is provided

### batch.yaml
Configuration variant for batch systematic scans.
- **Full model grid**: Multiple co_cuts, ndbd_cutoffs, and co_shifts for comprehensive scanning
- **Batch output paths**: Results organized by type
- **Usage**: Use with `--toys batch` argument

## Configuration Schema

All toy configs share the same structure:

```yaml
# Bias study signal injection ranges (counts)
bias:
  min: 5        # Minimum injected signal counts
  max: 200      # Maximum injected signal counts (null = 2x posterior mean)

# Output configuration
output:
  csv: "path/to/output.csv"      # CSV output path
  logs: "path/to/logs.log"        # Log file path
  append: true                    # Append to existing files (false = fail on exist)

# Model parameter grid
model_grid:
  co_cuts:           # Co-60 calibration cut selections
    - gamma10_v150   # Combined gamma-10 + v150 cut
  ndbd_cutoffs:      # Signal region cutoff selections  
    - u1740          # u1740 cutoff
  co_shifts:         # Co-60 energy shifts (keV)
    - 1.0            # Energy shift in keV
```

## Usage Examples

### Interactive Toy Fits (using base.yaml)
```bash
# Uses default base.yaml config
python -m qpym2.run.run_toy_batch --type sens --n 5

# Run full grid with run_toy_fits.py
python -m qpym2.run.run_toy_fits --type all --n 200
```

### Batch Scans (using batch.yaml)
```bash
# Use batch config for full parameter grid
python -m qpym2.run.run_toy_fits --toys batch --type all --n 200

# Single batch job with batch config
python -m qpym2.run.run_toy_batch --toys batch --type bias --n 10
```

### SLURM Batch Jobs
```bash
# Submit array job using base config (default)
sbatch --array=0-99 scripts/submit_toys.sh sens 5

# Submit using batch config
sbatch --array=0-99 scripts/submit_toys.sh sens 5 --toys batch

# Override specific parameters
sbatch --array=0-49 scripts/submit_toys.sh bias 10 --toys batch --bias-min 10 --bias-max 500
```

### Sequential Runner
```bash
# Run full grid sequentially with batch config
scripts/run_toy_jobs_sequential.sh -c batch -n 200 -t all

# Quick test with base config
scripts/run_toy_jobs_sequential.sh -n 10 -t sens
```

## Command-Line Argument Pattern

Toy configs follow the same pattern as other QPyM2 config tiers:

```python
# Hydra configs (environments, analysis, hyperparams)
get_config(environment='nersc', analysis='tech_note', hyperparams='bias_study')

# Toy configs (toys)
load_toy_config('batch')  # Or None for base.yaml
```

Command-line arguments:
- `--environment <name>` → loads `config/environments/<name>.yaml`
- `--analysis <name>` → loads `config/analysis/<name>.yaml`
- `--hyperparams <name>` → loads `config/hyperparams/<name>.yaml`
- `--toys <name>` → loads `config/toys/<name>.yaml` (default: base)

## Creating Custom Toy Configs

To create a new toy config variant:

1. Copy `base.yaml` or `batch.yaml` as a template
2. Modify the parameters as needed
3. Save as `config/toys/<your_name>.yaml`
4. Use with `--toys <your_name>`

Example custom config for a focused systematic study:
```yaml
# config/toys/co_systematic.yaml
bias:
  min: 5
  max: 200

output:
  csv: "data/toy_results/co_systematic.csv"
  logs: "logs/co_systematic.log"
  append: true

model_grid:
  co_cuts:
    - gamma10_v150
  ndbd_cutoffs:
    - u1740
  co_shifts:         # Focus on Co-60 energy shifts
    - 0.0
    - 0.5
    - 1.0
    - 1.5
    - 2.0
```

Usage:
```bash
python -m qpym2.run.run_toy_fits --toys co_systematic --type all --n 500
```

## File Existence Handling

The `append` parameter controls file behavior:
- `append: true` (default): Appends to existing files, creates if missing
- `append: false`: Fails if output files already exist (prevents accidental overwrites)

Override from command line:
```bash
# Force non-append mode (fail if files exist)
python -m qpym2.run.run_toy_batch --type sens --no-append
```

## Parameter Override Priority

Parameters are resolved in this order (highest to lowest priority):
1. Command-line arguments (e.g., `--bias-min 10`)
2. Toy config file (e.g., `batch.yaml`)
3. Built-in defaults

Example:
```bash
# bias_min=5 from batch.yaml
python -m qpym2.run.run_toy_batch --toys batch --type bias

# bias_min=10 (CLI overrides config)
python -m qpym2.run.run_toy_batch --toys batch --type bias --bias-min 10
```
