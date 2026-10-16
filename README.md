# QPyM2 — CUORE M2 Analysis Framework

Framework for the CUORE Multiplicity-2 (M2) analysis: staging raw ROOT data,
building histogram fit inputs, running PyMC Bayesian fits, and producing
systematics, toy/sensitivity/bias studies, and plots.

## Overview

The pipeline is driven by per-stage parameter grids (TSV files) whose rows are
executed on NERSC Perlmutter via SLURM array jobs (or locally for testing):

```
raw ROOT ─► staging ROOT ─► input.h5 ─► fit.nc/.h5 ─► systematics ─► toys/sens/bias ─► plots
```

- Package: `qpym2/` (Python). Config: Hydra + OmegaConf. CLIs: `argparse`.
- HPC: NERSC Perlmutter via `scripts/slurm/*.slurm` + `scripts/run_job.sh`.

## Setup

Create the conda environment and install the package:

```bash
mamba env create -f environment.yml
conda activate xqpy
```

See [`docs/setup.md`](docs/setup.md) for validation steps, and note that
**PyROOT/ROOT**, **PyMC/PyTensor/ArviZ**, and PyTables (`tables`) are required.

## Quickstart

```bash
conda activate xqpy

# Generate a fit parameter grid
python scripts/gen_params.py fit --environment m2v5-1 --analysis tech_note \
    --chains 8 --draws 20000 --prior-width 2syst

# Run one stage locally
python -m qpym2.run.run_create_staging -e m2v5-1 -a tech_note
python -m qpym2.run.run_create_input   -e m2v5-1 -a tech_note
python -m qpym2.run.run_single_fit -e m2v5-1 -a tech_note \
    --co-cut None --ndbd-cutoff None --co-shift 0.0 --chains 8 --draws 20000

# Submit to SLURM
bash scripts/submit.sh fit --chains 8 --draws 20000
```

See [`scripts/README.md`](scripts/README.md) for the full SLURM command
reference.

## Testing

```bash
pytest -m "not integration"   # pure tests (no data/ROOT required)
pytest -m integration         # requires configured env + NERSC data
```

Some tests are currently known-broken; see the testing section of
[`AGENTS.md`](AGENTS.md).

## Documentation

- [`AGENTS.md`](AGENTS.md) — architecture, commands, conventions, gotchas
  (also nested in `qpym2/`, `config/`, `scripts/`, `notebooks/`).
- [`docs/setup.md`](docs/setup.md) — installation and environment validation.
- [`config/README.md`](config/README.md) — configuration layout.
- [`scripts/README.md`](scripts/README.md) — parameter grids and job submission.
