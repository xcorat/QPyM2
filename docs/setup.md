# Setup/Installation steps

## Setup conda environment

The local environment is managed using **micromamba** with packages installed under `/data/cuore/conda`.

### 1. Create micromamba environment
```bash
export MAMBA_ROOT_PREFIX=/data/cuore/conda
export XDG_CACHE_HOME=/data/cuore/conda/cache
/data/cuore/bin/micromamba create -f environment.yml -y
```

### 2. Activate environment
In each new shell:
```bash
source ./activate.sh
```

Or manually:
```bash
export MAMBA_ROOT_PREFIX=/data/cuore/conda
export XDG_CACHE_HOME=/data/cuore/conda/cache
eval "$(/data/cuore/bin/micromamba shell hook --shell bash)"
micromamba activate qpym2
```

The package itself is installed editable (`pip install -e .`).

### 3. Validate
```bash
python -c "import ROOT, pymc, qpym2; print('PyROOT works! ROOT version:', ROOT.gROOT.GetVersion(), 'PyMC:', pymc.__version__)"
```

## Required dependencies

The core analysis code uses dependencies declared in `environment.yml`:
**PyTables** (`tables`), **xarray**, **netcdf4**, **h5py**, and **PyMC / PyTensor / ArviZ**.

## Running tests

```bash
# Pure tests — no external data or ROOT required (67 passing tests)
pytest -m "not integration"

# Integration tests — require the configured environment and NERSC data
pytest -m integration
```

## Smoke test

Run one stage locally before submitting to SLURM (activate environment first):

```bash
source ./activate.sh

python -m qpym2.run.run_create_staging -e local
python -m qpym2.run.run_create_input   -e local
python -m qpym2.run.run_single_fit -e local \
    --co-cut None --ndbd-cutoff None --co-shift 0.0 --chains 8 --draws 20000
```

For the HPC workflow see [`scripts/README.md`](../scripts/README.md).
