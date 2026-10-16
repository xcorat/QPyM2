# Setup/Installation steps

## Setup conda environment

1. Create conda environment using the environment.yml file
    `mamba env create -f environment.yml`

    *If the base environment doesn't have mamba, create and install mamba prior, and then update,*
    `conda create -n xqpy mamba -y`
    `conda activate xqpy`
    `mamba env update -f environment.yml`

2. Activate conda environment
    `conda activate xqpy`

    The package itself is installed editable via `environment.yml`
    (`pip: -e .`). To install manually: `pip install -e .`.

3. Validate
    + `python -c "import ROOT; h = ROOT.TH1F('test', 'test', 10, 0, 10); print('PyROOT works!')"`

## Required dependencies

The core analysis code additionally needs packages that are not all declared in
`setup.py`: **PyTables** (`tables`, used by `pandas.to_hdf/read_hdf`),
**xarray**, **netcdf4**, and **h5py**. PyMC/PyTensor/ArviZ are required for
fitting. If an import fails at runtime, check these first.

## Running tests

```bash
# Pure tests — no data or ROOT required
pytest -m "not integration"

# Integration tests — require the configured environment and NERSC data
pytest -m integration
```

Some tests are currently known-broken (see the Testing section of
[`AGENTS.md`](../AGENTS.md)); a full test run may error at collection until
those are fixed.

## Smoke test

Run one stage locally before submitting to SLURM (activate `xqpy` first):

```bash
python -m qpym2.run.run_create_staging -e m2v5-1 -a tech_note
python -m qpym2.run.run_create_input   -e m2v5-1 -a tech_note
python -m qpym2.run.run_single_fit -e m2v5-1 -a tech_note \
    --co-cut None --ndbd-cutoff None --co-shift 0.0 --chains 8 --draws 20000
```

For the HPC workflow see [`scripts/README.md`](../scripts/README.md).
