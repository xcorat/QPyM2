# AGENTS.md — QPyM2

Guidance for AI coding agents working in this repository.

## What this is

QPyM2 is the CUORE Multiplicity-2 (M2) analysis framework. It converts raw ROOT
data into staged ROOT files, builds histogram fit inputs, runs PyMC Bayesian
fits, and produces systematics / toy / sensitivity / bias studies and plots.

Language: Python (`qpym2` package, ~20k LOC). Config: Hydra + OmegaConf.
CLIs: stdlib `argparse`. HPC: NERSC Perlmutter via SLURM.

## Environment

- Code runs inside the conda env defined by `environment.yml` (name `xqpy`).
  SLURM scripts instead activate `/global/common/software/cuore/sachi/qpym2`
  (env name `qpym2`). Be aware of this mismatch.
- `python` is **not** guaranteed to be on `PATH` in a bare shell — activate the
  env first. The package is installed editable (`pip install -e .`).
- Hard requirements: **PyROOT/ROOT**, **PyMC >= 5.25 + PyTensor + ArviZ**,
  `pandas` (PyTables/`tables`), `h5py`, `netcdf4`, `matplotlib`, `scipy`.
  `tables` and `xarray` are imported by the code but are not declared in
  `setup.py` — do not assume a clean install has them.
- No GPU/CUDA/JAX. All CPU. `numba` is declared but currently unused.
- Data paths point at `/global/cfs/...`. `data/` is a symlink to
  `../../MC/tech_note/`, and `data/`, `output/`, `notebooks/` are gitignored.

## Pipeline

The pipeline is **TSV-parameter-grid driven**: a generator writes a TSV of
per-task commands; each row is run by a SLURM array task or a local pool.

```
scripts/gen_params.py ─► scripts/params/<type>_params.tsv ─► SLURM array (run_job.sh)
                                                          └► python -m qpym2.run.<stage>
```

| Stage | Entry point | Output |
|---|---|---|
| Params grid | `scripts/gen_params.py <staging\|input\|fit\|toys\|sens\|bias\|syst_lss\|syst_mc>` | `scripts/params/*.tsv` |
| Staging | `qpym2/run/run_create_staging.py` | staged `uvtree` ROOT |
| Fit input | `qpym2/run/run_create_input.py` | `input.h5` |
| Fit input (storage) | `qpym2/run/run_create_storage_inputs.py` | `fit_data/.../input/*.h5` |
| Fits | `qpym2/run/run_single_fit.py` (atomic), `qpym2/run/run_fits.py` (local pool) | `fit.nc`, `fit.h5` |
| Fits (storage) | `qpym2/run/run_single_storage_fit.py` | `fit_data/.../output/fit_<UUID>.h5` |
| Studies | `qpym2/run/run_study.py` | `studies/<name>/{results.csv, .sqlite, bias_calibration.json}` |
| Systematics | `qpym2/run/run_syst_lss.py`, `run_syst_mc.py`, `run_ndbd_shift_syst.py` | `.nc` / CSV |
| Toys/sens/bias | `qpym2/run/run_toy_batch.py` (SLURM task), `run_toy_fits.py` (library) | per-task CSV |
| Aggregate | `scripts/collect_results.py` | merged CSV |
| Plots | `qpym2/plotting/*.py` (library only, called from notebooks) | figures |

### Canonical commands

```bash
# Activate env first
conda activate xqpy

# 1. Generate a parameter grid (writes scripts/params/<type>_params.tsv)
python scripts/gen_params.py fit --environment m2v5-1 --analysis tech_note \
    --chains 8 --draws 20000 --prior-width 2syst
python scripts/gen_params.py toys --n-fits 200

# 2. Run one stage locally
python -m qpym2.run.run_create_staging -e m2v5-1 -a tech_note
python -m qpym2.run.run_create_input   -e m2v5-1 -a tech_note
python -m qpym2.run.run_create_storage_inputs -e local -a base -m tech_note_200 --co-shifts 0.0 0.8 1.0 1.2
python -m qpym2.run.run_single_storage_fit -e local -a base -m tech_note_200 --co-shift 0.8 --save-trace
python -m qpym2.run.run_study --study-dir studies/example_co_shift --check-inputs
python -m qpym2.run.run_study --study-dir studies/example_bias_study --check-bias
python -m qpym2.run.run_study --study-dir studies/example_sens_study --check-sens
python -m qpym2.run.run_single_fit -e m2v5-1 -a tech_note \
    --co-cut None --ndbd-cutoff None --co-shift 0.0 --chains 8 --draws 20000
python -m qpym2.run.run_fits -p scripts/params/fit_params.tsv --dry-run

# 3. Submit to SLURM (array) / packed
bash scripts/submit.sh fit --chains 8 --draws 20000
bash scripts/submit_packed.sh staging
bash scripts/submit_syst.sh lss
```

CLI conventions: every runner is `parse_args() -> main()` with
`-e/--environment` and (usually) `-a/--analysis`. `run_syst_lss.py` uses
argparse subparsers (`sscl` / `dq`).

## Configuration

Hydra + OmegaConf loaded **programmatically**, not via `@hydra.main` and not
exposed as CLI flags. `qpym2.cfg.hydra_adapter_new.get_config(environment=, analysis=, model=, ...)`
calls `compose()` over `config/config.yaml`.

- Master defaults (`config/config.yaml`): `environments: local`,
  `staging: default`, `analysis: base`, `models: tech_note`.
- Tiers: `config/environments/*.yaml` (infrastructure, paths, data_root),
  `config/staging/*.yaml` (raw-to-staged transform contract; default `broad_preselection: "TotalEnergy >= 2390 && TotalEnergy <= 2570"`),
  `config/analysis/*.yaml` (baseline physics cuts dictionary & component groupings),
  `config/models/*.yaml` (named statistical fit models, ROI & binning).
- Legacy pre-v0.6 configs are backed up in `config_old/`.
- Actual layout uses Hydra `defaults:` composition and `.yaml` extensions. See
  `config/README.md` and `config/config_guide.md`; the live YAML is the source of truth.
- `qpym2/cfg/config_builders_new.py` adapts `DictConfig` into the plain dicts the
  runners expect (`build_staging_config`, `build_input_config`, `get_h2uv_model`,
  `get_comp_cfg`, `get_artifact_store_dir`, ...). Path resolution is rooted at
  `infrastructure.data_root` (or `QPYM2_DATA_DIR`).

## Testing

`pyproject.toml` configures pytest: `testpaths = ["qpym2"]`, the `integration`
marker is registered, and the default `addopts` is `-m 'not integration'`. Tests
split into:

- Pure (no data/ROOT): `qpym2/tests/test_staging_transforms.py`,
  `test_make_fik_digitize.py`, `test_fit_syst.py`, `test_hydra_adapter.py`.
- Integration (`@pytest.mark.integration`): need real config + `/global/...`
  data — `qpym2/tests/test_run_fits.py`.

```bash
pytest -m "not integration"      # pure tests (safe, no data)
pytest -m integration            # requires configured env + NERSC data
```

The suite was pruned so that only tests importing cleanly against the current
tree remain.

## Code conventions

- **Logging**: two coexisting styles. Prefer the custom
  `qpym2.utils.logging` API (`log(level, ...)`, `debug`, global `DEBUG`,
  `start_timer`/`log_elapsed`; levels `< 0` are skipped). A few newer modules use
  stdlib `logging.getLogger(__name__)`. Runners also use `print()`.
- **Docstrings**: Google-style (`Args:` / `Returns:`), with reST `:func:` refs in
  newer modules. Older core modules have minimal docstrings.
- **Types**: inconsistent. Newer modules (`cfg/`, `utils/pool.py`, `io/csv.py`,
  `syst/`) are typed; `fit/`, `models/`, most `io/` are not. No mypy.
- **Errors**: `ValueError` / `FileNotFoundError` / `RuntimeError`; some broad
  `except`. No custom exception hierarchy.
- Do not add comments unless asked. Follow existing module patterns.

## Lint / format

`pyproject.toml` configures ruff (conservative ruleset: `E4,E7,E9,F`). The
`ruff` and `pytest` executables live in the `xqpy` conda env.

```bash
ruff check .
ruff format --check .      # apply with: ruff format .
pytest -m "not integration"
pytest -m integration
```

Ruff ignores backup variants (`*.bak`, `*.old`, `*.remove`). Do not do a
repo-wide reformat as part of unrelated changes.

## Legacy cleanup

Stale/broken code was removed from this repository. Do not re-introduce backup
or duplicate variants (`foo copy.py`, `foo.py.old`, `foo.bak`, `*DO_NOT_USE*`);
use git branches or a scratch directory instead.

Canonical modules to prefer: `run_create_staging.py`, `run_create_input.py`,
`run_single_fit.py`, `run_fits.py`, `io/hdf.py`, `io/rdf.py`,
`io/staging.py`, `fit/create_m2uv.py`, `fit/params.py`, `fit/fit.py`,
`models/m2.py`, `cfg/hydra_adapter.py`, `cfg/config_builders.py`,
`storage/paths.py`, `storage/backend.py`, `storage/hasher.py`, `storage/inputs.py`, `storage/study.py`.

## Repository state

`main` has many uncommitted working-tree changes (including deletions under
`tech_note/`). Check `git status` / `git diff` before assuming HEAD reflects the
working tree. Never commit unless explicitly asked.

## Nested guidance

See `qpym2/AGENTS.md`, `config/AGENTS.md`, `scripts/AGENTS.md`, and
`notebooks/AGENTS.md` for area-specific details.
