# qpym2/ — package layout

See the root `AGENTS.md` for env setup, the pipeline, config system, and the
canonical test/lint commands. This file maps the package.

## Module map

| Dir | Purpose | Canonical entry points |
|---|---|---|
| `cfg/` | Hydra config loading + adapters | `hydra_adapter.py` (`get_config`), `config_builders.py`, `params_grid.py`, `model_utils.py`, `toy_config.py`, `syst_config.py` |
| `run/` | Pipeline stage runners (argparse `parse_args()`/`main()`) | `run_create_staging.py`, `run_create_input.py`, `run_create_storage_inputs.py`, `run_single_storage_fit.py`, `run_study.py`, `run_single_fit.py`, `run_fits.py`, `run_toy_batch.py`, `run_toy_fits.py`, `run_syst_lss.py`, `run_syst_mc.py`, `run_ndbd_shift_syst.py` |
| `io/` | ROOT / HDF5 / CSV I/O | `rdf.py` (ROOT RDataFrame), `hdf.py` (pandas HDF), `staging_transforms.py`, `csv.py`, `toys.py` |
| `fit/` | Fit construction, priors, results | `create_m2uv.py`, `params.py`, `fit.py`, `results.py`, `bias.py`, `syst.py`, `syst_lss.py` |
| `models/` | PyMC model definitions | `m2.py` |
| `syst/` | Systematic-uncertainty helpers | `core.py`, `rate.py`, `smear.py`, `io.py`, `plots.py` |
| `plotting/` | Plot helpers (library only, no CLI) | `posterior_plots.py`, `projection_plots.py`, `post2.py`, `syst_posterior.py`, `data_plots.py` |
| `utils/` | Shared utilities | `logging.py`, `hists.py`, `hist_model.py`, `pool.py`, `mathy.py`, `hw.py`, `TempConfig.py` |
| `core/` | Lower-level ROOT/hist primitives | `qroot.py`, `hists.py` |
| `storage/` | Data storage layout, hashing, HDF5 metadata, input generation, fit assembly, and study indexing | `schema.py`, `hasher.py`, `paths.py`, `backend.py`, `inputs.py`, `fits.py`, `study.py` |
| `studies/` | Study folder config inheritance, assembly, input validation, fit execution, bias toy fits, systematics, and rate marginalisation | `config.py`, `point.py`, `study.py`, `assembly.py`, `toys.py`, `syst_runner.py`, `marginalize.py` |
| `tests/` | Test suite | see root `AGENTS.md` "Testing" |

### `studies/` layout

- `point.py`: `StudyPoint` dataclass, path resolvers, point-level runner/collector delegates.
- `study.py`: `Study` container, study manifest and input builders, campaign delegates.
- `assembly.py`: Study directory inspection and instantiation (`assemble_study`).
- `config.py`: Hierarchical OmegaConf study config inheritance (`study.yaml`).
- `toys.py`: Bias & sensitivity toys generation, execution, progress tracking, and calibration result collection.
- `syst_mc.py`: Background model Monte Carlo systematic fits and posterior sample extraction with prior re-weighting.
- `syst_lss.py`: Lineshape (SSCL & dQ) simulations, completion checks, and delta parameter extractions.
- `syst_bm.py`: Background model variant systematics (threshold 200 vs 250 keV, binning alternatives).
- `syst_template.py`: Template likelihood binning and smoothing variant fits and delta extractions.
- `syst_runner.py`: Dispatcher and aggregator for multi-point systematic campaigns, delta tables, and reporting.
- `marginalize.py`: Prior smearing, split-Gaussian integration, PSA marginalisation, and limit extraction.
- `_common.py`: Shared atomic file writes, task status inspectors, posterior mode finders, summary stats extractors.


## Conventions

- `run/__init__.py` documents each runner. `run_fits` is a working local
  `multiprocessing.Pool` runner (not a deprecated shim).
- Logging: prefer the custom `qpym2.utils.logging` API (`log`, `debug`, global
  `DEBUG`, `start_timer`/`log_elapsed`). Some newer modules use stdlib logging.
- Runners are dual-use: importable library functions (`run_fits()`,
  `setup_and_run()`, `run_toy_fit()`) plus a thin CLI. Keep that separation.
- No custom exception hierarchy; use `ValueError` / `FileNotFoundError` /
  `RuntimeError`.

## Do not reintroduce

Do not re-add backup/duplicate variants (`*.bak`, `*.old`, `* copy*`,
`*DO_NOT_USE*`) or unimportable legacy runners; use git branches or a scratch
directory instead.
