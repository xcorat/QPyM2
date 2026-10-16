# QPyM2 Configuration Guide

This guide documents how configuration actually works in QPyM2. It reflects the
live Hydra layout in `config/`. (An earlier version of this file described an
aspirational `extends:` / `tools/*.yml` scheme that was never implemented.)

## Architecture

QPyM2 uses **Hydra + OmegaConf**, loaded **programmatically** via the compose
API — not through the `@hydra.main` decorator and not exposed as CLI config
flags. The master file is `config/config.yaml`, which composes four tiers:

```yaml
# config/config.yaml
defaults:
  - environments: base     # Tier 1: infrastructure (paths, system, qroot)
  - analysis: models_base  # Tier 2: physics (models, ROI, cuts, binning)
  - hyperparams: base      # Tier 3: systematic selections/overrides
  - toys: base             # Tier 4: toy-study grid
  - _self_
```

All tiers use `@package _global_` and the `.yaml` extension.

| Tier | Directory | Contains | Example names |
|---|---|---|---|
| Environment | `config/environments/` | Absolute paths, `qroot` read config, `system` | `base`, `tech_note`, `nersc`, `m2v5-1`, `m2v5-1-250` |
| Analysis | `config/analysis/` | `models.*` ROI/threshold/binning, `cuts.*`, `components` | `models_base`, `tech_note`, `m2`, `m2_2700` |
| Hyperparams | `config/hyperparams/` | `co_cuts`, `ndbd_cutoff`, `co_shift`, `ndbd_shift`, `blinded`, `skip_roi`, `smooth` | `base`, `bias_study`, `tech_note`, `m2uv` |
| Toys | `config/toys/` | Bias range, output paths, `model_grid` | `base`, `batch` (see `toys/README.md`) |

The tiers are selected by name. `environment` and `analysis` are nearly always
passed by the CLI runners; `hyperparams` and `toys` default to `base`.

## Loading configuration

The public interface is `qpym2.cfg.get_config` (a thin wrapper around
`HydraConfigManager.get_config`):

```python
from qpym2.cfg import get_config

cfg = get_config(
    environment="m2v5-1",
    analysis="tech_note",
    hyperparams=None,          # optional, default base
    overrides=["hyperparams.co_shift=0.5"],
)

staging_dir = cfg.infrastructure.output_paths.staging_dir
threshold   = cfg.models.tech_note.roi.threshold   # model key from analysis yaml
```

- `environment=`, `analysis=`, `hyperparams=` add Hydra overrides
  `environments=<x>`, `analysis=<x>`, `hyperparams=<x>`.
- `overrides=[...]` appends arbitrary Hydra-style overrides.
- `plots_dir` is derived from `fit_dir/plots` when not explicitly set
  (`_apply_output_paths_defaults`).
- Toy configs are loaded separately:
  `qpym2.cfg.load_toy_config("batch")` (simple YAML, see `toys/README.md`).

Helpers in `HydraConfigManager`: `list_options()`, `get_default_selections()`,
`get_config_dict()`, `print_config()`, `set_config_dir()`.

### From the CLI

Runners expose `-e/--environment` and `-a/--analysis` and forward them to
`get_config`. They do **not** expose the rest of the config as flags; change
those in the YAML.

```bash
python -m qpym2.run.run_create_input -e m2v5-1 -a tech_note
```

## Adapters for runners

`qpym2/cfg/config_builders.py` converts the `DictConfig` into the plain dicts
the pipeline code expects. Prefer these over reaching into `cfg` directly:

| Function | Purpose |
|---|---|
| `build_staging_config` / `build_mc_config` / `build_ndbd_config` / `build_data_config` | Stage input dicts |
| `build_input_config` / `build_write_config` | Input creation + output |
| `get_h2uv_model` / `get_comp_cfg` | Model + component config |
| `get_fit_dir` / `get_plots_dir` / `get_syst_dir` / `get_toy_results_dir` / `get_data_dir` | Output paths |
| `get_jags_path` / `get_data_names` / `get_corrmat_path` | Input data paths/names |

`qpym2/cfg/params_grid.py` is the single source of truth for grid axes
(co_cuts, ndbd_cutoffs, co_shifts). `qpym2/cfg/model_utils.py` builds model
directory names and parses them back.

## Model definitions

Models live under `models:` in the analysis tier, keyed by name. Example
(abridged from `config/analysis/models_base.yaml` + a model variant):

```yaml
cuts:
  co_cuts:
    v150: "v>150"
    gamma10: "abs(e1-1334)>10 && abs(e2-1173)>10"
    gamma10_v150: "abs(e1-1334)>10 && abs(e2-1173)>10 && v>150"
  ndbd_cutoff:
    u1740: "u>1740"
    u1775: "u>1775"

models:
  tech_note:
    roi:
      esum_min: 2400
      esum_max: 2570
      threshold: 200
    binning:
      bin_size: 1
```

- Rotated-axis (`u`/`v`) limits are **derived in code** from ROI
  (`esum_min`, `esum_max`) and `threshold`; do not configure them directly.
- Set only `models.*.binning.bin_size` (defined on the `e1`/`e2` axes).
- `hyperparams.co_cuts` selects a key from `cuts.co_cuts`; `hyperparams.ndbd_cutoff`
  selects a key from `cuts.ndbd_cutoff`.
- Systematics are produced by varying these hyperparameters (`co_shift`,
  `ndbd_shift`, cut selections) — see `qpym2/cfg/params_grid.py`.

## Best practices

- Keep **paths** in `environments/` and **physics** in `analysis/`.
- Select config by name; avoid editing the master `config.yaml` for one-off runs.
- Treat the YAML as the source of truth.
