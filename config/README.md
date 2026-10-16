# Configuration Settings for the M2 Analysis.

This directory contains the Hydra/OmegaConf configuration that drives the M2
analysis. The master file is [`config.yaml`](config.yaml), which composes four
tiers via Hydra `defaults:`:

```
environments  →  paths, system settings, qroot read config
analysis      →  models (ROI, threshold, binning), cuts, components
hyperparams   →  systematic selections/overrides (co_cuts, co_shift, ...)
toys          →  toy-study grid (bias range, output paths)
```

1. **Environments** (`environments/`): computational environment — data/output
   paths, file patterns, naming conventions, ROOT read config. `base.yaml` is the
   default; select another with `environment=<name>`.

2. **Analysis** (`analysis/`): physics content — model definitions (`models.*`
   ROI/threshold/binning), cut definitions, and background component groupings.
   `models_base.yaml` is the default; select another with `analysis=<name>`.

3. **Hyperparameters** (`hyperparams/`): systematic variation selections applied
   on top of a model (`co_cuts`, `ndbd_cutoff`, `co_shift`, `ndbd_shift`,
   `blinded`, `skip_roi`, `smooth`). Default `base` (no variation).

4. **Toys** (`toys/`): parameters for toy/sensitivity/bias runs, loaded via
   `qpym2.cfg.load_toy_config`. See [`toys/README.md`](toys/README.md).

Files are selected by name through `qpym2.cfg.get_config(...)`. Add a new variant
by copying the closest existing file and registering/selecting it by name — do
not edit `config.yaml` for one-off runs.

The shared matplotlib style is [`qm2style.mplstyle`](qm2style.mplstyle).

Refer to the [configuration guide](config_guide.md) for details and
[AGENTS.md](AGENTS.md) for agent-oriented notes.

For active `h2uv` models, set only `models.*.binning.bin_size` (defined on the
`e1`/`e2` axes). Rotated-axis (`u`/`v`) limits are derived from ROI (`esum_min`,
`esum_max`) and `threshold` in code, and should not be configured directly.
