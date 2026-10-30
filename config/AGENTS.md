# config/ — Hydra configuration (v0.6)

See the root `AGENTS.md` and `config/README.md` for how config is loaded programmatically.

## Architecture & Layout

`config/config.yaml` is the master entry point and composes four decoupled tiers via Hydra `defaults:`:

```yaml
defaults:
  - environments: local   # Tier 1: System infrastructure, data roots, output paths
  - staging: default      # Tier 2: Raw-to-staged data transformation contract
  - analysis: base        # Tier 3: Physics baseline (component groupings, cut definitions)
  - models: tech_note     # Tier 4: Named statistical models (ROI, threshold, binning)
  - _self_
```

- **`environments/*.yaml`** — machine-specific paths and infrastructure only: `infrastructure.data_root`, `infrastructure.raw_data_dir`, `output_paths.*`, `qroot`. Never place physics cuts or model parameters here.
- **`staging/*.yaml`** — data transformation contract for staging raw ROOT to `uvtree`. Staging default: `broad_preselection: "TotalEnergy >= 2390 && TotalEnergy <= 2570"`. Staging cuts must remain decoupled from downstream model ROI and fit variations.
- **`analysis/*.yaml`** — physics baseline conventions: component groupings (`components`), cut definitions (`cuts.co_cuts`, `cuts.ndbd_cutoff`).
- **`models/*.yaml`** — named fit models (`models.<name>`): ROI (`esum_min`, `esum_max`, `threshold`), binning (`bin_size`).

## Programmatic Access

Load configurations using `qpym2.cfg.hydra_adapter_new` (or `qpym2.cfg.hydra_adapter`):
```python
from qpym2.cfg.hydra_adapter_new import get_config

cfg = get_config(environment="local", analysis="base", model="tech_note")
```

Adapt to pipeline dicts using `qpym2.cfg.config_builders_new`:
```python
from qpym2.cfg.config_builders_new import (
    build_staging_config,
    build_input_config,
    get_h2uv_model,
    get_artifact_store_dir,
)
```

## Editing Rules

- Keep paths (`environments/`) and physics (`analysis/`, `models/`) strictly separated.
- Never add derived $u/v$ rotated coordinate limits to YAML; $u/v$ limits are derived in code from `esum_min`, `esum_max`, and `threshold`.
- Only configure binning via `models.<name>.binning.bin_size` on the $e_1, e_2$ axes.
- Relative paths in `environments/*.yaml` automatically resolve against `infrastructure.data_root` (or `QPYM2_DATA_DIR`). Explicit absolute paths pass through untouched.
- Legacy configuration is archived in `config_old/`. Do not re-introduce monolithic or deprecated YAML files.
