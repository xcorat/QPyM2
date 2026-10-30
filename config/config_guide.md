# QPyM2 Configuration Guide (v0.6)

This guide documents the modern configuration architecture in QPyM2. It reflects the 4-tier Hydra layout in `config/` designed to decouple execution environments, data contracts, and physical modeling.

---

## 1. Architecture

QPyM2 uses **Hydra + OmegaConf**, loaded **programmatically** via `qpym2.cfg.hydra_adapter_new.HydraConfigManager` (or `qpym2.cfg.hydra_adapter`). The entry point is `config/config.yaml`, which composes four orthogonal tiers:

```yaml
# config/config.yaml
pipeline_version: "0.6.0"

defaults:
  - environments: local   # Tier 1: System infrastructure, data roots, output paths
  - staging: default      # Tier 2: Staging data contract (cuts, input patterns)
  - analysis: base        # Tier 3: Physics baseline (cuts dictionary, component groupings)
  - models: tech_note     # Tier 4: Named statistical fit models & ROI
  - _self_
```

All tier files use `@package _global_` and standard `.yaml` extensions.

| Tier | Directory | Scope | Example Targets |
|---|---|---|---|
| **Environment** | `config/environments/` | Machine paths, data roots, storage paths, system settings | `base`, `local`, `nersc` |
| **Staging** | `config/staging/` | Raw data matching, filtering, preselection cut contract | `default` |
| **Analysis** | `config/analysis/` | Baseline physics definitions, standard cut alternatives | `base` |
| **Models** | `config/models/` | Named model configurations (ROI, threshold, bin size) | `tech_note` |

---

## 2. Tier Details

### Tier 1: Environments (`config/environments/`)
Contains machine-specific directory structures and path settings.
- `base.yaml`: Abstract layout defining standard directory variables, file naming conventions, `qroot` ROOT tree definitions, and default output paths (e.g., `output_paths.artifact_store: "artifact_store"`).
- `local.yaml`: Configures `infrastructure.data_root: "/data/cuore/Analysis/MC"` for local workstation development.
- `nersc.yaml`: Configures `infrastructure.data_root: "/global/cfs/cdirs/cuore/..."` for NERSC Perlmutter execution.

### Tier 2: Staging Contract (`config/staging/`)
Specifies the transformation contract for producing staged `uvtree` ROOT files from raw data.
- Staging default preselection:
  ```yaml
  cuts:
    broad_preselection: "TotalEnergy >= 2390 && TotalEnergy <= 2570"
    base: "Multiplicity == 2"
    multiplicity: "TimeDiff >= 0 && TimeDiff <= 100e-3"
    data_quality: "BaselineQuality == 1 && PulseQuality == 1"
  ```
- **Decoupling principle**: Staging preselection must be broad enough to encompass all downstream model ROIs, but never dependent on individual model thresholds, cutoffs, or energy shifts.

### Tier 3: Analysis Conventions (`config/analysis/`)
Defines shared component classifications and standard cuts dictionary:
- `cuts.co_cuts`: Standard selection variants (`v150`, `gamma10`, `gamma10_v150`).
- `cuts.ndbd_cutoff`: Standard $0\nu\beta\beta$ cutoff selections (`u1740`, `u1775`).
- `components`: Standard background groupings (`co60`, `u238`, `th232`, `others`).

### Tier 4: Models & ROI (`config/models/`)
Defines concrete fit configurations keyed by name:
- `tech_note_200`: `threshold: 200`, `bin_size: 2.5`, `esum_min: 2400`, `esum_max: 2570`.
- `tech_note_250`: `threshold: 250`, `bin_size: 1.0`, `esum_min: 2400`, `esum_max: 2570`.
- `tech_note_150`: `threshold: 150`, `bin_size: 1.0`, `esum_min: 2400`, `esum_max: 2570`.
- `tech_note_3comp`: Component variation merging Uranium and Thorium.

---

## 3. Path Resolution Engine

Paths are resolved using deterministic rules through `qpym2.cfg.hydra_adapter_new.resolve_data_path`:
1. **Explicit absolute paths**: Any path starting with `/` is treated as an absolute path and returned as-is.
2. **Relative paths**: Any relative path is resolved relative to `cfg.infrastructure.data_root`.
3. **Environment variable override**: If `QPYM2_DATA_DIR` is set in the shell environment, it overrides `data_root` as the base directory for relative paths.

---

## 4. Programmatic Usage

Loading configuration in Python:

```python
from qpym2.cfg.hydra_adapter_new import HydraConfigManager, get_config

# Load default composition
cfg = get_config()

# Override environment or model parameters
cfg = get_config(
    environment="nersc",
    model="tech_note",
    overrides=["models.tech_note_200.binning.bin_size=1.0"],
)
```

Bridging to pipeline dictionaries:

```python
from qpym2.cfg.config_builders_new import (
    build_staging_config,
    build_input_config,
    get_h2uv_model,
    get_artifact_store_dir,
)

staging_dict = build_staging_config(cfg)
input_dict = build_input_config(cfg, model_key="tech_note_200")
h2uv_model = get_h2uv_model(cfg, model_key="tech_note_200")
store_path = get_artifact_store_dir(cfg)
```
