# QPyM2 Configuration System (v0.6)

QPyM2 uses [Hydra](https://hydra.cc/) + [OmegaConf](https://omegaconf.readthedocs.io/) loaded programmatically via `qpym2.cfg.hydra_adapter.HydraConfigManager`.

Configuration is structured into four explicit, decoupled tiers that compose together into a complete analysis configuration.

---

## Configuration Tiers

```
config/
├── config.yaml              # Master entry point & default composition
├── environments/            # Tier 1: Infrastructure, paths, and system settings
│   ├── base.yaml            # Abstract baseline (paths, output directories, qroot conventions)
│   ├── local.yaml           # Local workstation environment profile (data_root: /data/cuore/Analysis/MC)
│   └── nersc.yaml           # NERSC Perlmutter environment profile (data_root: /global/cfs/...)
├── staging/                 # Tier 2: Staging data contract
│   └── default.yaml         # Raw-to-staged transform contract (cuts, source mappings)
├── analysis/                # Tier 3: Physics baseline & conventions
│   └── base.yaml            # Baseline background component groupings and standard cuts
└── models/                  # Tier 4: Physics models & ROI
    └── tech_note.yaml       # Named analysis models (tech_note_200, tech_note_250, etc.)
```

### 1. Environments (`config/environments/`)
Controls system-level paths, data roots, raw input directories, and output locations.
- **`base.yaml`**: Declares standard directory variables, file naming conventions, `qroot` ROOT tree definitions, and default output paths (including `output_paths.artifact_store: "artifact_store"`).
- **`local.yaml`**: Sets `infrastructure.data_root` to local storage (e.g. `/data/cuore/Analysis/MC`).
- **`nersc.yaml`**: Configured for NERSC Perlmutter execution paths under `/global/cfs/...`.

### 2. Staging (`config/staging/`)
Defines the data transformation contract to produce staged ROOT `uvtree` files from raw data.
- **`default.yaml`**: Contains preselection cuts:
  - `broad_preselection`: `"TotalEnergy >= 2390 && TotalEnergy <= 2570"`
  - `base`: `"Multiplicity == 2"`
  - `multiplicity`: `"TimeDiff >= 0 && TimeDiff <= 100e-3"`
  - `data_quality`: `"BaselineQuality == 1 && PulseQuality == 1"`
  And source file match patterns (`bkg_mc`, `ndbd_mc`, `data`).
- *Decoupling principle:* Staging cuts do not depend on downstream model-specific ROI, `co_cut`, `ndbd_cutoff`, binning, or energy shifts.

### 3. Analysis (`config/analysis/`)
Defines the baseline physical definitions and standard cut options:
- **`base.yaml`**: Baseline background MC component groupings (`co60`, `u238`, `th232`, `others`) and standard cuts dictionary (`co_cuts: {v150, gamma10, gamma10_v150}`, `ndbd_cutoff: {u1740, u1775}`).

### 4. Models (`config/models/`)
Defines statistical fit models and energy regions of interest (ROI):
- **`tech_note.yaml`**: Contains named model points:
  - `tech_note_200`: `threshold: 200`, `bin_size: 2.5`, `esum_min: 2400`, `esum_max: 2570`.
  - `tech_note_250`: `threshold: 250`, `bin_size: 1.0`, `esum_min: 2400`, `esum_max: 2570`.
  - `tech_note_150`: `threshold: 150`, `bin_size: 1.0`, `esum_min: 2400`, `esum_max: 2570`.
  - `tech_note_3comp`: Merged background grouping variation (`u_th` component).

---

## Path Resolution Engine

Paths in configuration are resolved using predictable rules:
1. **Explicit absolute paths**: If a path starts with `/`, it is treated as absolute and passed through unchanged.
2. **Relative paths**: Relative paths are resolved against `infrastructure.data_root` (or the `QPYM2_DATA_DIR` environment variable if set).
3. **Artifact Store**: Resolves to `<data_root>/<output_paths.artifact_store>` by default.

---

## Programmatic Loading

QPyM2 loads configurations programmatically without `@hydra.main` decorators:

```python
from qpym2.cfg.hydra_adapter import HydraConfigManager, get_config

# Load default composition (environments=local, staging=default, analysis=base, models=tech_note)
cfg = HydraConfigManager.get_config()

# Load with environment and CLI-style overrides
cfg = HydraConfigManager.get_config(
    environment="nersc",
    overrides=["models.tech_note_200.binning.bin_size=1.0"]
)
```

## Helper Bridges (`qpym2.cfg.config_builders`)

`qpym2.cfg.config_builders` provides bridge functions adapting OmegaConf `DictConfig` structures to the plain dictionaries expected by runners:
- `build_staging_config(cfg)`: Extracts staging paths, tree names, and merged staging cuts.
- `build_input_config(cfg, model_key=...)`: Builds input creation configurations.
- `get_h2uv_model(cfg, model_key=...)`: Instantiates `TempConfig` histogram model with rotated coordinate bounds ($u, v$).
- `get_comp_cfg(cfg, model_key=...)`: Retrieves component groupings, allowing model-specific overrides.
- `get_artifact_store_dir(cfg)`: Resolves the centralized artifact store directory path.

---

## Legacy Config Preservation

Pre-v0.6 configuration files (which contained coupled `hyperparams/` and `toys/` tiers, as well as per-model environment YAMLs) are backed up under `config_old/` for historical reference.
