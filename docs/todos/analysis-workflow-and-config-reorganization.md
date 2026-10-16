# TODO: organize analysis campaigns and configuration

## Goal

Make it straightforward to define, launch, resume, and audit a collection of model fits, toys, and systematic studies without duplicating parameter choices and paths across runners and notebooks.

The working analysis pattern is:

```text
raw inputs
  → staged ROOT data (reused when compatible)
  → model-point input HDF5
  → fit artifacts
  → toys and systematic variations
  → combined result / plots
```

This plan builds on the current command-line and SLURM machinery. It does not require replacing Hydra, xarray, or TSV task grids in one step.

## What the current workflow tells us

- `run_create_staging.py` produces ROOT staging files that downstream model and systematic workflows read repeatedly.
- `run_create_input.py` expands combinations of `co_cut`, `ndbd_cutoff`, and `co_shift`; `create_m2uv.create_input()` reads staged events, applies model filters/transforms, makes histograms, and writes HDF5.
- Fits consume those inputs and produce `fit.h5` / `fit.nc` artifacts. Toy runs require the fitted artifacts; the bias runner also requires the input edges.
- Toy runs have their own design parameters (type, count, injection range, seed). These are conceptually separate from selecting the physics model or defining a systematic variation.
- Systematics have different dependency patterns. Co-shift and binning studies make related fit points; smoothing can reuse fixed inputs; LSS variations transform staged events; MC systematics vary MC/prior samples. They should share campaign tracking without being forced into one execution algorithm.
- Recent `notebooks/old/build_code/` examples demonstrate reusable-input comparisons:
  - `bin_size_syst_5.ipynb` compares bin sizes using the same staging and data selection, rebuilding templates and fits.
  - `smoothing_syst_2_5.ipynb` keeps input/model/prior settings fixed and varies smoothing at fit time.
  - `matrix_bias_test.ipynb` checks 72 bias-matrix points across fit family, bin size, and Co-60 shift.
  - `qpym2/run/run_bin_size_fits.py` is a useful prototype: it configures bin-size and shift axes, isolates each point in an output directory, and reuses nominal-bin priors. It still hardcodes config-group choices and is specific to that matrix.

## Priority 1 — version the staging layer and make reuse explicit

### Separate staging configuration from model configuration

Add a staging config group describing only the transformation and selection needed to produce reusable staged data:

- input data / MC sources and file patterns;
- input tree and branch mapping;
- output tree, columns, and coordinate transforms;
- common quality/event-selection cuts;
- broad preselection limits needed to keep stage files practical;
- staging schema/version and expected output set.

Do not make a model's fit ROI, Co-60 fit cut, NDBD cutoff, bin size, smoothing, or toy settings part of the staging identity unless they change the stored staged events or columns.

**Compatibility issue to resolve:** `build_staging_config()` currently selects a model and derives a QROOT cut from `model.roi` (`qpym2/cfg/config_builders.py:45-78`). Consequently, a stage made with one ROI/threshold may not be reusable for a model needing events excluded by that stage. Decide on a sufficiently broad reusable stage preselection, and apply model-specific ROI/threshold in input construction. If a broad stage is not viable, make compatibility explicit and require a new stage version when the model lies outside its preselection.

### Give staged artifacts immutable identities and provenance

- Compute a `staging_id` from the resolved staging config, raw input identity (paths plus available size/checksum/version information), staging code/schema version, and relevant source revision.
- Write to a versioned/immutable location such as `staging/<staging_id>/`; never silently overwrite a stage used by completed downstream runs.
- Store a machine-readable manifest beside the output with the resolved config, input files, cuts/transforms, branch schema, creation time/source revision, and per-file event counts.
- Validate an existing stage's manifest and schema before reuse. Expose an explicit rebuild/new-version operation rather than treating path existence as proof of compatibility.
- Make each input manifest record the exact `staging_id` it consumed.

## Priority 2 — separate configuration by responsibility

The current tier boundaries mix concerns: `analysis` contains shared physics definitions and model variants; `hyperparams` mixes a selected model point, study axes, and fit options; and `toys` mixes toy execution defaults with grids of models. Reorganize gradually into these responsibilities:

| Config group | Owns | Examples |
|---|---|---|
| `environments/` | Where/how the job runs and where source/output data live | NERSC/local profile, raw data paths, output roots, ROOT tree/branch names, runtime backend settings |
| `staging/` | Versioned raw-to-staged transformation contract | branch mapping, stored columns, quality cuts, broad preselection, staging schema |
| `analysis/` | Physics conventions common across a set of models | shared cut definitions, component grouping, constants and analysis-wide selection rules |
| `models/` | A named model definition | ROI, threshold, bin size, model-specific component/cut overrides |
| `fit_profiles/` | How an input is fitted | sampler defaults, chains/draws/tune, prior mode/scale, smoothing settings, fit-level switches |
| `toys/` | Toy experiment design and execution | bias/sensitivity/toy type, repetitions, injection distribution/range, seed policy, toy sampler settings |
| `systematics/` | A named systematic's method and variation definition | axis/values, transform or refit method, dataset set, reference-point requirements, uncertainty/delta calculation |
| `studies/` or `campaigns/` | A concrete run plan that composes the profiles above | selected models, nominal point, requested fit/toy/systematic tasks, sweep axes, output campaign name |
| `plots/` (optional) | Presentation choices independent of numerical analysis | style, figure defaults; backend remains environment-specific |

Environment-specific Matplotlib backend choices can stay with the environment; figure style and publication defaults should not be tied to the compute environment.

### Prefer explicit model points and study definitions over a generic hyperparameter bag

Treat a **model point** as the resolved values that determine one model input/fit (for example `co_cut`, `ndbd_cutoff`, `co_shift`, model key, and any applicable ROI/binning override). Treat a **study** as the instruction to vary a subset of those values and compare the resulting artifacts. This makes it explicit that a setting such as `co_shift` affects input construction, while `smooth` can affect fitting without requiring new ROOT input.

Example of the intended separation (illustrative, not a committed schema):

```yaml
# studies/bin-size-co-shift.yaml
campaign:
  name: bin-size-co-shift
  environment: m2v5_mu
  staging: m2v5-loose
  analysis: tech_note
  model: tech_note_200
  fit_profile: nominal
  reference:
    co_cut: gamma10_v150
    ndbd_cutoff: u1740
    co_shift: 1.0
    bin_size_keV: 1.0
  sweep:
    co_shift: [0.0, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5]
    bin_size_keV: [1.0, 2.5, 5.0, 10.0]
  tasks:
    - inputs
    - fit
    - bias_toys
  toys:
    bias:
      count: 200
      injection: {distribution: uniform, min: 5, max: 200}
```

In practice, systematics may need a specialized procedure rather than a generic sweep. The study should name that procedure/profile and its required dependencies, rather than pretending every variation is a Cartesian model grid.

## Priority 3 — resolve a campaign once, then plan the dependency graph

Create a planning layer that composes the selected environment, staging, analysis, model, fit/toy profiles, and study definition into a validated resolved campaign. It should:

1. Validate that selected profile names and parameter values exist and are mutually compatible.
2. Classify which artifacts each varied parameter invalidates.
3. Build only required tasks and dependencies; avoid expanding irrelevant axes into one huge Cartesian product.
4. Reuse compatible staging and identical model inputs, with cache keys based on resolved parameters rather than notebook-local assumptions.
5. Emit a human-readable dry-run summary and an immutable task manifest before execution.
6. Generate the current TSV command rows (or a compatible scheduler input) from that manifest, retaining SLURM array and packed execution.
7. Record task state, seed, command, input artifact IDs, output paths, failures, and completion status; support safe resume and explicit force/rebuild behavior.

The manifest should make the dependency graph visible, for example:

```text
stage:<staging_id>
  └── input:<model_point_id>
        └── fit:<fit_id>
              ├── toys:bias:<toy_profile_id>
              ├── toys:sensitivity:<toy_profile_id>
              └── systematic:<systematic_profile_id>
```

For systematic methods that start from staged ROOT rather than nominal inputs, the graph should declare that directly.

## Priority 4 — standardize artifact layout and provenance

Use campaign/model-point/artifact IDs as the primary identity, with parameter values still present as metadata. For example:

```text
<results>/<campaign_id>/
  campaign.yaml                 # submitted intent
  manifest.json                 # resolved config, IDs, tasks, seeds and status
  models/<model_point_id>/
    input.h5
    fit.h5
    fit.nc
    toys/bias/...
    toys/sensitivity/...
  systematics/<systematic_id>/...
  summaries/...
```

The exact paths can follow existing storage conventions, but each artifact should carry or reference its complete resolved parameters and upstream IDs. Avoid relying on parsing parameter values from directory names, or collecting results into a shared CSV without a campaign/reference identifier.

For systems such as the bin-size matrix that intentionally reuse priors from a nominal point, record that prior source ID in each fit's provenance. The same principle applies to bias correction, nominal posterior references, and systematic delta summaries.

## Priority 5 — make notebooks consumers of campaigns

- Keep notebooks for exploratory validation, comparing outputs, plots, and interpretation.
- Provide one reusable notebook helper to load a campaign, select model points, and resolve artifact paths from the manifest.
- Move repeated parameter declarations, run loops, output path assembly, and copy/pasted deltas into campaign definitions or library functions.
- Have notebooks report missing/incomplete tasks from the manifest rather than reconstructing expected matrix points from hardcoded lists.
- Keep exploratory study definitions easy to create; once a study is used for production/results, save the exact config and resolved manifest with outputs.

## Suggested implementation sequence

1. **Write and validate the staging contract.** Decide broad preselection and version fields; add a manifest and compatibility check. Test reuse across the currently relevant model ROIs/thresholds.
2. **Pilot on the existing bin-size × Co-shift matrix.** Move its hardcoded axes/profile choices into a study definition while preserving the current numerical behavior and prior-reuse rule.
3. **Introduce a campaign planner and manifest.** First generate the current runner commands; do not change the SLURM execution backend yet.
4. **Separate analysis-wide settings from model definitions.** Migrate one representative analysis, then make runner APIs require an explicit model key instead of selecting the first model implicitly.
5. **Add fit and toy profiles.** Separate fit sampler/prior choices from toy design, and bind them to a selected model point through the campaign.
6. **Unify systematic definitions and outputs.** Start with parameter-grid studies and then bring specialized MC/LSS procedures into the same manifest/provenance model.
7. **Migrate notebooks as consumers.** Replace hardcoded paths and point lists in the active build-code workflows; leave clearly identified historical notebooks available as records.

## Acceptance criteria

- A staging output can be identified, inspected, and safely reused without relying on its directory name alone.
- Model ROI, binning, or systematic-parameter changes do not trigger staging again unless the staging compatibility check says they do.
- A study can declare its reference point, varied parameters, requested tasks, seeds, and fit/toy settings in one reviewable place.
- A dry run shows the expected number of inputs/fits/toys/systematic tasks and why each task is needed.
- Re-running a campaign resumes valid completed tasks and detects stale or incompatible artifacts.
- Every final summary can trace back to the resolved campaign, model point, fit/toy/systematic artifacts, and relevant seeds.
- Notebooks no longer need to reconstruct campaign outputs using hardcoded path and parameter lists.
