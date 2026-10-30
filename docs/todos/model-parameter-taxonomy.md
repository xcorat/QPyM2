# Model parameters, artifact identity, and study variation

- **Status:** Implemented in v0.6.0
- **Date:** 2026-10-17
- **Reference:** [`docs/todos/data-management-and-artifact-store.md`](data-management-and-artifact-store.md), [`docs/ai/reviews/2026-10-17-data-management-review.md`](../ai/reviews/2026-10-17-data-management-review.md)
- **Related:** [`docs/todos/analysis-workflow-and-config-reorganization.md`](analysis-workflow-and-config-reorganization.md), [`docs/todos/v0.6-implementation-list.md`](v0.6-implementation-list.md), [`docs/todos/fit-data-layer.md`](fit-data-layer.md)

> **Fit-data follow-up:** Object-level input sharing, fit/result sets, terminology, and reversible profile-based IDs are specified in [`fit-data-layer.md`](fit-data-layer.md) and [`docs/data-layer-spec.md`](../data-layer-spec.md). That follow-up refines fit-input/result identity; staging remains outside its scope.

---

## 1. Problem

Today "model" is overloaded across three different concepts:

| Sense | Where it lives | What it controls |
|---|---|---|
| Statistical/physics model | `config/analysis/*.yaml`, `config/models/*.yaml` → `models.*`, `components`, `roi`, `binning` | Histogram/template construction and likelihood shape |
| Model point | `hyperparams.*` / CLI flags | Selected cuts, shifts, and other numerical choices |
| Folder label | `model_utils.build_model_dir` → `model_<co_cut>_<ndbd_cutoff>` | Path identity |

Two problems follow:

- `run_single_fit.main()` calls `get_h2uv_model(cfg)` without a `model_key`; `_get_model_key` then returns `next(iter(cfg.models))` (`qpym2/cfg/config_builders.py`). A production task should carry an explicit model key rather than silently choosing the first one.
- The effective fit specification is assembled as an untyped `mcfg` dict from CLI arguments and configuration (`qpym2/fit/params.py:model_cfg`). It is not a complete, validated record of the resolved scientific request or its provenance.

The v0.6 artifact store needs explicit parameter schemas and a clear rule for deriving each artifact identity from its actual inputs and operation settings.

## 2. Design decision: group settings by operation, not by comparison axis

Each parameter must answer two questions:

1. **Operation dependency:** which operation consumes this value, and which outputs can it change?
2. **Study intent:** is this field deliberately varied in this study, or held fixed?

These questions do not require a second taxonomy of axes. An axis enum would classify parameters, but would not define artifact dependencies: fields assigned to one axis may feed different artifacts, while fields assigned to different axes may feed the same operation. The hasher still needs the operation's effective inputs.

Use typed, operation-owned specifications. Each operation's identity is derived from its canonical resolved specification and named upstream artifact identities. A study lists the exact field paths it varies, for example `vary=("fit.smooth",)`. Comparison validation checks the declared fields and protected invariants; it does not infer intent from a category name.

Configuration values may be authored in Hydra YAML. The Python structured schema is the single source of truth for field names, types, defaults, and validation. Avoid a second parameter-registry YAML or parallel enum/dataclass/hasher lists that must be manually synchronized. The resolved, versioned specifications are serialized into the campaign and artifact manifests.

### Data context, layers, and metadata

Three buckets, each with a different role:

1. **Data context (separate data-roots).** `blinded`, the blinding method and the staging coverage define a *data-root*. Data from different roots is never interchangeable, so these are not comparison fields and no study can vary them. Each root has a `context.yaml` and every artifact under it records the `context_id`. Separation by directory makes cross-context comparison impossible by construction, without a validator rule.
2. **Parameters.** Everything else that changes a scientific product (selection, binning, components, smoothing, template transforms, priors, study-method variations). All are typed, comparable, and can be varied against a nominal point.
3. **Run metadata.** Sampler, chains, draws, tune, posterior seed, host, git, SLURM. Not parameters. Stored in the manifest next to the artifact in the same data-root.

Data layers inside one data-root:

| Layer | Contents | Notes |
|---|---|---|
| `staging` | staged ROOT | Derived from raw sources and the data context. |
| `input` | observed-data objects, source-level templates, and input-set descriptors | Derived from fixed staged inputs and input-affecting settings; a monolithic `input.h5` may be a compatibility materialization. |
| `output` | fit results, toy experiments/fits, systematic variations, summaries, plots | Fit results reference input sets and fit specifications; labels are associations, not payload owners. |

Fits, systematic variations and toys are all fit-like outputs. They share the fit contract but are saved under distinct `output/<kind>/` folders (`fit`, `syst`, `toys`, `plots`) and tagged with a `kind` and `role` (`nominal`, `variation:<name>`, `toy:<i>`) in the manifest, so they can be listed and cleaned independently.

## 3. Operation-owned parameter groups

These groups describe current and anticipated operation boundaries. The group is not itself an artifact identity: each artifact identity includes the fields actually consumed by its operation and its named upstream dependencies.

```python
@dataclass(frozen=True)
class StagingParams:
    schema_version: int
    source_set: ...
    cuts: ...
    coverage: ...
    blinding: ...

@dataclass(frozen=True)
class InputParams:
    histogram_type: ...
    bin_size: float
    esum_min: float
    esum_max: float
    threshold: float
    co_cut: str | None
    ndbd_cutoff: str | None
    skip_roi: tuple[float, float] | None
    template_transforms: ...

@dataclass(frozen=True)
class FitParams:
    components: ...
    smooth: ...
    likelihood: ...
    priors: ...

@dataclass(frozen=True)
class RunProfile:
    """Metadata only: recorded in the manifest, never hashed or varied."""
    sampler: ...
    chains: int
    draws: int
    tune: int
    seed: int | None

@dataclass(frozen=True)
class StudySpec:
    nominal: ...
    vary: tuple[str, ...]
    variations: ...
    outputs: ...

@dataclass(frozen=True)
class ExecutionMetadata:
    execution_id: str
    host: str
    git_commit: str | None
    git_dirty: bool | None
    slurm_job_id: str | None
    started_at: str
```

The ellipses are illustrative types, not a proposed open-ended dictionary schema. Define concrete typed structures for source identities, cuts, transformations, priors, likelihoods, and variation methods as their contracts are implemented.

### Parameter placement and effect

| Parameter | Operation/group | Current effect and identity consequence |
|---|---|---|
| Raw source identities, staging cuts, output columns/schema | `StagingParams` | Change the staged artifact. Include the resolved stage contract and source identities in `staging_id`. |
| Staging coverage / broad preselection | `StagingParams.coverage` | Changes which events exist downstream; include it in `staging_id`. It is not the analysis ROI. |
| `blinded` and blinding injection specification | `StagingParams.blinding` / data context | Selects the data-root and is also hashed into `staging_id`. Not a comparison field. Include the blinding method, inputs, and randomization identity as applicable. |
| `esum_min`, `esum_max`, `threshold` | `InputParams` | Define analysis selection and histogram geometry. Include in the histogram/input derivation that consumes them. |
| `co_cut`, `ndbd_cutoff` | `InputParams` | Select events/templates in input construction; include the resolved cut definitions and scope, not just labels. |
| `skip_roi` | `InputParams` | Excludes an energy interval during input selection. Include in the affected input identity; protect against accidental co-plotting with unskipped results by study validation. |
| `bin_size`, `histogram_type` | `InputParams` | Affect histogram bin edges/representation and therefore histogram/input products. `histogram_type` may be fixed by a profile, but remains identity-relevant. |
| `co_shift`, `ndbd_shift` | `InputParams.template_transforms` | Transform selected background or signal templates during histogram construction. Include transform parameters and target scope in the affected template/input derivation. |
| `components` | `FitParams` | Groups source-level MC histograms into likelihood components during fit construction. It changes the fit problem, not the current source-granularity `input.h5`. |
| `smooth` | `FitParams` | Current code smooths assembled background and signal histograms during fit construction (`qpym2/models/m2.py`, `qpym2/fit/fit.py`). It changes the fit problem, not the current source-level `input.h5`. If smoothing is moved upstream later, revise the operation contract then. |
| Likelihood/model definition | `FitParams` | Changes the statistical fit problem and its fit-problem identity. |
| Priors and prior source/reference | `FitParams.priors` | Resolved distributions, parameters, and reference/source affect the fit problem. A label such as `5sig` alone may not identify the realized prior. |
| Sampler, chains, draws, tune, target acceptance, initialization, posterior-sampling seed | `RunProfile` (metadata) | **Not parameters.** Recorded in the manifest run block; never hashed, never varied in a study. Cache reuse checks them as a minimum-quality policy (section 5). |
| Host, git/dirty state, SLURM IDs, timestamps | `ExecutionMetadata` | Record where/how an attempt ran. Keep out of scientific derivation identity unless an explicit reproducibility policy makes a code/environment identity an operation dependency. |

### Advisory groups (clarity only)

Fields may carry a non-identity `group` label so they can be displayed, indexed and varied together. A group never affects hashing or validation.

| Group | Fields |
|---|---|
| `structure` | `bin_size`, `histogram_type`, `components` |
| `selection` | `esum_min`, `esum_max`, `threshold`, `co_cut`, `ndbd_cutoff`, `skip_roi` |
| `transforms` | `co_shift`, `ndbd_shift`, `smooth`, study-method variations (LSS-DQ, LSS-scaling, background normalization) |
| `priors` | `priors` |

`vary` always lists exact field paths. A helper may expand `group:transforms` into its field paths, but the expanded list is what is validated and recorded. A group is a convenience label, not an axis.

## 4. Systematics and study variations

“Systematic” describes the purpose of a study, not one operation boundary. Do not force every systematic method into one `SystParams` object or one generic refit implementation. A study variation declares its method, parameters, reference, and dependencies; it can invoke input building, template transformation, fitting, or aggregation as appropriate.

Examples based on the current code:

| Study method | What it consumes/changes | Identity boundary |
|---|---|---|
| Co/NDBD energy shift | Transforms selected templates before or during histogram construction (`get_fit_transformers`) | Affected template/input product; downstream fit refers to that product |
| LSS-DQ / LSS-scaling | Reads staged MC events, transforms event coordinates, rebuilds histograms, then fits with nominal priors (`fit/syst_lss.py`) | Variation-specific transformed-template derivation and fit result; the observed data/staging may be reused, but this is not simply a refit of an unchanged complete input |
| Background-model normalization sampling | Samples JAGS source normalizations, optionally perturbed by the configured covariance, and rebuilds grouped templates for fits (`run_syst_mc.py`, `io/bkg_model.py`) | Variation/sample specification and fit result; record each sampled vector or its deterministic seed and generation recipe |
| Toy generation | Generates a distinct pseudo-dataset from a specified model/design | Toy artifact includes its generating model/input, toy design, index, and generation seed; its subsequent fit is a separate fit artifact |
| Systematic summary | Aggregates a reference result and named variations | Summary identity includes the exact reference/variation artifact IDs and aggregation method/version |

The existing `mc_stat` name is ambiguous. The current `run_syst_mc.py` samples background-source normalizations from a JAGS chain, optionally with covariance perturbations; it is not the usual bin-wise finite-MC template-statistics nuisance model. Prefer a precise method name such as `background_normalization` or `bm_normalization`. Reserve `finite_mc_template_stats` for an implementation that models finite histogram/template statistics.

## 5. Artifact identity and cache policy

Derive IDs from operation contracts, not from broad parameter categories:

```python
staging_id = derive_id(
    operation="staging",
    params=resolved_staging_params,
    upstream={"context": context_id, "raw": raw_source_identities},
)

histogram_id = derive_id(
    operation="histogram_build",
    params=resolved_input_params,
    upstream={"staging": staging_id},
)

fit_problem_id = derive_id(
    operation="fit_problem",
    params=resolved_fit_params,
    upstream={"histograms": histogram_id, "observations": observation_id},
)
# The fit output (trace + manifest) is stored under fit_problem_id.
# The run profile is recorded in the manifest, not hashed.
```

This is the v0.6 hash-based identity example. For fit-data objects and sets, use the reversible, versioned identity-profile/default contract in [`docs/data-layer-spec.md`](../data-layer-spec.md); the identity still has to distinguish any setting or member reference that changes the scientific product.

These IDs are illustrative; where products are packed together, the manifest can describe multiple logical products without requiring one file per ID. Identity follows the actual operation graph and the effective resolved settings.

Distinguish:

- **Fit-problem identity:** likelihood, observations/templates, components/grouping, smoothing, priors, and other settings that define the posterior target.
- **Fit output identity:** fit-problem identity only. The output folder is addressed by `fit_problem_id`; it holds the trace and its manifest.
- **Run profile (metadata, not parameters):** sampler, chains, draws, tune, target acceptance, initialization, and the posterior-sampling seed. They describe how the trace was produced. They never enter any ID, are not comparison fields, and cannot appear in `vary`. They are stored in the manifest's run block next to the artifact, in the same data-root.
- **Execution-attempt metadata:** host, job ID, start time, retry number, git state. Same treatment as the run profile.

Stochastic *scientific* draws are different from the posterior-sampling seed. A seed that realizes a blinding injection, prior generation, a transform/variation sample, or a toy dataset is part of that operation's recipe and therefore of its identity, because it creates a different input. Only the posterior-sampling seed is run-profile metadata. Name seeds by role; a single generic `seed` field is not sufficient provenance.

Cache reuse is a policy check against the manifest, not an identity check. An existing output satisfies a request when its contract, schema and checksums match and its recorded run profile meets the request's minimum (draws, chains, convergence diagnostics). If the profile is insufficient, run again (or pass `--force`) and replace or version the trace under the same `fit_problem_id`.

This follows the data-management review's rule: a fit-only change must not invalidate upstream data, while any effective selection, transform, or representation that changes a stored product must produce the appropriate downstream identity (`docs/ai/reviews/2026-10-17-data-management-review.md`, §2 and §5).

## 6. Study comparison and protected choices

A study names the exact fields it intentionally varies:

```python
study = StudySpec(
    name="smoothing_scan",
    nominal=NominalPoint(
        model_key="tech_note_200",
        staging=...,
        input=...,
        fit=...,
    ),
    vary=("fit.smooth",),
    variations=...,  # explicit values/points
    outputs=("fit_summary",),
)
```

The planner validates the **resolved task graph**, not a taxonomy enum:

- Every requested field path exists and is typed.
- Every point differs from the nominal only in declared fields and their downstream derived consequences.
- Non-varied scientific inputs are fixed.
- Required upstream artifacts, references, and prior sources are explicit.
- All points live in one data-root. Blinded and unblinded data cannot be mixed because they are in different roots. `skip_roi` remains a protected parameter: it must match across ordinary comparisons, and an exceptional study may request an explicit policy override. A `group:` expansion must not silently authorize varying a protected field.
- Runtime metadata may differ without making points scientifically incomparable.

The `vary` list expresses study intent and provides plot/grouping fields. It does not control artifact hashing; operation specifications and upstream IDs do that.

## 7. Consequences for v0.6

- Define typed operation specs for staging, input/histogram construction, fit problem, sampling, and study/variation requests in `qpym2/storage/schemas.py` or a small set of focused schema modules.
- Use the Python structured schema as the source of truth; use Hydra YAML for values and named profiles, not as a duplicate parameter registry.
- Implement canonical serialization and operation-specific derivation IDs from resolved specs and named upstream artifact IDs.
- Store complete resolved specs and lineage in manifests. Index common scalar fields for queries; retain nested transformations and priors in typed JSON rather than building a generic column for every nested value.
- Treat each systematic as a named study method with explicit dependencies and outputs. Reuse the same artifact/task contracts without forcing LSS, background-normalization studies, binning scans, smoothing scans, and toy batches into one algorithm.
- Keep model keys explicit in all model-dependent task records. Model names remain labels/associations, never the sole artifact identity.

## 8. Open decisions

- **Staging coverage:** confirm whether the current broad preselection is adequate for all supported analysis selections and shifts; record its resolved contract in staging identity. Current staging config uses `broad_preselection`, not the older `core/qroot.get_roi_cuts` path.
- **Blinding boundary:** blinding currently creates a separate staged ROOT via `create_blinded_rdf`, but fit-input writing/loading has additional blinded-data behavior. Specify whether blinding is a staging operation, an observation-construction operation, or both, and key each resulting product accordingly.
- **Fit vs input smoothing boundary:** current code smooths in fit construction. If smoothing/template assembly becomes a persisted reusable input product, move its identity boundary with that operation and update the contract.
- **Reuse policy:** specify which run-profile fields and diagnostics (draws, chains, R-hat/ESS thresholds) make an existing trace reusable for a request, and whether an insufficient trace is replaced or versioned under the same `fit_problem_id`.
- **Prior provenance:** specify whether prior generation is a separately identified artifact (with JAGS source, sample selection/seed, and method) or is resolved and embedded in `FitParams` with an explicit source reference.
- **Template variation persistence:** decide which transformed/rebuilt LSS and normalization-sampled templates are persisted as artifacts versus materialized transiently inside a systematic task. Either way, retain their complete derivation specification and sample identity.
- **Data-root granularity:** confirm that one data-root is keyed by (`blinded`, staging coverage, staging schema) and that nothing else needs to partition the store.
- **`skip_roi` placement:** currently an input parameter with a protected-comparison flag. Alternative: promote it to the data context if robustness runs should never share an `input` folder with nominal ones.
- **Output kind naming:** confirm the folder/manifest `kind` values (`fit`, `syst`, `toys`, `plots`) and whether `role` (`nominal`, `variation:<name>`, `toy:<i>`) is a manifest field only.
- **Group names:** confirm the advisory groups (`structure`, `selection`, `transforms`, `priors`) in section 3.
