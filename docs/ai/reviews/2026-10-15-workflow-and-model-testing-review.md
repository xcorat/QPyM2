# Scientific workflow and multi-model testing review

**Review date:** 2026-10-15  
**Scope:** Multi-step analysis pipelines, parameter/model studies, testing patterns, provenance, and HPC workflow tools; compared with the current QPyM2 working tree and stated need for fast, dynamic, modular experimentation.  
**Interpretation:** “Multi-model” means evaluating multiple QPyM2 physics/model configurations (for example ROI, threshold, bin size, or component definitions) and parameter points through shared pipeline stages. This review does not assume that QPyM2 is doing generic machine-learning model selection.

## Executive summary

QPyM2 already has useful foundations for an adaptable analysis system: independent pipeline runners, importable library functions, Hydra configuration groups, explicit model definitions, a central parameter-axis module, local multiprocessing, TSV/SLURM execution, and some stage-level caching. These are worth preserving.

The main gap is that the *campaign*—the declaration of which models, parameter points, variations, and downstream products should be produced—is not yet a first-class object. The current runners and grids encode pieces of the plan separately. For example, `run_create_input.py` selects a single model through `get_h2uv_model(cfg)` and hard-codes the default cut/shift slices, while `run_fits.py` task rows include cut and shift values but not a model key. Consequently, a multi-model study is not consistently represented as one validated, inspectable, resumable run plan.

The strongest fit to the requirements is an incremental, **QPyM2-native campaign planner and artifact graph** that composes the existing Hydra configs and emits the existing TSV/SLURM jobs. Borrow the dynamic task construction and cache discipline of Pydra/Parsl, file-contract and dry-run discipline of Snakemake/CWL, and the campaign-spec ideas in Maestro/REANA. Keep ROOT/PyMC analysis code as domain libraries and treat any external workflow engine as an optional execution adapter, not as the source of physics semantics.

This is especially compatible with QPyM2’s own workflow reorganization note, which already proposes named models/model points/studies, explicit dependencies, immutable artifact identities, manifests, dry runs, and incremental migration while retaining TSV and SLURM. The external research supports that direction; it does not suggest a need to replace the current pipeline wholesale.

## Requirements inferred from QPyM2

The design should optimize for these concrete requirements:

1. **Fast experimentation:** add a new model or systematic method without editing a large central runner or duplicating an entire pipeline script.
2. **Dynamic plans:** a study may have a simple Cartesian sweep, a filtered set of points, a nominal/reference point, a specialized systematic algorithm, or downstream tasks conditional on upstream artifacts.
3. **Correct reuse:** re-use staging and fit inputs when their *scientific dependencies* are unchanged; rebuild when a relevant selection, transformation, input dataset, or implementation changes.
4. **Multiple configurations as first-class inputs:** model keys and fit/toy/systematic profiles should travel with each task and result, rather than being implied by whichever model a helper selects by default.
5. **HPC-friendly execution:** retain local dry runs, SLURM arrays/packed jobs, and the option to run a small matrix on a workstation. Avoid requiring a persistent orchestration service for ordinary work.
6. **Testability without NERSC data:** test configuration resolution, task planning, artifact contracts, and a representative small end-to-end pipeline using synthetic/fixture data; leave full-data checks as a separate integration tier.
7. **Auditability:** for any result, recover the resolved parameters, input identities, code/environment, upstream artifact IDs, random seed, and execution outcome.

The priorities are not “use the newest workflow engine” or “run every model combination all the time.” Dynamicity should come from a declarative study plus small reusable analysis functions, and test breadth should be proportional to cost and scientific risk.

## External framework patterns

| Framework / practice | Useful pattern | Fit for QPyM2 | Limits / caveats |
|---|---|---|---|
| **Snakemake** | File-oriented rules infer a DAG from declared inputs/outputs; configuration, modules, profiles, dry runs, DAG visualization, resource settings, caching, and reports. Can generate small per-rule tests from a successful run. | Strong conceptual match for staging → input → fit → toys/systematics → summaries, especially when artifacts are files and independent points can run concurrently. Its modules and config profiles illustrate code/workflow reuse separately from cluster settings. | Its generic file DAG does not understand QPyM2 physics compatibility. Must still encode precise artifact keys and model-point semantics. NERSC recommends Snakemake primarily for single-node use on Perlmutter rather than submitting every workflow rule as a Slurm job; avoid overwhelming the scheduler with fine-grained submissions. |
| **Nextflow** | Composable named workflows and processes connected by dataflow channels; supports scatter/gather, conditional paths, modules, cached tasks, resume, and per-task trace/report metadata. | Useful pattern for modular workflows and fan-out over many models/points followed by aggregation. Its “process as a function with explicit inputs/outputs” is a good contract to emulate. | Groovy/JVM, separate runtime model, file staging and work-directory management add operational complexity. NERSC documents scheduler/polling and filesystem caveats; better to pilot only if QPyM2 needs its execution features rather than adopting it for a naming convention. |
| **Parsl** | Python functions become asynchronous apps; ordinary Python constructs a dependency graph dynamically from futures and files. Executors separate workflow code from local, Slurm, or high-throughput execution. | Very relevant for a Python-centric, dynamic planner and for fan-out/fan-in across model points. It can run on NERSC and can express workflows that do not fit one static parameter grid. | Distributed apps require explicit inputs and robust serialization; a long-lived driver/orchestrator and Slurm configuration need care. Caching/checkpointing must be designed around deterministic task inputs and should not substitute for QPyM2’s domain-level artifact manifest. |
| **Pydra (Nipype ecosystem)** | Python tasks and nested workflows share one abstraction; split/combine (map/reduce), conditional dynamic construction, type checking, per-task working directories, input-hash caching, and multiple workers including SLURM. | Closest design reference for “Python functions as modular steps” with fast local debugging and model-point fan-out. Good examples of composition, hashing, and testable typed interfaces. | The current documentation/site exposes 1.0 alpha versions and recent versions target Python 3.11+, so check stability and compatibility with QPyM2’s supported Python environment before any adoption. QPyM2 can borrow the design without taking a runtime dependency. |
| **Maestro (LLNL; NERSC-supported)** | Human-readable YAML study specs declare ordered steps, dependencies, Slurm resources, and parameter sweeps; status is tracked per step/parameter instance. | Particularly close to the *campaign* concept and useful as a reference for readable study definitions that can expand to task matrices. NERSC documents native Slurm integration. | Requires the `maestro run` driver to remain alive and has polling delay. A generic study spec still needs domain-specific artifact compatibility and fit/model semantics. Could be evaluated as an execution backend later, not as the campaign data model by default. |
| **REANA / CERN Analysis Preservation** | Separates analysis inputs, code, software environment, workflow steps, and outputs; runs containerized workflows and is aimed at preserving/reusing particle-physics analyses. | Strong HEP precedent for recording enough information to rerun a published analysis, and reinforces that configuration, code, environment, input data, and execution steps are all part of the analysis record. | REANA is a platform that needs infrastructure and containerized execution. QPyM2’s Perlmutter/data constraints make its preservation principles more immediately useful than deploying REANA itself. |
| **ROOT RDataFrame** | Lazy, composable filter/define/action graph; graph visualization, cut-flow reports, implicit parallelism, and `Vary`/variation extraction for systematic columns. | A good model for modular event-level transformations and for computing multiple related histograms in one event loop. It aligns with QPyM2’s ROOT-based raw/staging path and systematic variations. | RDataFrame is an event-processing layer, not a multi-stage campaign orchestrator. It should not own fit/toy task tracking. Its ROOT expressions also reinforce the value of isolating format-specific operations behind a staging I/O contract. |
| **Coffea (HEP Python)** | Separates an analysis processor from execution; columnar event operations with Awkward/NumPy and executors ranging from local to distributed systems. | Useful HEP precedent for making analysis transformations reusable and execution-independent, and for testing the analysis code locally before scaling. It supports QPyM2’s direction to reduce ROOT coupling downstream. | Coffea targets collider-style columnar event processing, not QPyM2’s PyMC fit orchestration or NERSC artifact lifecycle. Wholesale adoption would not solve campaigns, fit provenance, or systematics dependencies. |

### What the frameworks agree on

- **Represent dependencies explicitly.** A step consumes identified inputs and produces named outputs. Independent work can run in parallel; downstream work starts only when dependencies are satisfied.
- **Separate scientific logic from execution.** The same operation should be callable directly in tests/local mode and by a pool, task array, or scheduler adapter.
- **Make parameter sweeps data, not handwritten loops.** A point has an identity and resolved values; expansion, filtering, dry-run counts, and scheduler-row generation happen centrally.
- **Resume from verified outputs, not mere path existence.** Cache/reuse decisions need task code, parameters, environment, and input/output checks, with explicit force/rebuild behavior.
- **Keep intermediate products inspectable.** Workflow systems use files and metadata to make step boundaries durable. In QPyM2, ROOT/HDF5/NetCDF outputs can remain the immediate artifacts while their schemas and lineage become explicit.
- **Allow special workflows.** Real analyses include conditional steps, reference fits, and specialized systematics; forcing all studies into a single generic Cartesian product is a poor abstraction.

## Standards and recommendations

These standards are complementary; QPyM2 does not need to implement all of them to benefit from their concepts.

### FAIR data principles

FAIR means Findable, Accessible, Interoperable, and Reusable. Especially applicable principles are persistent identifiers, rich machine-readable metadata, explicit links to the data being described, interoperability through shared formats, and detailed provenance. For QPyM2 this translates into a stable ID and sidecar manifest for every campaign/artifact, with each result pointing to its input artifacts and resolved configuration. FAIR is a guiding framework, not a requirement to publish CUORE data openly; access restrictions and controlled paths remain compatible with FAIR metadata.

### W3C PROV

PROV-O provides a general provenance vocabulary: **Entities** (input data, configs, artifacts), **Activities** (staging, input building, fit, aggregation), and **Agents** (software, analysis team, execution service). A QPyM2 manifest can implement this simple structure without RDF. The useful invariant is that an output artifact is generated by a named activity that used known upstream entities and was associated with a particular code/software agent.

### CWL and RO-Crate

CWL is a vendor-neutral specification for tools and workflows with declared inputs/outputs and DAG dependencies. It validates input/output schemas and allows independent steps to run concurrently; it also has scatter/gather for repeated work. RO-Crate and Workflow Run RO-Crate package workflow definitions, inputs, outputs, code, and run provenance in interoperable metadata. These are valuable **export/interoperability targets** once QPyM2’s internal manifest settles. Adopting CWL immediately could add descriptor boilerplate and constrain Python-native dynamic study algorithms; start by making QPyM2 artifact/task contracts CWL-like.

### Testing practices

pytest supports parameterized inputs and readable case IDs, which naturally suits a small suite of selected model/profile combinations. Scientific Python’s development guide recommends separate fast and slow suites, user-visible behavior tests, integration tests that exercise component interactions, and unit tests for focused logic. Snakemake’s generated-rule tests make a further useful point: pipeline step tests should use small data and compare meaningful outputs rather than try to store or rerun production-sized data.

For numerical science, use several complementary oracles:

- **Analytic/invariant tests:** event counts are conserved where expected; arrays/histograms have expected shapes and finite values; probability priors are normalized; IDs are stable under ordering changes.
- **Differential tests:** compare ROOT and alternate staging backends, or new and reference transformations, on a tiny representative sample with explicit tolerances.
- **Golden/reference tests:** store compact summary quantities (selected counts, histogram bins, fit diagnostics) and update them deliberately when the scientific result intentionally changes.
- **Parameterized model tests:** validate each declared model’s config and build path; run a deliberately small representative set of model/point combinations end to end.
- **Full-data integration tests:** keep separate and optional/controlled because they need real data, ROOT/PyMC resources, and configured NERSC paths.

## QPyM2: what already matches and what is missing

### Good foundations to retain

- `qpym2.run` has stage runners and importable functions, giving a natural seam between science code and command-line/HPC execution.
- Hydra already separates environments from analysis/model definitions, while `config_builders.py` centralizes much of the compatibility with existing plain-dict APIs.
- `params_grid.py` makes several parameter axes explicit and ordered.
- The pipeline is already broken into durable artifact stages: staged ROOT, input HDF5, fit NetCDF/HDF5, toy/systematic results, and summaries.
- Local pools and TSV-driven SLURM arrays are familiar, debuggable execution mechanisms. The workflow refactor note correctly avoids requiring their immediate replacement.
- The input runner caches repeated reads across Co-shifts within a task, showing that reuse is already valuable at sub-task level.
- Tests are marked so the default suite can run without NERSC data; pure transforms and configuration adapter code already have test coverage.

### Main friction against dynamic multi-model work

1. **Model selection is not an explicit grid axis in the shared task schema.** `run_create_input.py` resolves one `h2uv_model` before constructing its cut/cutoff/shift xarray grid. `run_fits.py` reconstructs a model config from environment plus row parameters, but the TSV job schema does not include the model key or analysis/profile identity. Model variants can therefore require separate/manual runner invocations and risk inconsistent model selection between input and fit stages.
2. **Several defaults are code slices rather than named study definitions.** The input runner selects fixed slices of `ALL_CO_CUTS`, `ALL_NDBD_CUTS`, and `ALL_CO_SHIFTS`; a developer must edit Python to change that particular default experiment. This is fast for one author but fragile for review, reuse, and reproducibility.
3. **Not every variation has the same dependency graph.** A smoothing setting can be fit-time only; a Co-shift may change input construction; LSS variations can start from staged events; bias toys need fitted artifacts and input edges. A global Cartesian grid either rebuilds too much or hides meaningful dependencies.
4. **A path does not prove an artifact is compatible.** The workflow TODO already notes ROI-related staging compatibility and recommends immutable staging IDs/manifests. Without upstream IDs/schema/config identity, silent reuse or overwrite is difficult to distinguish from a valid cache hit.
5. **Test tiers do not yet verify a small complete campaign.** Current tests include pure functions and an integration test requiring the configured environment/data, but the desired middle tier—a small synthetic multi-model, multi-stage smoke campaign—is not yet a standard contract.

The preceding config-system review also identifies related needs: real config composition tests, explicit model keys, consistent configuration precedence, and fewer hidden defaults. Campaign planning would help bring these fixes together, but it should not hide unresolved config errors behind another layer.

## Recommended target design

### 1. Name the domain concepts

Keep these concepts separate in configuration and APIs:

- **Analysis:** shared physics conventions, cuts, components, constants.
- **Model:** one model definition (ROI, threshold, binning, components/overrides).
- **Model point:** a fully resolved input/fit choice, such as model key + Co cut + NDBD cutoff + Co shift + any input-affecting variations.
- **Fit profile:** sampler, chains/draws/tune, prior/smoothing settings and fit-only options.
- **Toy profile:** toy type/count/injection distribution/seed policy.
- **Systematic method:** a named algorithm and its parameters, e.g. LSS, bin-size matrix, smoothing, or MC variation; it can declare its own reference and dependency requirements.
- **Study/campaign:** a user-authored request composing the above, choosing points or sweep axes, and selecting outputs/tasks.
- **Artifact:** one immutable stage, input, fit, toy batch, systematic result, or summary with an ID and manifest.

This grows naturally from the current Hydra groups. A model key should be explicit in any model-dependent runner/task; the first model should never be an implicit default for a production study.

### 2. Resolve a campaign to a small, inspectable plan

Build a Python planning API that composes the selected Hydra configuration once, validates it, expands only the requested points, and produces an explicit task graph. The graph builder should encode algorithm-specific dependencies (not assume every task is “one row of one giant grid”). It should support direct Python construction for notebook experimentation and a YAML study file for repeatability.

An illustrative—not committed—study could look like:

```yaml
study:
  name: threshold-and-bin-validation
  environment: m2v5-1
  analysis: tech_note
  models: [tech_note_200, tech_note_250]
  fit_profile: quick_validation
  points:
    co_cut: [gamma10_v150]
    ndbd_cutoff: [u1740]
    co_shift: [0.8, 1.0]
  tasks: [input, fit, posterior_checks]
  outputs: [fit_summary]
```

Specialized procedures should be able to declare a reference model point, non-Cartesian variations, and required artifacts rather than pretending every systematic is the same sweep. Study specifications select named methods/profiles; the actual transformations and statistics remain ordinary tested Python functions.

### 3. Keep an explicit artifact/task contract

Each task should declare:

- stable task kind and version/schema;
- typed/validated parameters and upstream artifact IDs;
- expected outputs, format/schema, and validation checks;
- resource needs (cores, memory, walltime) separately from scientific configuration;
- code identity and software environment identity;
- random seed policy where stochastic operations occur.

Each output artifact should be written to a versioned/unique location and have a manifest. A practical first manifest can be JSON and include `artifact_id`, `artifact_type`, `schema_version`, canonical resolved config or its full snapshot, config hash, upstream artifact IDs, raw input file list/identity, source revision and dirty-tree indicator, package/dependency versions, command/task ID, random seed, timestamps, runtime/status, output paths/checksums, event counts/schema, and error/log references.

Use layered identities rather than one hash for the whole campaign:

- **staging ID:** raw dataset identity + staging transform/selection/schema version + relevant implementation identity;
- **input ID:** staging ID + model key/point values affecting selection/histograms + input builder identity;
- **fit ID:** input ID + fit profile/prior identity + sampler/seed/software details where they affect results;
- **toy/systematic ID:** exact reference artifact IDs + variation method/profile + seed and design.

Do not hash irrelevant settings into upstream identity: a fit-only smoothing variation should not invalidate staging or inputs. Conversely, do not reuse a staging artifact when its actual selection excluded events needed downstream. For very large raw files, record the canonical source list and available file metadata/checksums under a documented policy; avoid an unplanned full-data hashing pass on every run.

Use immutable output IDs and explicit `force`/`rebuild` behavior. Check manifest/schema and required files before reuse, write output to a temporary location and finalize atomically where practical, and record failed attempts separately from successful artifacts. Preserve scientific result identities separately from execution-attempt IDs so a retry does not look like a new model.

### 4. Preserve existing execution while making plans portable

The planner should initially emit:

1. a readable dry-run summary and dependency graph;
2. a campaign manifest with task IDs and dependencies;
3. the existing TSV rows / `scripts/gen_params.py`-compatible command rows for task arrays;
4. a local execution plan for small validation matrices.

Existing `run_job.sh`, packed SLURM, local pools, and stage runners remain execution adapters. Later, a backend interface could support NERSC Slurm arrays, local multiprocessing, or a workflow manager. The workflow plan—not any one backend—remains the record of what the analysis requested. Keep resource settings in an execution profile (NERSC/local) and physics/sweep settings in the study/model configs.

### 5. Make extension points modular but bounded

New model variants should mostly be declarative config. New scientific operations should be tested library functions with a small, explicit input/output API. For systematic methods, use a registry/factory keyed by a stable method name (for example `bin_size_matrix`, `lss_dq`, `mc_shape`) whose implementation declares how it maps inputs to outputs and how its parameters expand.

This preserves fast code changes without a general-purpose plugin framework. Avoid dynamically importing arbitrary code paths or building a universal workflow DSL before several real studies prove the common contracts. Put reusable logic in `fit/`, `syst/`, `io/`, or a narrow workflow-planning module, not inside CLI parsing or notebook cells.

### 6. Build a testing ladder around a tiny campaign

Recommended tiers:

| Tier | Example coverage | Runtime/data |
|---|---|---|
| Unit / pure | model config validation, dependency classification, canonical hashing, sweep expansion, path builders, transform invariants | seconds; no ROOT or NERSC data |
| Component | build an input with synthetic event arrays; validate model construction; fit a tiny synthetic histogram with a deliberately small sampler; read/write artifact manifests | seconds to minutes; temp directory, no production data |
| Pipeline smoke | select 2 model configs × 1–2 points, plan and run stage → input → lightweight fit/check → summary; assert every output points to correct upstream IDs | small checked-in/generated fixture data; routine developer/CI target |
| Differential/scientific regression | ROOT vs Parquet staging; old/new histogram or systematic summaries; fixed-seed posterior diagnostics with justified tolerances | small controlled sample; run after relevant changes/nightly |
| Full integration / production validation | real configured dataset, complete fit/toy/systematic behavior, resource/runtime observations | Perlmutter/analysis data; opt-in or scheduled |

Use pytest parameter IDs containing model and point names, not opaque array indices. Test that config errors fail during plan resolution, before spending scheduler time. Ensure the fast test suite can import without PyROOT where possible. Commit compact fixtures/reference summaries rather than production-sized ROOT/HDF5 files. Do not compare Bayesian samples byte-for-byte: test deterministic preprocessing exactly where appropriate and statistical properties/diagnostics or stable summaries for stochastic fits.

## Prioritized migration plan

The repository’s existing `docs/todos/analysis-workflow-and-config-reorganization.md` is already a good implementation plan; this review suggests the following order and clarifies what should be delivered at each step.

1. **Close the identity and schema gap first.** Define schemas for staging and fit-input artifacts, validate their manifests, and implement explicit staging/input IDs. Confirm broad staging preselection or make staging compatibility checks authoritative. Use one existing analysis on a small sample.
2. **Add explicit model keys to planning/task records.** Model identity must flow through generated inputs, fits, outputs, TSV rows, and collected results. Add a planner-level test that catches input/fit model mismatch.
3. **Pilot a real multi-model study.** Migrate the existing bin-size × Co-shift matrix or a small two-threshold study. Include its nominal/reference fit and prior-reuse rules explicitly. Dry-run output should explain each task and artifact dependency before execution.
4. **Generate current TSVs from a campaign manifest.** Keep current submit scripts and pool runner. Add deterministic task IDs independent of row order; record seeds based on campaign/task identity; make reruns idempotent and incompatible artifacts visible.
5. **Add the small pipeline smoke tier.** Use synthetic/fixture inputs to exercise at least two model variants through planning and representative steps. Keep expensive PyMC sampling tiny or replace it in smoke tests with model construction/diagnostics where statistically appropriate.
6. **Migrate remaining study types incrementally.** Add toys and then specialized systematics by declaring their actual graph patterns. Do not force LSS/MC/smoothing/binning into one execution function just to make the API look uniform.
7. **Decouple staged-event I/O.** Implement a backend-neutral staging reader/writer and compare ROOT and Parquet on schema, selected events, event counts, histograms, and relevant systematic outputs. This improves test speed and modularity, but should follow or accompany stage identity—not precede it as an isolated format change.
8. **Export portable provenance when useful.** Keep a simple JSON manifest internally. If external reuse/preservation becomes a requirement, map it to W3C PROV/RO-Crate and consider CWL descriptors for stable command-line steps. Evaluate REANA or an external scheduler only with a concrete deployment need.

## Adoption guidance: what to borrow and what not to adopt wholesale

**Borrow now:**

- Pydra/Parsl: Python-native task functions, dynamic composition, mapping/combining, local debug execution, explicit inputs and outputs.
- Snakemake/CWL: dry-run graph, named file contracts, output validation, resource profiles separated from analysis parameters, reusable step modules, small pipeline tests.
- Maestro: a concise campaign/study specification with dependencies and named parameter sweeps.
- REANA/FAIR/PROV/RO-Crate: complete reproducibility record and provenance chain, kept as lightweight manifests initially.
- ROOT RDataFrame/Coffea: event-processing APIs should be composable and execution/storage boundaries should not be fused into the statistical fit layer.

**Do not adopt yet:**

- A full Nextflow/Snakemake/Parsl/REANA migration before a QPyM2 artifact/task contract and pilot campaign exist.
- A database, persistent workflow server, or generic plugin marketplace for a workflow that currently runs through files and Slurm arrays.
- One giant Cartesian grid for every model and systematic.
- Cache reuse based only on output filename/existence, or cache keys that omit code, model, input, fit profile, or seed dependencies.
- Moving every numerical transformation into YAML; YAML should describe choices, while reviewed/tested numerical algorithms remain Python.
- Production-sized end-to-end test data or exact byte-comparisons of stochastic fits.

## Conclusion

QPyM2’s requirements point to a **thin domain-specific planning/provenance layer**, not an immediate replacement workflow manager. A declarative campaign should make experiments easy to change, while Python task implementations stay modular and testable, and existing local/SLURM mechanisms execute the resolved plan. The most important enabling changes are explicit model identity, dependency-aware artifact IDs, small deterministic plan generation, and a meaningful synthetic multi-model smoke test. Once these contracts are stable, QPyM2 can export to or adopt an established engine selectively without losing scientific semantics or slowing down experimental code changes.

## Sources

External documentation checked 2026-10-15. Framework pages evolve; use the linked official documentation for current details.

### Workflow and HEP frameworks

- [Snakemake: workflows, rules, and DAGs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
- [Snakemake: modularization](https://snakemake.readthedocs.io/en/master/snakefiles/modularization.html)
- [Snakemake: profiles and CLI configuration](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
- [Snakemake: automatically generating unit tests](https://snakemake.readthedocs.io/en/stable/snakefiles/testing.html)
- [Snakemake: reporting and provenance](https://snakemake.readthedocs.io/en/v9.18.2/snakefiles/reporting.html)
- [Nextflow: workflows and dataflow composition](https://www.nextflow.io/docs/latest/workflow.html)
- [Nextflow: caching and resuming](https://www.nextflow.io/docs/latest/cache-and-resume.html)
- [Parsl: workflow dataflow](https://parsl.readthedocs.io/en/stable/userguide/workflows/workflow.html)
- [Parsl: memoization and checkpointing](https://parsl.readthedocs.io/en/latest/userguide/workflows/checkpoints.html)
- [Pydra: project documentation](https://nipype.github.io/pydra/)
- [Pydra: hashing and caching](https://nipype.github.io/pydra/explanation/hashing-caching.html)
- [NERSC: workflow tools overview](https://docs.nersc.gov/jobs/workflow/)
- [NERSC: Snakemake recommendations](https://docs.nersc.gov/jobs/workflow/snakemake/)
- [NERSC: Parsl guidance](https://docs.nersc.gov/jobs/workflow/parsl/)
- [NERSC: Maestro](https://docs.nersc.gov/jobs/workflow/maestro/)
- [NERSC: Nextflow caveats](https://docs.nersc.gov/jobs/workflow/nextflow/)
- [REANA specification reference](https://docs.reana.io/reference/reana-yaml/)
- [CERN Analysis Preservation: reusable analyses and REANA](https://analysispreservation.cern.ch/docs/general/workflows.html)
- [ROOT RDataFrame reference and user guide](https://root.cern/doc/master/classROOT_1_1RDataFrame.html)
- [Coffea documentation](https://coffea-hep.readthedocs.io/en/latest/)

### Standards and software testing

- [GO FAIR: FAIR Guiding Principles](https://www.gofair.foundation/fair-principles)
- [W3C Recommendation: PROV-O](https://www.w3.org/TR/prov-o/)
- [Common Workflow Language v1.2.1: Workflow specification](https://www.commonwl.org/v1.2/Workflow.html)
- [CWL: workflow features and scatter](https://www.commonwl.org/features/)
- [Workflow Run RO-Crate](https://www.researchobject.org/workflow-run-crate/)
- [pytest: parametrizing fixtures and test functions](https://docs.pytest.org/en/stable/how-to/parametrize.html)
- [Scientific Python Development Guide: testing recommendations](https://learn.scientific-python.org/development/principles/testing/)

### QPyM2 sources

- [`AGENTS.md`](../../../AGENTS.md): current pipeline, configuration, environment, and test conventions.
- [`config/config_guide.md`](../../../config/config_guide.md): current Hydra layout and model configuration.
- [`docs/todos/analysis-workflow-and-config-reorganization.md`](../../todos/analysis-workflow-and-config-reorganization.md): proposed campaign, artifact identity, planner, provenance, and migration sequence.
- [`docs/todos/staging-parquet.md`](../../todos/staging-parquet.md): staged-data schema and backend-decoupling proposal.
- [`qpym2/run/run_create_input.py`](../../../qpym2/run/run_create_input.py), [`qpym2/run/run_fits.py`](../../../qpym2/run/run_fits.py), [`qpym2/cfg/params_grid.py`](../../../qpym2/cfg/params_grid.py): current model/grid/task handling.
- [`2026-10-08-config-system-review.md`](2026-10-08-config-system-review.md): configuration review findings relevant to explicit model selection and testing.
