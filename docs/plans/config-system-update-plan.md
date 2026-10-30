# Configuration System Update Plan

**Status:** Phase 1 & 2 Completed (v0.6 4-tier layout and documentation implemented)  
**Scope:** Foundation schema, programmatic adapters (`*_new.py`), tests, and documentation.

## Goal

Rework the configuration setup toward the clearer v0.6 structure, while doing
the minimum useful setup first. Establish a small, composable configuration
foundation now; document the remaining migration and quality work afterward.

## Documentation findings that shape this plan

- The v0.6 guide and the existing config-system review support the initial
  four-group split: environments, staging, analysis conventions, and named
  models.
- `docs/todos/analysis-workflow-and-config-reorganization.md` describes a
  broader eventual design: fit profiles, toy profiles, systematic methods,
  and study/campaign definitions. Those are not prerequisites for the basic
  configuration setup.
- The workflow and model-parameter docs explicitly distinguish model
  definitions from model points and fit-time choices. Avoid putting every
  setting into a generic `hyperparams` bag: input-affecting choices and
  fit-only choices have different owners and dependencies.
- The storage reviews discuss artifact IDs, layouts, and study records. This
  plan does not choose or implement a storage layout, identity scheme, or
  campaign planner. Config should provide data roots and settings without
  making storage paths into scientific identities.
- The current config README, root `AGENTS.md`, nested `config/AGENTS.md`, setup
  guide, and configuration review describe different stages of the design.
  Phase 2 must reconcile those summaries rather than add another competing
  guide.

The guiding separation is:

1. **Environment** — machine-specific roots, input locations, and output paths.
2. **Staging** — the shared raw-to-staged data contract.
3. **Analysis** — shared physics conventions such as cuts and components.
4. **Models** — named fit models, including ROI and binning.

These are configuration boundaries, not a requirement to migrate every runner
or notebook in the first pass.

## Phase 1 — Basic configuration setup

Keep this phase deliberately small and independently reviewable.

1. **Agree on the minimal schema and defaults.**
   - Use the four groups above, following the v0.6 organization.
   - Provide a composable baseline with a named model. A no-argument load
     should be valid for inspection/model construction, but must not imply that
     data paths are usable on every machine.
   - Keep machine-specific paths out of analysis/model definitions.
   - Keep staging selections independent of model ROI and fit variations.
   - Keep model selection explicit for production work; do not rely on
     whichever model happens to be first in a config mapping.
2. **Establish the master composition.**
   - Make `config/config.yaml` compose an environment baseline, staging
     contract, analysis baseline, and representative model set.
   - Avoid making a developer-specific absolute data path the portable default.
     A deployment profile or `QPYM2_DATA_DIR` can supply the real data root.
   - Preserve existing config selections/API where needed for current callers;
     compatibility removal belongs to a later migration, not this setup pass.
3. **Establish predictable path handling.**
   - Choose one data root as the base for relative data/output paths.
   - Define that absolute paths pass through and relative paths resolve from
     that root.
   - Do not add derived-path policies that silently override explicit paths.
   - Test path resolution with a temporary root; config composition itself
     should not require NERSC or analysis data to exist.
4. **Keep the first code touch limited to the config boundary.**
   - Update only config YAML and the minimum loader/builder code needed to
     compose the basic setup and resolve its paths.
   - Treat staging as a reusable config contract only. Stage manifests,
     artifact identities, storage layout, and format changes are separate work.
   - Do not migrate runner workflows, fit internals, notebooks, or storage
     behavior in this phase.
5. **Add narrow verification for the setup.**
   - Verify default composition and an alternate environment/model selection
     without requiring ROOT or production data.
   - Verify an explicitly selected model can be built by the model helper.
   - Verify relative paths, explicit absolute paths, and missing data-root
     behavior as specified.
   - Verify an invalid group selection produces a clear error.

### Phase 1 completion criteria

- The default config composes successfully and contains a usable named model.
- Environment paths and physics/model settings are separated.
- Staging configuration does not depend on a model ROI.
- Model-dependent production configuration can name its model explicitly.
- The basic path rules are consistent and covered by tests, while loading the
  config does not require production data to be mounted.
- Existing uncommitted work is preserved; no unrelated files or behavior are
  changed.

## Phase 2 — Document the resulting setup

After Phase 1 is agreed and implemented, update the existing documentation to
describe the actual YAML and loader behavior. Keep one canonical detailed
configuration guide and make other docs point to it. At minimum review:

- `config/README.md` and the current config guide (or decide which should be
  canonical if the guide is absent/obsolete).
- `config/AGENTS.md` and the configuration section of root `AGENTS.md`.
- `README.md` and `docs/setup.md` examples that show config selection or
  pipeline commands.
- `docs/ai/reviews/2026-10-08-config-system-review.md` status/recommendations
  if it is still presented as current guidance.

The guide should cover:

- The four groups, their ownership, and composition order.
- How to select an environment, analysis, and model, including the distinction
  between a model definition and a model point.
- The data-root and path-resolution rules, including precedence for explicit
  paths.
- Which builders are the supported boundary between `DictConfig` and legacy
  runner dictionaries.
- A short example of composing config and explicitly selecting/building a
  model, with a separate example for supplying a local/NERSC data root.
- Any deliberate defaults or limitations of the initial setup.

Keep the broader target in `docs/todos/analysis-workflow-and-config-reorganization.md`:
document fit profiles, toys, systematic methods, campaigns, and planner work
there as future phases rather than implying the basic four groups already
implement those concepts. Documentation should describe implemented YAML and
loader behavior, not intended behavior as if it already exists.

## Phase 3 — Follow-up cleanup and migration

Track these as separate follow-up work, not prerequisites for the basic setup:

- Decide the supported home and precedence for smoothing and other systematic
  selections. Put fit-time behavior in a fit-profile concept when that phase is
  undertaken; do not silently default stale/unsupported config fields.
- Decide how existing `hyperparams` and `toys` consumers transition, and remove
  old API options only after their callers have migrated.
- Resolve remaining path-helper precedence and optional-path behavior.
- Remove stale fallbacks that reload unrelated default configuration.
- Establish the relationship between Python `params_grid.py` axes, CLI
  defaults, and any future study/campaign definitions.
- Review Hydra global-state lifecycle and inert Hydra settings.
- Migrate runners incrementally to explicit environment/model/model-point and
  fit-profile selection through shared builders.
- Migrate or clearly label legacy notebooks; replace hardcoded paths and
  obsolete schemas over time.
- Expand tests to cover composition, builder outputs, path behavior, and
  invalid selections; later add a small synthetic multi-model pipeline smoke
  test as recommended by the workflow review.

## Out of scope for the initial update

- Rewriting all runner interfaces or the fit API.
- Bulk-editing historical notebooks.
- Building fit profiles, campaign/study planning, or complete toy/systematic
  schemas.
- Redesigning parameter grids, artifact identity/storage, or the proposed
  `fit_data/` layout.
- Changing staged-data formats or introducing Parquet.
- Unrelated code or working-tree changes.

## Review checkpoints

1. Approve the four-group starting schema, data-root behavior, and compatibility
   boundary before implementation.
2. Review Phase 1 as a focused config setup change, including its tests.
3. Reconcile existing docs after the basic setup is settled; avoid duplicating
   the broader campaign/storage roadmaps.
4. Prioritize later schema and runner migrations separately based on
   correctness impact and user friction.
