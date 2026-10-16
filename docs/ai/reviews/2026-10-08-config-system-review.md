# QPyM2 config system — review & quality analysis

**Review date:** 2026-10-08  
**Scope:** Configuration docs and YAML, `qpym2.cfg`, pipeline runners, tests, and real-world usage in `notebooks/old/`.

## Executive summary

QPyM2 has a useful configuration foundation: Hydra composes environment, analysis, and hyperparameter YAML, while builder functions adapt the resulting `DictConfig` into the plain dictionaries expected by existing pipeline code. The model geometry helper also centralizes nontrivial rotated-axis derivations.

The main quality concern is that this configuration boundary is not yet consistently enforced. Several default and fallback paths can fail or silently use the wrong config; output-path precedence is inconsistent; grid parameters have multiple sources of truth; and notebooks mix helpers, direct attribute access, hardcoded paths, and obsolete schemas. Tests currently provide little coverage of real config composition.

## How configuration currently flows

1. `qpym2.cfg.get_config(...)` is a convenience wrapper around `HydraConfigManager.get_config` in `qpym2/cfg/hydra_adapter.py`. It uses Hydra's programmatic `compose()` API, not `@hydra.main`.
2. `config/config.yaml` declares four defaults: `environments`, `analysis`, `hyperparams`, and `toys`. Environment YAML holds infrastructure and paths; analysis YAML holds models, cuts, and components; hyperparameters hold systematic selections; toy YAML is also exposed through the separate `load_toy_config()` helper.
3. `qpym2/cfg/config_builders.py` adapts the composed `DictConfig` into staging/input dictionaries and provides model/component/path helpers.
4. Runners usually parse `--environment` and `--analysis`, load the config once, then pass derived dictionaries or a `TempConfig` histogram model into lower-level code. Many runners use `analysis = args.analysis or args.environment`.
5. Grid construction is not driven solely by YAML: `qpym2/cfg/params_grid.py` provides the `ALL_*` axes consumed by `scripts/gen_params.py` and runners. `config/hyperparams/base_params.yaml` also contains list-valued grid settings, but no Python consumer was found.
6. In `notebooks/old/`, `get_config()` is common (82 notebooks in the explored tree). Notebooks often use config builders for paths and staging, but also directly access nested config, manually construct `TempConfig`, hardcode paths, or retain imports for removed modules.

## Strengths

- The intended separation of infrastructure paths and physics definitions is clear in the docs and YAML layout.
- `get_h2uv_model()` centralizes u/v geometry and guards against deprecated `binning.u_range`/`v_range` settings.
- The builder layer is a practical migration boundary between Hydra and older plain-dictionary APIs.
- `params_grid.ALL_*` makes complete, ordered grid axes explicit and index-stable.
- The current config guide describes the live Hydra layout rather than the previously aspirational `extends:` scheme.

## Findings

### Correctness and API issues

1. **Duplicate adapter methods.** `HydraConfigManager.list_options()` and `get_default_selections()` are each defined twice in `qpym2/cfg/hydra_adapter.py` (first at lines 150/170, again at 228/252). Python uses the latter definitions, leaving the earlier implementations dead and making it unclear which behavior is intended.

2. **The master default is not usable by model-dependent helpers.** `config/config.yaml` selects `analysis: models_base` by default, but `config/analysis/models_base.yaml` defines `base_model`, `components`, and `cuts`, not `models`. Consequently, `get_config()` with no arguments returns a config for which `get_h2uv_model()` and `build_staging_config()` fail because no model is defined. The default-selection test currently confirms `models_base` as the default without testing that downstream helpers can use it.

3. **Wrong fallback model name.** In `qpym2/cfg/config_builders.py:415`, `get_h2uv_model()` falls back to `f'hmu_{model_key}'`, even though it resolves a default key earlier as `resolved_key`. If the model has no explicit `name` and no `model_key` argument was passed, the generated name is `hmu_None`.

4. **Smoothing comes from two configuration locations.** `get_h2uv_model()` reads the enabled flag from `cfg.hyperparams.smooth`, but kernel and iteration values from `cfg.base_model.smoothing` (`config_builders.py:397-409`). These sources can disagree; the behavior should have a single documented precedence or a single source.

5. **Explicit `plots_dir` may be ignored.** `get_plots_dir()` prioritizes `data_dir/plots` over the explicit `plots_dir` (`config_builders.py:498-512`). Environments such as `m2v5-1.yaml` set both `data_dir` and `plots_dir`, so the helper returns a derived path even if the explicit path differs. Separately, `hydra_adapter._apply_output_paths_defaults()` derives `plots_dir` from `fit_dir`, giving the field multiple defaulting policies.

6. **A runner shadows the loaded config variable.** `run_create_input.py:166` uses `cfg` as the loop variable for model configuration dictionaries, replacing the `DictConfig` loaded earlier in `main()`. It may not fail in the current code path, but it makes subsequent config access error-prone and obscures the distinction between the two config types.

7. **`model_cfg()` can silently load unrelated defaults.** In `qpym2/fit/params.py:89-102`, missing cut dictionaries trigger a parameterless `get_config()`. That loads master defaults rather than the environment/analysis used by the caller. If those dictionaries were omitted, cut lookup may use the wrong definitions or fail; requiring explicit cut dictionaries would make the dependency visible.

8. **A deprecation notice points to a nonexistent helper.** `get_data_names()` says to use `get_data_name` (`config_builders.py:242-244`), but no `get_data_name` implementation exists in the package or scripts.

### Structure and maintainability

9. **Hydra global state is cleared on each load.** Both `get_config()` (`hydra_adapter.py:130-138`) and `load_toy_config()` (`toy_config.py:63-74`) call `GlobalHydra.instance().clear()` before composing. This makes config loading dependent on process-global mutable state and can interfere with concurrent or interleaved uses in one process. It also repeats composition for repeated loads.

10. **There are multiple grid definitions.** Grid axes live in `params_grid.py`, CLI defaults, and `hyperparams/base_params.yaml`. The latter duplicates grid-like information using list values, while `hyperparams/base.yaml` defines scalar selections; no code reference to `base_params.yaml` was found. This is an easy source of confusion and drift.

11. **Environment and analysis selection are coupled by naming convention.** Runners commonly default `analysis` to the environment name. This is convenient for paired names such as `tech_note`, but environment and analysis are separate tiers and can be independently composed. The convention is not centrally validated, and a mismatched or missing pair can fail late during composition or use.

12. **Hyperparameter configuration is not consistently applied by runners.** Some runners expose `--hyperparams`, while others such as `run_single_fit.py` accept individual systematic parameters but do not accept a hyperparameter config name. This leaves the effective source of systematic choices dependent on the entry point.

13. **Optional output keys are inconsistent across environments.** For example, `systematic_dir` is set by `tech_note.yaml`, but not by `m2v5-1.yaml`; `data_dir` is the intended derivation base in some configs. Direct attribute access to a missing key raises an OmegaConf struct error. Helpers such as `get_syst_dir()` provide clearer behavior, but notebook consumers do not consistently use them.

14. **Some master Hydra settings are inert in the current loading path.** `hydra.job.chdir` and `hydra.run.dir` are declared in `config/config.yaml`, but config is loaded through `compose()` rather than Hydra's decorated application lifecycle. Their presence suggests behavior that this loader does not provide. Other fields such as `pipeline` and some file-pattern keys also appear only partially consumed.

### Documentation and schema drift

15. **The package docstring says three tiers.** `qpym2/cfg/__init__.py:5-8` documents three tiers, while the current master config and config documentation define four.

16. **Notebooks retain references to removed modules and schemas.** Examples in `notebooks/old/` include imports from `qpym2.cfg.input`, `qpym2.cfg.base_tech_note`, `qpym2.cfg.staging_1ton`, `qpym2.cfg.vars3`, and `qpym2.run.run_syst_new`; these modules are absent from the current tree. Some notebooks also probe old top-level layouts such as `cfg.environment` or `cfg.output_paths` instead of `cfg.infrastructure...`.

17. **Notebook configuration access is only partly migrated.** The notebooks use `get_config()` widely, but still construct `TempConfig` directly, read deep config paths rather than using builders, and contain hardcoded paths. These are practical signs that helper adoption and schema migration are incomplete.

### Real-world notebook friction

18. **Repeated setup boilerplate and import variation.** Notebook headers commonly repeat autoreload, sys.path/site setup, plotting settings, random seeds, and environment/analysis selection. `get_config` is imported through both `qpym2.cfg` and `qpym2.cfg.hydra_adapter`; helper imports are split between the package and `config_builders`.

19. **Paths diverge between config and manual overrides.** Examples include toy-results directories hardcoded under `data/toy_results` while helper-based config resolves `data_2syst/toy_results`, and notebooks manually setting fit paths instead of calling `get_fit_dir()`. Such overrides make it difficult to know which output tree a notebook is reading or writing.

20. **Notebook execution depends on hidden kernel state.** At least one notebook (`notebooks/old/plots/04_validation_plots.ipynb`) references `cfg` and helper names without defining/importing them in the notebook, with captured `NameError` output. Other cells use clone-specific `sys.path` additions and captured errors from removed modules.

### Test coverage

21. **Current tests do not validate real config composition.** `qpym2/tests/test_hydra_adapter.py` checks option listing and default selections, but not composing a real environment/analysis pair, building a model, resolving path precedence, or handling absent optional keys. The working tree also shows deletions under `qpym2/cfg/tests/` (fixtures and tests); those deletions are existing working-tree changes and are not assessed here as intentional or accidental.

## Severity overview

| Priority | Findings | Impact |
|---|---|---|
| High | Default analysis lacks `models`; `model_cfg()` silently loads unrelated defaults; lack of composition coverage | Default/helper failures or silently incorrect configuration can reach analysis code. |
| Medium | Path precedence; split smoothing source; repeated Hydra global-state clearing; duplicated grid definitions; inconsistent optional keys and runner interfaces | Behavior depends on entry point/config combination and is harder to reason about. |
| Low | Duplicate methods; fallback model name; variable shadowing; stale deprecation/docstrings | Confusion, latent errors, and maintenance overhead. |

## Recommended order of work

1. **Address correctness first:** decide on a usable master analysis default; make cut dictionaries explicit rather than silently reloading defaults; resolve `plots_dir` precedence; unify smoothing configuration; fix the model-name fallback.
2. **Reduce ambiguity:** remove duplicate adapter methods, repair the `get_data_names()` deprecation, rename the runner loop variable, and settle whether grid definitions belong in `params_grid.py` or YAML.
3. **Add config-focused tests:** cover representative environment/analysis composition, a successful model build, path precedence, optional path behavior, and invalid selection errors.
4. **Migrate notebooks incrementally:** standardize imports and a small setup cell; replace hardcoded result paths and obsolete imports; avoid relying on kernel state. Retain legacy notebooks as historical artifacts where useful, but distinguish them clearly from executable current examples.
5. **Clarify lifecycle and CLI semantics:** either remove inert Hydra settings or document their purpose under programmatic composition; validate environment/analysis selections and make hyperparameter handling consistent across runners.

## Review limitations

`notebooks/old/` is a local working directory and is gitignored per its `AGENTS.md`; its notebook outputs and hardcoded paths describe real usage but may reflect older environments or captured failures. The repository also has pre-existing uncommitted changes, including changes to config files and deletions under `qpym2/cfg/tests/`; this review describes the current working tree and does not attribute those changes to this review.
