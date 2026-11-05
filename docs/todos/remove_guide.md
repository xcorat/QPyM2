I’ll inventory `qpym2/` against the new storage layer and the retirement scope in the supplied docs. I’ll only report recommendations—no files will be changed—and I’ll distinguish genuinely obsolete orchestration from reusable algorithms and storage primitives.

The key constraint is that storage is not yet a full replacement boundary: the current study point still calls `run_single_storage_fit`, while staging identity/manifests and some reuse validation are still pending in the plan. I’m treating runner removal as phased/conditional, not as safe immediate deletion.

## Recommendation

**Don’t remove the legacy pipeline wholesale yet.** The storage layer handles input and fit artifacts, but it is not yet a complete replacement boundary: `studies.point` still invokes `run_single_storage_fit`, and staging identity/manifests and reuse validation remain unfinished in the v0.6 plan. Treat most runner removals as **post-migration work**, not safe deletions today.

## Likely removable after a consumer check

These files have no import/call sites elsewhere in the repository scan, or duplicate a more actively used implementation. Check for external consumers before deleting:

- `qpym2/cfg/model_utils_new.py` — duplicate of `model_utils.py`; no consumers found.
- `qpym2/core/qroot.py` — no consumers found.
- `qpym2/core/hists.py` — no consumers found; overlaps with `qpym2/utils/hists.py`.
- `qpym2/utils/math.py` — appears to serve only `core/hists.py`; remove with that module.
- `qpym2/plotting/post2.py` — no consumers found; appears to be a prototype and uses the older fit-I/O interface.
- `qpym2/plotting/read_hists.py` — no consumers found; uses the older input/histogram construction path.

`qpym2/utils/grid_parser.py` is also unreferenced elsewhere, but it is a standalone CLI. Retire it only when the toy-grid/launcher workflow has been replaced by validated study definitions.

## Planned runner retirement

The following are migration/compatibility surfaces from the sundry TODO. Retire them only after the supported planner/API owns the behavior and parity is checked:

- **Legacy input and fit grid:** `qpym2/run/run_create_input.py`, `run_single_fit.py`, `run_fits.py`
- **Legacy toy orchestration:** `run_toy_fits.py`, `run_toy_batch.py`
- **Specialized fit matrices:** `run_bin_size_fits.py`, `run_bin_size_bias_fits.py`
- **v0.6 CLI wrappers:** `run_create_storage_inputs.py`, `run_single_storage_fit.py`, `run_study.py` — move their execution logic behind the supported API first; retain thin CLIs if they remain a supported interface.
- **Staging:** `run_create_staging.py` — do not retire until staging has an identity, manifest, publication/reuse validation, and planner integration.
- **Systematics:** `run_syst_lss.py`, `run_syst_mc.py`, `run_ndbd_shift_syst.py` — migrate orchestration only after their behavior and consumers are accounted for. Their algorithms are explicitly out of scope for the systematics redesign.
- `qpym2/run/__init__.py` — remove after the runners it catalogs are retired.

## Refactor or consolidate

- **Configuration boundary (COMPLETED):** Merged `cfg/hydra_adapter.py` with `cfg/hydra_adapter_new.py`, and `cfg/config_builders.py` with `cfg/config_builders_new.py`. Promoted modern 4-tier adapters to canonical module names; migrated all callers and removed duplicate `*_new.py` files (`hydra_adapter_new.py`, `config_builders_new.py`, `model_utils_new.py`).
- **Old configuration helpers:** Retire or replace `cfg/model_utils.py`, `cfg/params_grid.py`, and `cfg/toy_config.py` as old fit paths and TSV/toy launchers are migrated. Keep the scientific LSS definitions in `cfg/syst_config.py`, but move path resolution and defaults into the validated config/method boundary rather than deleting the scientific settings with the launcher.
- **Study planning and execution:** Refactor `studies/config.py`, `assembly.py`, `point.py`, and `study.py`. Separate validated study schema/planning, point specifications, task execution, and result collection. `point.py` currently calls a runner directly, and `study.py` mixes status checks, execution, toy algorithms, and collection.
- **Keep the storage/study boundary distinct:** `storage/study.py` handles persistence and indexing; `studies/` handles study definitions and execution. Coordinate their manifest responsibilities, but don't merge them into one undifferentiated module.
- **Fit/input orchestration:** Consider splitting identity/hash resolution from execution in `storage/inputs.py` and `storage/fits.py`. Both currently combine several responsibilities, and the fit assembly still resolves component config without an explicit model key.
- **Old versus new file formats:** `io/fit.py` is a legacy fit-file API and can be retired after its consumers migrate. **Do not delete `io/hdf.py` wholesale**: it still provides shared helpers used by toy/systematic I/O and supports older products. Keep a compatibility reader/writer until those consumers and files are handled.
- **Plotting:** Keep `plotting/legacy.py` for now; `posterior_plots.py` calls it. Port those plots before considering its removal.
- **`TempConfig`:** `utils/TempConfig.py` is still used by the histogram-model path. Replace it with a validated model/profile as part of config cleanup; it is not currently dead code.
- **Tests (COMPLETED FOR CONFIG):** Rewrote `test_hydra_adapter.py` and promoted `test_config_composition.py` to test the canonical 4-tier configuration system. Removed `test_hydra_adapter_new.py` and `test_config_composition_new.py`.

## Keep as reusable functionality

The storage package’s core files—`schema.py`, `hasher.py`, `paths.py`, `metadata.py`, `backend.py`, `inputs.py`, `fits.py`, and `study.py`—are the replacement foundation, not cleanup targets. Also retain the scientific/model and I/O functionality: `fit/`, `models/m2.py`, the staging/ROOT/JAGS readers and transforms, systematic algorithms, and active plotting and utility modules. In particular, `fit/create_m2uv.py` and `fit/create_evtrees.py` still support existing input or systematic workflows; removing their callers does not by itself make those algorithms obsolete.

**Storage identity is a removal gate:** current input generation uses a version string such as `cfg.staging.version` (with a `"v1"` fallback), and existence checks can reuse files without validating a manifest. Complete the staging identity/publication and artifact-reuse design before treating the old pipeline as safely replaced.

Finally, `qpym2/io/fit_data/sidecar.py` is **untracked working-tree content**. I would not classify it for deletion; review and align it with the planned artifact-manifest design before deciding whether its responsibilities belong in `storage/metadata.py` or elsewhere. No files were changed for this inventory.
