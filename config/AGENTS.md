# config/ — Hydra configuration

See the root `AGENTS.md` for how config is loaded (`qpym2.cfg.get_config`, not
CLI flags).

## Live layout (source of truth)

`config/config.yaml` is the master and composes four tiers via Hydra `defaults:`:

```
defaults:
  - environments: base     # @package _global_ — paths, system, qroot
  - analysis: models_base  # @package _global_ — models/ROI/cuts/binning
  - hyperparams: base      # co_cuts, ndbd_cutoff, co_shift, blinded, ...
  - toys: base             # toy-run grid (loaded separately)
  - _self_
```

- `environments/*.yaml` — infrastructure only: `infrastructure.directories`,
  `output_paths`, `data_files`, `naming_conventions`, `qroot`, `system`.
  No physics. Common: `base`, `tech_note`, `nersc`, `m2v5-1`, `m2v5-1-250`.
- `analysis/*.yaml` — physics only: `models.*` (ROI `esum_min/esum_max`,
  `threshold`, `binning.bin_size`), `cuts.*` (base/mc/data, `co_cuts`,
  `ndbd_cutoff`), `components`. Common: `models_base`, `tech_note`, `m2`.
- `hyperparams/*.yaml` — systematic selection keys: `co_cuts`, `ndbd_cutoff`,
  `co_shift`, `ndbd_shift`, `blinded`, `skip_roi`, `smooth`.
- `toys/*.yaml` — loaded via `qpym2.cfg.load_toy_config` (see `toys/README.md`).

Files use the `.yaml` extension and are selected by name, e.g.
`environments=m2v5-1 analysis=tech_note`.

## Docs

`README.md` and `config_guide.md` describe the live Hydra layout. If they ever
disagree with the YAML, the YAML wins. (A stale aspirational `extends:` /
`tools/*.yml` guide was removed.)

## Editing rules

- Keep paths (environments) and physics (analysis) separated.
- Rotated-axis (`u/v`) ROI limits are derived in code from `esum_min/esum_max`
  and `threshold`; do not add `u/v` limits to config.
- Only set `models.*.binning.bin_size` (defined on the `e1/e2` axes).
- `qm2style.mplstyle` is the shared matplotlib style.
