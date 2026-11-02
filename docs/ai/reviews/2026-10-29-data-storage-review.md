# QPyM2 data storage — proposed layout

**Date:** 2026-10-29  
**Scope:** Proposed storage layout only. This document does not implement the layout.

## Layout

```text
staging/<component_name.root>

fit_data/
  roi_2400_2570/thr200/h2uv_b2.5/
    input/
      data.h5
      data_blinded.h5
      <component_name>_<input_UUID>.h5
    output/
      fit_<output_UUID>.h5
      fit_<output_UUID>.nc             # optional: traces retained for selected fits only
  studies/
    <study-name>/
      <study-name>.manifest
      <study-name>.sqlite
      results.csv                      # aggregated summary statistics across study points
      plots/                           # comparative study plots, pulls, systematic scans
```

`staging/` holds the staged ROOT files. The parameter subtree contains input products and fit outputs for the specified ROI, threshold, coordinate system, and binning. The coordinate system and binning are combined in one directory name (for example, `h2uv_b2.5`); other coordinate-system/binning combinations get their own subtrees and do not share coordinate-specific inputs.

Binned observed-data histograms (`data.h5` and `data_blinded.h5`) are stored alongside the component templates in the parameter-level `input/` directory, matching the exact ROI and binning geometry.

The path records the selected ROI, threshold, coordinate system, and binning. Other input parameters distinguish input products via `input_UUID`, which is a deterministic hash derived from the applied parameter set (cuts, energy shifts, staging contract version) independently of the component name. Sources evaluated under the same parameters share the exact same `input_UUID` tag. Thus multiple iterations of the same component can coexist (for example `Co60_<input_UUID>.h5` across shifts), while invariant sources visibly share the common baseline tag. Metadata in each HDF5 remains the complete record of resolved settings and provenance; the UUID in the filename provides a reproducible parameter address.

Outputs are flat files directly under `output/`, named `fit_<output_UUID>.h5`, where `<output_UUID>` is a deterministic hash of the point and method settings affecting that fit. Posterior traces (`fit_<output_UUID>.nc`) are written only for selected fits (e.g., nominal baseline or diagnostic points) rather than all fits, avoiding inode and metadata overhead on shared filesystems.

Study records, aggregations, and comparative artifacts live under `fit_data/studies/<study-name>/`. A study references input templates and output fits, but output products are not owned or duplicated by a study. The `<study-name>.manifest` describes the study definition and parameter sweep; `<study-name>.sqlite` indexes the mapping between coordinates and input/output artifacts; `results.csv` collects summary statistics; and `plots/` houses comparative figures (pull distributions, parameter scans, systematic deltas).

## Coordinate and product decisions

- **Directory coordinates:** include ROI, threshold, coordinate system, and binning. Combine coordinate system and binning in the directory name, such as `h2uv_b2.5`. Other parameters distinguish products through deterministic UUIDs, not extra directory levels.
- **Binned data in parameter input:** keep the binned observed-data histograms in the parameter-level `input/` folder alongside templates (`data.h5`, `data_blinded.h5`), avoiding ambiguous global input locations.
- **Deterministic UUIDs:** `<input_UUID>` and `<output_UUID>` are deterministic hashes generated from canonical parameter strings or matrix keys, guaranteeing cache hits on identical configurations. `<input_UUID>` identifies the parameter set (cuts and coordinate transforms) independently of the component name, allowing invariant background sources to share a common UUID tag.
- **Smoothing:** persist the original, unsmoothed template histograms. Apply smoothing during fit construction, not while creating the reusable input/template products; smoothing is a fit-method setting and is not a directory coordinate.
- **`skip_roi` and blinding:** these are data-side settings for now. Keep them with the selected data's metadata; do not use them to distinguish template inputs.
- **Bias injection count:** injected signal count `n` may distinguish bias-study points and should be recorded as a study coordinate. Bias-fit `xmin` is a method setting, not a coordinate.
- **Staging:** remains flat (`staging/<component_name.root>`).

## Study records and concurrent results

The study SQLite database is a collected index, not a shared write target for fit workers. Workers write their own outputs independently; after the study tasks finish, a collection/aggregation step creates or updates the SQLite database, publishes the study manifest, merges summary statistics into `results.csv`, and generates figures under `plots/`. Analysis/query consumers can then treat the completed study records as read-only.

For toy/systematic summaries where full traces are not retained, each task writes a small result file to a unique path (for example, one CSV per task or point). The collector validates and merges these files into `results.csv`, then records the resulting paths and summary metrics in the study database. Do not have workers append concurrently to one shared CSV: concurrent line writes can interleave or lose data. CSV is appropriate for tabular results and quick inspection; SQLite is useful as the collected index for point coordinates, statuses, and links across products.

## Metadata

Each HDF5 file should include all resolved parameter values, including those not represented in its path or UUID, along with relevant input provenance and product-specific settings. The study database and manifest organize study membership and references; they do not own output products or replace the HDF5 metadata.

## Implementation note

The storage backend is implemented in `qpym2/storage/`:
- `paths.py`: Directory hierarchy (`roi_{min}_{max}/thr{val}/{coord}_b{bin_size}/input|output`) and study routing.
- `schema.py`: Canonical schema normalization for template inputs and fit outputs.
- `hasher.py`: Deterministic 12-character SHA-256 UUID generation via canonical JSON serialization.
- `metadata.py`: Full metadata and provenance serialization into HDF5 file-level attributes.
- `backend.py`: Core I/O routines for composite templates, binned observed data, and fit results (`.h5` and optional `.nc`).
- `study.py`: Post-run SQLite indexing, manifest publication, and task CSV collection into `results.csv`.

