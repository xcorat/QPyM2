# QPyM2 data storage — proposed layout

**Date:** 2026-10-29  
**Scope:** Proposed storage layout only. This document does not implement the layout.

## Layout

```text
staging/<component_name.root>

fit_data/
  h2uv/
    input/
      1ton_data.h5
      1ton_data_blinded.h5
    roi_2400_2570/thr200/b2.5/smoothed/
      input/
        <study_UUID_component_name>.h5
      output/
        <study_UUID>/...
      studies/
        <study-name>.sqlite
        <study-name>.manifest
```

`staging/` holds the staged ROOT files. The `fit_data/h2uv/` subtree contains all fit data, templates, outputs, and study records for the `h2uv` coordinate system. Other coordinate systems get their own subtree and do not share inputs with `h2uv`.

The coordinate system, ROI, threshold, binning, and smoothing directories represent coordinates that define distinct template-histogram families. Inputs from different coordinate systems are never shared. Other coordinates are not added to this directory path; the study UUID distinguishes them, including cases where configurations share some template data. For example, a component file may be named `Muons_e4553df.h5`. The UUID is constructed from the relevant coordinates not represented by the directory structure.

The `output/<study_UUID>/` directory contains outputs for that study point. The study SQLite database records all relevant point IDs and their coordinates, and may include references to fit, toy, systematic, or plot results. The accompanying manifest provides a readable description of the study.

## Coordinate decisions

- **Coordinate system:** include it as a directory level above all fit data, such as `fit_data/h2uv/`. Never share input files across coordinate systems.
- **`skip_roi` and blinding:** these are data-side settings, not coordinates (for now). Keep them with the selected data's metadata; do not use them to distinguish template inputs.
- **Bias injection count:** injected signal count `n` may distinguish bias-study points and should be recorded as a study coordinate. The bias-fit `xmin` is a method setting, not a coordinate.

## Metadata

Each HDF5 file should include all resolved parameter values, including those not represented in its path, along with relevant input provenance and product-specific settings. The study database and manifest organize points and results; they do not replace the HDF5 metadata.

## Implementation note

The storage strategy is intended to be implemented in a new `qpym2/storage/` module, with new run files to create and consume these products. These are planned additions, not current code.
