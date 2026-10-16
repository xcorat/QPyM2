# Decouple staging storage and add Parquet

## Goal

Add Parquet as a staging-data format and isolate storage-format details from the
rest of the analysis. Analysis code should consume staged event columns through
a format-independent interface rather than depending directly on ROOT trees,
RDataFrames, or ROOT-specific filter expressions.

## Work

- Define and document the staging schema, including column names, data types,
  and any required metadata. The current staged columns are `u0`, `v0`, `ch1`,
  `ch2`, `ds`, and `t`.
- Add Parquet staging readers and writers, and make the staging format and
  filenames configurable. Retain support for existing ROOT staging files
  during migration.
- Establish a format-independent staging I/O boundary. Keep ROOT and Parquet
  serialization and format-specific access in their own backend modules.
- Refactor staging consumers to use that boundary, including event-list and
  histogram construction, fit-input creation, and systematic workflows. Move
  filtering and event transformations into format-independent Python/NumPy
  code where practical instead of passing ROOT C++ expressions downstream.
- Replace ROOT-based histogram smoothing (`TH2D.Smooth` in
  `qpym2/utils/hists.py`) with a NumPy/SciPy implementation. Check its output
  against current ROOT smoothing behavior before adopting it.
- Inventory other ROOT-specific operations on the staging/analysis path and
  either move them into the ROOT backend or replace them with Python/NumPy
  implementations. Keep unrelated ROOT use, such as raw ROOT input and
  ROOT-based plotting, out of scope unless separately agreed.

## Validation

- Compare ROOT and Parquet staging schemas, event counts, and selected events
  on representative data and MC samples.
- Compare histograms and relevant systematic outputs produced from each
  staging format, including smoothing behavior.
- Verify that downstream fit-input creation works from both formats during the
  compatibility period and that read/write failures are reported rather than
  silently treated as empty data.
- Benchmark file size and representative read/filter time to confirm the new
  format is useful for this workload.

## Scope note

Using Parquet for staging will not by itself remove ROOT from the full
pipeline: raw inputs and some separate systematic or plotting paths still use
ROOT. The immediate goal is to decouple staging storage and downstream
analysis from ROOT where feasible.
