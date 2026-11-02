# TODO: Shared fit-data objects and input/result sets

- **Status:** In Progress (Core Storage & I/O Modules Complete)
- **Scope:** Derived input data, fit setups, fit results, and toy-fit results. Staged files are treated as fixed upstream inputs.
- **Specification:** [`docs/data-layer-spec.md`](../data-layer-spec.md)
- **Related:** [`model-parameter-taxonomy.md`](model-parameter-taxonomy.md), [`data-management-and-artifact-store.md`](data-management-and-artifact-store.md)

## Decisions and terminology

- Call event samples from the experiment **observed data** (short object kind: `observation`). “Measurements” is less suitable for event-level arrays and can imply already-aggregated measured quantities.
- Call histogram payloads consumed by the likelihood **templates**. Record whether a template is an MC background template or an NDBD signal template. “Simulation” describes the MC source/sample, not the derived histogram itself.
- Use **fit** for inference on observed data. Use **toy experiment** or **pseudo-experiment** for generated pseudo-data, and **toy fit** for a fit to that pseudo-data. A toy experiment and its fit result are separate products.
- Keep one store root with folders by object kind; keep observations/templates, input sets, and results independently addressable. Model/campaign names are associations or views, not payload owners.
- For HDF5 objects, embed the canonical manifest in the file. An input-set HDF5 contains references/metadata rather than copies of member arrays. A fit bundle embeds its manifest in `fit.h5`; `fit.nc` carries the matching IDs. Sidecars are optional conveniences or fallbacks for formats that cannot carry metadata.
- Add new fit-data I/O objects alongside the current `qpym2.io` modules. Leave existing I/O implementations intact; migrate callers incrementally through the new API and a compatibility adapter.
- Object IDs must round-trip between an ID and its normalized JSON/key-value identity. IDs need not be hashes. A versioned identity profile supplies defaults for omitted ID keys; defaults are resolved before encoding and decoding.
- Only keys that identify the relevant object belong in its ID. Additional resolved information—including component lists, actual bin-edge arrays, and provenance—can be stored in metadata. If metadata describes a value that changes a scientific payload, the identity must still include a key or profile reference that distinguishes that change.

## Proposed product boundaries

1. **Observation objects:** selected event-coordinate arrays, separately identified for observed/blinded variants and for selections that produce different arrays.
2. **Template objects:** independently addressable source-level background histograms and NDBD signal histograms. Preserve source-level granularity so component regrouping does not force template rebuilding.
3. **Input sets:** manifests that compose an observation, template references, and the histogram geometry/selection metadata needed by a particular input setup. Sets with identical members and effective settings share an ID.
4. **Fit setups:** an input-set reference plus fit-level choices such as component grouping, smoothing, likelihood, and priors. Fit-only changes reuse input objects.
5. **Results:** fit-specific outputs in a results area, linked to the fit setup/input set. Toy experiment data and toy-fit results remain distinct and link to their generator/reference fit as applicable.

Initially keep small metadata (component bookkeeping, edge arrays, source summaries) in input-set or fit-result metadata. Avoid one physical file per small array: object identities may resolve to datasets within immutable HDF5 payloads/shards.

## Implementation plan

1. [x] Inventory input/output datasets and their current readers. Confirm which consumers require the selected event array in `fit.h5`; avoid duplicating it in results if an input-set reference is sufficient.
2. [x] Define versioned identity profiles by object kind (`qpym2.storage.identity_profiles`). Specify the ID keys, canonical value normalization, profile defaults, and reversible JSON/key-value encoding. Preserve old profile decoders when adding keys or changing defaults; never reinterpret an old ID using new defaults.
3. [x] Add a single-root path resolver (`qpym2.storage.fit_data_store.FitDataStore`) and HDF5/NetCDF embedded-metadata adapter (`qpym2.io.fit_data.embedded_metadata`); keep the current sidecar store available for v0.6 objects.
4. [x] Add new observation/template/prior-input object readers and writers (`qpym2.io.fit_data.objects`). Start with observations and source-level templates; keep small metadata with the set unless measured reuse justifies separate objects.
5. [x] Add input-set composition and an optional HDF5 materializer for the unchanged legacy reader (`qpym2.io.fit_data.input_set`). A materialized file is shared by composition ID, not copied per model label.
6. [x] Add a fit-result bundle reader/writer for `fit.h5`/`fit.nc` (`qpym2.io.fit_data.fit_result`). Record the exact input-set and fit-setup IDs consumed.
7. [x] Import existing input bundles with the new legacy adapter (`qpym2.io.fit_data.legacy.decompose_legacy_input`, `verify_legacy_roundtrip`).
8. [ ] Route fit runners and result readers through the new objects.
9. [ ] Validate sharing with model/fit variations before migrating toy-experiment and toy-fit I/O.

## Acceptance criteria

- Encoding an ID and decoding it returns the same normalized identity; sparse JSON input resolves through documented, versioned defaults.
- Different model labels with identical effective input products refer to the same input set.
- A variation creates new identities only for the affected observation/template objects and the resulting input set.
- Component grouping or smoothing changes do not duplicate input objects but do distinguish fit setups/results when scientifically relevant.
- A fit result records the exact input set and fit setup used; readers can resolve them without reconstructing a model directory path.
- The HDF5 compatibility materializer round-trips through existing fit-input readers.
- Large reused payloads are not copied for each input set or model view; physical packing is measured to avoid excessive small-file counts.
- Existing `qpym2.io` modules remain unchanged during initial implementation; new object APIs and compatibility tests cover their current read/write contracts.
