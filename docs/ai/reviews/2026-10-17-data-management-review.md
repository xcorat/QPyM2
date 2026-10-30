# QPyM2 data management and artifact storage review

- **Review date:** 2026-10-17
- **Scope:** Persistent storage and discovery of staged data, fit inputs, fits, toys, and systematic outputs; artifact reuse, metadata/catalog options, and Python access patterns.
- **Context:** QPyM2 runs on shared NERSC filesystems and submits independent tasks through SLURM arrays as well as local pools. The goal is a small, adaptable storage layer that keeps the underlying scientific files directly accessible.

## Requirements and decision criteria

These requirements are the basis for the recommendations below. They separate what the native conversion must provide from integrations that can be added later.

| Priority | Requirement | Design consequence |
|---|---|---|
| **Core** | Multiple models/campaigns may use the same staged data, arrays, histograms, templates, fit inputs, or other products. Avoid copying identical payloads just to preserve a model folder layout. | Model names are aliases/associations, never storage keys. Reuse derivations from effective inputs/settings and store each exact payload once by content hash. |
| **Core** | Reuse should be granular: share an event array or histogram even if a model's full fit-input composition differs. | Represent stable, expensive/reused data products independently; let fit-input manifests compose references. Keep small/rare intermediates packed where separate files would overload filesystem metadata. |
| **Core** | A missing derived file should be straightforward to regenerate. The project should know what produced it and what it needs. | Store operation settings, exact input identities, schema, implementation/environment identity, and seed policy; provide an explicit `ensure/materialize` path that reuses valid ancestors and reports unavailable dependencies. |
| **Core** | Preserve scientific provenance and prevent stale/incompatible cache hits. | Separate a pre-compute derivation key from post-write content IDs; validate schema/checksum and publish outputs immutably. Never overwrite a completed result silently. |
| **Core** | People should be able to identify and browse data without opening a Python session, while Python users retain direct access to raw files. | Keep readable model/campaign views, manifests, and TSV/CSV indexes; resolve to ordinary paths so ROOT, h5py, NumPy, xarray, and pandas can read the payloads directly. |
| **Core** | Fit tasks use one node apiece, with multiple independent fits running concurrently across nodes as SLURM capacity and project limits permit. | Preserve the current TSV/SLURM array or packed-job execution during conversion; the campaign plan must map task indices to stable artifact IDs and support a concurrency cap. |
| **Core** | The same analysis jobs must also run locally without SLURM or a cluster, even if serial or much slower. | Keep task functions and artifact resolution executor-independent. Support a direct serial/local-pool backend as well as the SLURM array adapter; both consume the same resolved task records and publish the same artifact format/identity. |
| **Core** | Keep the architecture minimal and layered; existing scientific file formats and lower-level I/O remain usable. | Filesystem is the payload layer, manifests carry durable provenance, and SQLite is a rebuildable search index—not a binary database or multi-node write service. Do not require a service for local/ordinary analysis. |
| **Later** | Align with common ML and scientific workflow tools without coupling the data model to one engine. | First convert the existing repository onto native artifact/task contracts. Add Snakemake for DAG/execution and MLflow for run tracking later as optional adapters over those contracts. |

## Executive summary

QPyM2 currently locates related artifacts through model-named directories and fixed filenames. That makes a path serve two jobs: it is both a storage address and an implicit declaration of which model/parameter choice produced the data. This is fragile: model identity is not a stable property of a directory name, the naming scheme can omit relevant settings, and shared upstream files can be copied into multiple model folders.

The recommended direction is a **single logical artifact store root** on the existing filesystem. A model label is a many-to-many association, not a storage key: when different model points need the same effective stage/input, they should resolve to the same artifact record rather than copy its data. Separate the identity of a *derivation* from the identity of its actual bytes: a derivation key enables reuse before rerunning work, while a content hash stores each exact payload once even if it came from different derivations. Keep the original format files (ROOT, HDF5, NetCDF, CSV, and potentially Parquet) as ordinary files. A compact metadata catalog maps model/campaign relationships to shared artifact records and blob paths; it should not contain the large binary data or become an opaque replacement for them.

The best minimal fit is a **filesystem artifact store + immutable JSON manifests + SQLite search/index catalog**. On one machine or a single-writer process, SQLite is a good embedded database. For SLURM jobs that may write concurrently from multiple compute nodes, do not have every task update one shared SQLite file directly. Let workers publish their own complete artifact directories and manifests, and have a controlled indexer update/rebuild the SQLite catalog. If concurrent, central metadata writes later become a real requirement, move only the catalog backend to a client/server database such as PostgreSQL.

This borrows the separation used by MLflow (metadata store vs artifact store) and content-addressed reuse used by DVC, while avoiding the need to deploy either as QPyM2 infrastructure. It complements the campaign/artifact identity design in the 2026-10-15 workflow review: that document defines *what an artifact means and depends on*; this one defines *where it lives, how to discover it, and how to avoid storing/recomputing it twice*.

## Current QPyM2 behavior and implications

- Environment config currently defines separate locations such as `staging_dir`, `fit_dir`, `systematic_dir`, and `data_dir`; e.g. `config/environments/tech_note.yaml` points stage and fit products at distinct trees.
- Fit runners resolve a model output directory and expect familiar names such as `input.h5`, `fit.h5`, and `fit.nc`. Toy and systematic runners reconstruct locations from `fit_dir`, the model configuration, and naming helpers.
- `qpym2/cfg/model_utils.py` derives a folder name from selected cut/cutoff/smoothing settings. `get_model_dir_from_config()` uses the base name `model`; model configuration identity is not itself the directory identity. Some run paths also store shifts in filename suffixes.
- This convention groups input, fit, toy, and other files under a selected model point, even when some objects are independent of that point or are reusable by another campaign/model point. Staging files are especially likely to be reused. Repeating a large stage or input file merely to maintain this folder layout wastes space and makes invalidation hard to reason about.
- A folder name can be renamed, truncated, manually edited, or collide when model definitions evolve. It cannot reliably distinguish “same human label, different resolved config” from “different label, same resolved config.” The model key and resolved config belong in metadata; the path should be a storage implementation detail.

## Established patterns worth borrowing

| Pattern | What it provides | QPyM2 fit |
|---|---|---|
| **MLflow tracking/artifact separation** | A metadata backend records run IDs, parameters, metrics, and tags, while a separate artifact store holds large output files. The metadata backend can be SQLite for local use or a server database for team use. | Strong conceptual precedent for keeping catalog records small and scientific files in regular storage. MLflow's tracking server and broader platform are unnecessary unless a concrete service/UI need appears. |
| **DVC content-addressed cache and run cache** | Data is held once in a cache by content hash, while lightweight references expose it in workspaces. Repeated data contents can share a cache object. Run-cache keys reflect dependencies/parameters/commands, making it possible to skip equivalent deterministic stages. | Useful reuse and indirection principles. QPyM2 needs domain-specific artifact identity and stochastic-fit handling; adopting the full DVC pipeline layer is optional and not required for this design. |
| **Scientific/analytics file formats** | ROOT, HDF5, NetCDF, and CSV remain self-contained files with mature domain readers. Parquet is a portable column-oriented format for event tables and analytical scans; Zarr is a chunked, compressed N-dimensional array store. | The catalog should identify and describe files, not force a format replacement. Parquet is a plausible staging backend already considered in `docs/todos/staging-parquet.md`; Zarr should only be considered if access patterns demand chunked array access. |
| **Data-lake table catalogs (Iceberg/Delta)** | Transactional snapshots, schema evolution, table-level metadata, and query planning for large, regularly queried collections of tabular files. | Too much machinery for a collection of stage products and fit bundles today. Revisit only if staged events become a shared, queried data lake with concurrent table writers and evolving schemas. |

### Cross-industry lessons

1. **Separate metadata from the large payload.** The catalog stores IDs, attributes, lineage, and paths; HDF5/NetCDF/ROOT/etc. hold numerical arrays and event records.
2. **Separate logical identity from human-readable names and physical paths.** Names are labels for people; a stable ID is what enables reliable lookup and reuse.
3. **Use a content-addressed or identity-addressed cache.** Reuse an existing artifact when its actual inputs and relevant settings match, not just because a conventional filename exists.
4. **Make files directly usable.** Storage/catalog helpers should return a `pathlib.Path` (or a URI later), so current readers such as PyROOT, h5py/PyTables, xarray, and pandas can be used without a proprietary wrapper.
5. **Keep a portable manifest with the object.** A database index is useful, but should not be the only way to learn what an artifact contains or recover its identity.

## Minimal target architecture

### 1. One logical store root, grouped by artifact type

Use one configured store root (for example, `artifact_store`) instead of requiring consumers to know separate model-folder conventions for every output class. Subdirectories may separate object types and catalog files; “one store” does not mean putting every file into one physical directory.

Illustrative layout:

```text
<artifact_store>/
  catalog.sqlite                 # searchable index; rebuildable from manifests
  artifacts/
    staging/ab/<artifact-id>/    # manifest + links to staged payload files
    event-array/ab/<artifact-id>/
    histogram/ab/<artifact-id>/
    fit-input/ab/<artifact-id>/  # manifest + optional shared input.h5 bundle
    fit/ab/<artifact-id>/        # manifest + fit.h5 / fit.nc links
    toy/ab/<artifact-id>/        # one batch/shard, as appropriate
    systematic/ab/<artifact-id>/
  blobs/sha256/ab/<digest>       # immutable payload bytes; one copy per exact content
  campaigns/<campaign-id>/
    campaign.yaml               # submitted intent
    manifest.json               # resolved plan and references to artifact IDs
```

The artifact-record path is generated from `artifact_type` and opaque `artifact_id`; it does not encode `model_name`, cut labels, or the current study name. Payload paths are generated from their `content_id`. Sharding by ID prefixes avoids accumulating a huge number of entries in one directory. A campaign manifest stores references to shared artifact records rather than copies of their data. Raw input files remain external, read-only inputs referenced by canonical source paths and whatever stable size/checksum/version data is available.

The exact root may be on NERSC Community Filesystem or another configured filesystem. Environment config chooses the root; artifact IDs and manifests do not depend on its absolute location, so copying the store or changing the mount does not redefine the science.

### 2. Model and artifact identities are explicit metadata

Use separate concepts rather than trying to make one folder name carry them all. Most importantly, a model name is an **association/label**, not an artifact key. The same artifact can be referenced by multiple model points or campaigns when the data-producing dependencies are equivalent.

- `model_key` / `model_name`: human-facing label, e.g. `tech_note_200`; may be many-to-many with data artifacts.
- `model_config_id`: digest of the canonical resolved model configuration. Model parameters should have a **static, declared schema** for a given software or major campaign version rather than arbitrary freeform keys. Static schemas ensure deterministic key ordering, simplify canonical serialization and hashing, catch configuration typos early, and allow physics parameters to map directly to indexed SQLite columns.
- `model_point_id`: identity of a resolved model plus the static input-affecting point settings (e.g. `co_cut`, `ndbd_cutoff`, `co_shift`, `ndbd_shift`).
- `artifact_id`: **scientific derivation/cache key** for one logical output bundle. Hash the artifact type/schema version, the named input payload content IDs (or canonical external-source identities), and the canonicalized scientific settings that directly affect the numerical computation. Upstream artifact IDs are kept in the manifest for provenance. Do not include arbitrary runtime labels or aliases in this key. This key is computable before execution and enables skip/reuse.
- `execution_id` / `run_metadata`: environment and execution tracking stored in the manifest and catalog. Software versioning, git commit, working tree dirty state, execution hostname, and Slurm task IDs belong here—generating execution run hashes (similar to random seeds or run attempt IDs) rather than invalidating the scientific derivation key. This avoids spurious cache misses during iterative analysis development when code edits do not alter numerical logic.
- `content_id`: content address for one payload file, e.g. `sha256:<digest-of-bytes>`, calculated after writing. Different derivation records that produce byte-identical files may reference the same immutable payload. The content hash also verifies integrity, but by itself cannot skip the first execution.
- `campaign_id`: identity for a requested study. Campaigns refer to artifacts; a new campaign can reuse objects produced for an earlier one.

The same model label with a changed config yields a different config ID. Two model labels resolving to the same canonical config can point to the same config identity while preserving both labels as aliases. More importantly, two different model configs can still share a stage or input artifact if the *effective settings relevant to that stage/input* are identical. An irrelevant fit-only option must not change a staging or input ID; a changed selection that affects stored rows must. A fit ID must differ when the model likelihood or other fit behavior differs, even if the two fits use the same HDF5 input.

These are two levels of addressing, not competing choices: use the derivation key for pre-compute reuse; use content-addressed blobs for exact-byte deduplication. A derivation key is a hash of the computation recipe and its effective input content, while a blob content ID is a hash of output bytes. Exact-byte deduplication will not combine numerically equivalent but byte-different outputs (for example, if files embed timestamps or stochastic samples differ), which is why dependency-aware derivation identity remains essential. Any side input or model detail that affects the numerical result must appear in the derivation descriptor even if it is not stored inside the main input file.

Each derivation manifest should minimally record the operation and its dependencies, but **not one model label as if it owned the payload**:

```json
{
  "artifact_id": "<derivation-id>",
  "artifact_type": "fit-input",
  "schema_version": 1,
  "upstream": {"staging": "<artifact-id>"},
  "operation_config": {"static_parameter_settings": "..."},
  "files": [{"path": "input.h5", "size": 123, "content_id": "sha256:<digest>"}],
  "execution_metadata": {
    "execution_id": "<run-uuid>",
    "code_revision": "<git-commit-or-digest>",
    "git_dirty": false,
    "seed": 12345
  },
  "created_at": "<UTC timestamp>"
}
```

The canonical resolved operation settings and upstream IDs are retained for auditability. Model labels/configs/points are stored as separate records and associations; they can refer to one artifact without being copied into its identity or payload. The manifest format should be versioned and permit future fields. It can later map to richer provenance standards without requiring a database migration first.

### 3. SQLite as catalog/index—not as a container for scientific data

SQLite is a reasonable first catalog because it is embedded, standard in Python, easy to back up, queryable with SQL, and requires no service. Store references and searchable metadata, not binary fit/staging payloads. A small schema could have:

- `artifacts(artifact_id, artifact_type, schema_version, status, created_at, metadata_json)`
- `artifact_files(artifact_id, relative_name, content_id, size_bytes)`
- `blobs(content_id, path, size_bytes, checksum_algorithm)` for immutable content-addressed payload files
- `artifact_edges(parent_id, child_id, role)` for lineage and reuse queries
- `models(model_config_id, canonical_config_json)` and `model_labels(model_key, model_config_id)`
- `model_points(model_point_id, model_config_id, point_json)` and `model_point_artifacts(model_point_id, artifact_id, role)` for many-to-many associations
- `campaigns(campaign_id, manifest_path, created_at, status)`
- `campaign_artifacts(campaign_id, artifact_id, role)` for campaign-to-object references

Index common query columns such as type, status, and content ID. Join through the model-point/campaign association tables to answer “which models use this artifact?” Preserve complete, typed operation configuration and provenance in JSON manifests rather than building a generic EAV database for every physics parameter. A small Python query API can provide filters for the frequent cases; SQL remains available for exploratory use.

#### SQLite and SLURM concurrency

The SQLite choice depends on where and how it is accessed. SQLite allows only one writer at a time. Its WAL mode specifically does not work across a network filesystem because WAL requires shared memory among participating processes. SQLite also cautions against direct simultaneous access from many computers over network filesystems where locking semantics may not be reliable. NERSC Global File Systems are shared across NERSC systems; actual behavior of the chosen mount must be checked rather than assumed.

For QPyM2, use this write pattern initially:

1. The planner assigns each task a derivation `artifact_id` and record path; the key is independent of model/campaign labels and includes only relevant dependencies/config.
2. A task writes under a unique temporary directory on the same filesystem.
3. It validates outputs, computes each file's content ID, and atomically installs a missing immutable blob with create-if-absent semantics; if another task already installed that hash, it verifies and reuses the existing blob. The artifact manifest and file-to-blob links are published only after validation.
4. A single process (for example, a campaign collector or explicit `catalog index` command) imports manifests into `catalog.sqlite`. Workers do not concurrently mutate the shared catalog.
5. The catalog can be rebuilt from manifests and directory contents if it is lost or stale.

For local interactive use, a user may update SQLite directly with short transactions. If multiple nodes need live concurrent catalog updates, use a supported client/server metadata service (e.g. PostgreSQL) rather than turn SQLite into an unverified distributed database. Test temporary-file publication/rename semantics and catalog access on the actual target filesystem.

### 4. A thin Python API over direct files

The storage API should hide path construction and validation, not the science formats. A narrow interface could look like:

```python
artifact = store.get(artifact_id)
path = artifact.file("input.h5")
with h5py.File(path, "r") as f:
    ...

inputs = store.query(kind="fit-input", model="tech_note_200", point={"co_shift": 0.5})
fits = store.query(kind="fit", model="tech_note_200", point={"co_shift": 0.5})
```

Useful first operations are `register/publish`, `get(id)`, `query(filters)`, `verify(id)`, `path(id, filename)`, and an explicit `ensure/materialize(id)`. A model-aware query joins through model-point/campaign associations to return artifacts; the returned artifact itself is not owned by that model label. `get()` returns a typed record containing its manifest and ordinary local paths. Existing HDF/ROOT/NetCDF readers continue to open those files directly. If a user needs lower-level access, they can read `artifact.path` and bypass catalog convenience calls without copying or unpacking data.

Avoid a universal storage plugin framework initially. A configured root and `pathlib` implementation is enough; a later filesystem backend can accept `Path` or URI through a small protocol if data moves to object storage. Keep data format adapters (`io/rdf.py`, `io/hdf.py`, future Parquet reader) below the artifact/catalog API.

### 5. Reuse and duplication policy

- **Use the derivation key before execution.** If an artifact with the same valid artifact ID already exists, reuse it after checking its complete manifest, schema, and linked blobs. The key is based on effective operation settings and input content IDs, so changing a model label or campaign does not trigger a new computation.
- **Store each exact payload once.** Hash output files when they are written and store them immutably by content ID. Different derivation IDs can share byte-identical payload files. A blob hash verifies bytes; it does not by itself say which model/configuration produced them.
- **Keep model membership many-to-many.** Model keys, config IDs, and model points link to artifact IDs through catalog records and campaign manifests. Do not put one `model_key` in a canonical artifact record as its owner or require one physical copy per model view.
- **Keep alias-only metadata out of shared payload bytes.** Store model labels and campaign names in manifests/catalog associations. If a ROOT/HDF5 file embeds those labels, its bytes will differ and byte-level deduplication will correctly treat it as a different blob; keep such fields in the sidecar when they are not part of the scientific payload.
- **Make stochastic identity explicit.** Toy and fit seeds or seed policies belong in the identity. A retry/execution attempt can be separately recorded without changing the scientific point; retries should not overwrite earlier completed bytes silently.
- **Keep bundles together.** `fit.h5` and `fit.nc` from one fit can be files in one fit artifact directory, with one manifest listing both. The directory is a bundle, not a hierarchy of independent copies for every campaign.
- **Use links to expose shared content.** A named file in an artifact record can be a symlink to its immutable content-addressed blob; a human-friendly `views/<campaign>/<model-label>/...` entry can point to that artifact record. Both levels preserve normal file access without copying payload bytes. The manifest/content ID remains authoritative, not the symlink name. Confirm filesystem link support before relying on hard links; never mutate a shared blob in place.
- **Make cleanup reference-aware.** Garbage collection must preserve artifacts referenced by retained campaigns or downstream artifacts. Retention/deletion should be explicit and report what will be removed.

### 6. Make derived products lazily regenerable

Treat a derived payload as a reusable cache entry backed by a durable recipe/lineage record. A path in a query result is a current location, not the identity. If its blob is gone, the artifact record still provides enough information to plan its regeneration:

- artifact type/schema and expected named outputs;
- canonical operation settings and input content IDs (plus logical upstream artifact IDs for provenance);
- raw source file identity/locations and the model-independent selection/transform settings used;
- code revision or source digest, dependency/container/environment identity, and random seed policy/seed;
- validators needed before publishing the regenerated output.

An explicit API/CLI such as `store.ensure(artifact_id)` / `qpym2 artifact ensure <artifact-id>` should follow this behavior:

1. Resolve the artifact manifest and verify each required blob. If present and valid, return ordinary paths immediately.
2. If a blob is missing, resolve the operation that produces it and recursively check its input artifacts. Reuse verified ancestors; schedule only missing or invalid nodes.
3. Run the resulting plan with the selected backend. Local access may materialize synchronously; a SLURM campaign should emit/submit the missing task records and return the task status rather than unexpectedly launching a long analysis inside `get()`.
4. Write to a temporary location, validate schema and outputs, compute content IDs, publish immutable blobs and the manifest, then refresh/rebuild the catalog index.
5. If required source files, code/environment, or non-file metadata cannot be resolved, report the exact blocked dependency rather than return a stale path or silently substitute another version.

This gives a simple “find, verify, and materialize if needed” experience while retaining lower-level control. Keep `get/path` as fast lookup; make `ensure/materialize` the explicit regenerating operation. Durable manifests and recipes are retained even when policy allows some derived blobs to be evicted. Raw data and other non-regenerable authoritative products need their own retention policy. A deleted SQLite index is also recoverable from manifests; it should not force recomputation of valid blobs.

To make recovery scientifically sound, the derivation key must include every behavior-affecting setting and the manifest must preserve the exact dependency/environment/seed identity. For stochastic jobs, record the seed and keep execution-attempt IDs separate; a rerun with different bytes under the same expected key is a conflict to diagnose, not permission to overwrite the old result.

## Options considered

| Option | Benefits | Costs / risks | Assessment |
|---|---|---|---|
| **Keep per-model folders, improve names** | Smallest change; compatible with current runners. | Still couples identity to mutable naming; shared stages/inputs duplicate; path parsers remain implicit APIs. | Good only as a temporary compatibility view. |
| **One folder per campaign, self-contained outputs** | Portable, easy to archive and understand as one run. | Duplicates large common stage/input files across campaigns; difficult global search/reuse unless a second catalog is added. | Use for manifests and summaries, not as the sole payload store. |
| **Shared typed artifact store + manifests + SQLite index** | Reuses computations by derivation key and stores each exact payload once by content hash; easy to query; no service; direct file access; simple migration path. | Needs stable keys, content-hash publication/validation rules, and serialized shared SQLite writes on HPC. | **Recommended.** |
| **DVC** | Mature CAS cache, references, remote storage, data versioning, pipeline caching, Git-friendly review. | Adds a tool and workflow/storage conventions; DVC run caching assumes stage determinism and does not define QPyM2 physics artifact identities by itself. HPC integration and shared-cache behavior need a pilot. | Evaluate if data versioning, remote transfer, or team collaboration needs justify it; borrow CAS ideas now. |
| **MLflow** | Mature experiment tracking and metadata/artifact-store separation, searchable run metadata, local SQLite and server-backed options. | Tracking server, client lifecycle, and its run/model abstractions exceed the immediate need; scientific stage reuse/lineage still needs QPyM2 semantics. | Reference architecture; consider if a shared UI and tracking service become explicit requirements. |
| **Store file payloads in SQLite** | One file, transactional metadata/payload bundle. | Awkward for large ROOT/HDF5/NetCDF datasets, independent format readers, partial access, and multi-node workflows. | Do not use for payloads. |
| **Delta Lake / Iceberg catalog and tables** | Strong for large tabular datasets with snapshots, transactions, schema evolution, and query optimization. | New runtime/catalog and operational model; not a natural fit for individual heterogeneous fit artifacts and present filesystem workflow. | Defer. |

## Data formats: keep the catalog decision separate

An artifact catalog does not imply a universal data file format. Preserve established files where they work:

| Data | Current/likely format | Guidance |
|---|---|---|
| Raw CUORE/MC source | ROOT and external source files | Treat as immutable external inputs; record canonical file list and whatever stable size/checksum/version data is available. Do not duplicate these in the artifact store. |
| Staged events | ROOT today; Parquet is a candidate | Keep a backend-neutral reader/writer contract. Benchmark Parquet on representative stage data and access patterns as proposed in `docs/todos/staging-parquet.md`; don't migrate merely to match the catalog design. |
| Histogram/template arrays and fit input | HDF5 today; `.npy` for simple arrays or Zarr for chunked arrays are candidates | Give reusable histogram/template objects their own derivation IDs and content IDs. A fit-input manifest can compose references; retain `input.h5` as one shared bundle per unique composition where current readers require it, not one copy per model. |
| Posterior/fit products | NetCDF and HDF5 | Bundle related files under one fit artifact ID, preserving use of the existing readers. |
| Toy/systematic products | CSV or current stage-specific formats | Keep small tabular results easy to inspect; define schema and link them to their source fit/input IDs. |
| Large chunked numerical arrays (if introduced) | Zarr or HDF5, based on access pattern | Zarr's chunking/compression and directory/object-store backends are useful where partial concurrent array access is needed. It is not a metadata database or automatic replacement for HDF5. |

Parquet is optimized for column-oriented event/table queries; it does not replace the metadata catalog. Likewise, choosing SQLite for metadata does not mean putting ROOT/HDF5 content in SQL blobs.

## A lightweight data-lake layer: reusable arrays, histograms, and fit bundles

Yes: it is useful to treat staged arrays and derived histograms as **separate reusable data products**, with fit inputs as a composition of references to those products. That is a better sharing boundary than treating the entire model directory (or each model's `input.h5`) as indivisible. The payload bytes live in the shared content-addressed blob layer; typed manifests say what each array/histogram means, and model points/campaigns point to the object IDs they use.

| Logical product | Example contents | Inputs that determine its derivation key | Likely sharing boundary |
|---|---|---|---|
| **Staged event data** | Named columns/arrays such as `u`, `v`, `ch1`, `ch2`, `ds`, and `t`, with run/sample partitions | Raw source identity, staging transform/version, common preselection, schema | Share broadly across models if the stage stores a sufficiently broad selection; model-specific cuts belong downstream unless they change what was staged. |
| **Histogram/template** | `counts`, `sumw2` (or other variance), axis edges, component/source ID, units and normalization metadata | Input data content IDs, component/sample, cuts, weights, bin edges, and transformations that actually affect the histogram | Reuse between model points when the effective histogram definition is the same. A different ROI, cut, binning, shift, or smoothing step gets a distinct key only if it changes the histogram being stored. |
| **Fit-input composition** | References to observed-data and template histograms plus the fit-input schema/config | Exact member artifact/content IDs and any input-assembly settings | Prefer a small manifest listing members. If existing fit code requires `input.h5`, materialize one HDF5 bundle per *unique composition*, then let all model/campaign views reference that bundle. |
| **Fit result (nominal/single)** | `fit.h5`, `fit.nc`, diagnostics and fit manifest | Fit-input content IDs, effective likelihood/prior/sampler settings, seed | One object per distinct actual fit computation. Full traces kept for nominal/detailed fits. |
| **Toy / systematic batch results** | Consolidated table (`results.parquet`/`.h5`), optional bundled trace file | Fit-input content IDs, batch specification, seed policy | Summary fit results (medians, CIs, diagnostics) saved in a single table per batch; traces optionally combined in a single file to prevent inode exhaustion. |

Illustrative logical structure (the file references can be symlinks to blobs):

```text
artifacts/
  event-array/<artifact-id>/manifest.json
  histogram/<artifact-id>/manifest.json
  fit-input/<artifact-id>/manifest.json  # lists histogram artifact IDs/content IDs
  fit/<artifact-id>/manifest.json
blobs/sha256/<prefix>/<digest>.npy        # a physical array payload stored once
blobs/sha256/<prefix>/<digest>.h5         # or ROOT, NetCDF, Parquet, etc.
model-points/<model-point-id>.json        # aliases/relationships, not payload copies
```

For example, the `model=tech_note_200` and `model=tech_note_250` views can both point to one fit-input manifest listing the shared observed-data histogram and component templates. If they require different signal templates or bin edges, only those products and the resulting fit-input composition get new artifact keys. When the complete fit-input composition is identical, it has one derivation key and, if materialized to HDF5, one HDF5 payload blob regardless of how many models or campaigns refer to it.

### Choose logical granularity without creating a small-file problem

Make each expensive or independently reused array/histogram an addressable logical object, but do **not** persist every transient NumPy intermediate or assume every histogram needs its own physical file. One file per small array can create excessive filesystem metadata traffic and inode use; on the other hand, one monolithic HDF5 per model prevents sharing objects between model points.

Keep logical object IDs/manifests separate from physical packing. Store a large/reused array or histogram as its own blob; pack many small, rarely reused values into a shard/container and record each member's group/path in the manifest. Benchmark representative reads and NERSC filesystem behavior before choosing shard sizes. `.npy` is a direct representation for a simple NumPy array; Zarr is useful for chunked/partial array access; HDF5 can retain named datasets/groups and can serve as a compatibility bundle; Parquet is a candidate for event tables. The catalog and content-addressing scheme need not force one format.

The model/point-to-product relation should be many-to-many. This permits an operator to browse a model-specific view while two model labels point at one physical histogram or fit-input blob. The model label and interpretation live in the association/manifest; the payload holds the numerical array or file. This is the data-lake-like property QPyM2 needs—shared, discoverable, versioned data products—without requiring a lakehouse runtime or central distributed database.

## Target query examples

The architecture should make common analysis questions straightforward without knowing folder spellings:

- Find the input artifact for model config `tech_note_200` at a selected cut/shift point.
- Find all fit artifacts produced with a particular fit profile and list their checksums/statuses.
- Find campaigns that consumed a particular staging artifact.
- Find the toys and systematics derived from a given fit artifact.
- Resolve an artifact ID to a direct path, verify its manifest/checksum, and load it in a notebook with the current HDF5/NetCDF/ROOT reader.
- Identify storage duplicates and unreferenced artifacts before cleanup.

Model labels and parameter values are query fields, not path parsers. Unknown model configs should be rejected at plan construction or explicitly recorded as legacy/unregistered; they should not silently resolve to whichever folder name happens to exist.

## Follow-up: fit multiplicity, readable paths, campaigns, and tracking

### Treat each fit as a fit artifact; use role/lineage to describe its purpose

The useful hierarchy is not “a model folder that contains one fit plus special toy/systematic files.” A model point may have one nominal fit and any number of other fit executions. Each fit output is the same broad artifact type, with explicit metadata describing its role and inputs:

```text
model config / model point
  ├── nominal fit
  ├── systematic fit: method=co_shift, value=0.5
  ├── systematic fit: method=co_shift, value=1.0
  ├── toy generation: design=bias, toy_index=42, seed=...
  │     └── fit of that toy
  └── toy generation: design=bias, toy_index=43, seed=...
        └── fit of that toy
```

Represent every actual fit output as `artifact_type=fit`, whether nominal, a systematic variation, or a fit to a generated toy. Put its effective fit profile, likelihood/model settings, input content IDs, implementation version, and seed in the derivation descriptor when they affect computation. Record `model_point_id`, `fit_role`, `input_artifact_id`, `systematic_method`, `variation`, `toy_batch_id`, and `toy_index` as association/task metadata as applicable; only computation-affecting values enter the artifact key. The catalog joins allow “all fits for this model point” to be queried without giving that model point ownership of a physical copy.

Keep genuinely different objects distinct: toy-generation output/design is not itself a fit; a systematic summary/delta is not the same object as the alternative fit trace that contributed to it. Both can link to their fit artifacts. 

**Consolidated storage for toys and systematics:** For toy studies and systematic scans, storing individual directories and traces for hundreds of fits causes severe inode and metadata overhead on shared HPC filesystems. Instead:
- **Resort to saving fit results only:** In the vast majority of toy and systematic studies, only summary posterior statistics (median, credible intervals, convergence metrics like $\hat{R}$ and ESS, pull/bias values) are required. These should be aggregated and stored in a **single consolidated file** per batch (e.g. Parquet, HDF5, or CSV).
- **Consolidated trace files when retained:** If full posterior chains must be retained for validation or diagnostic inspection, combine the traces into a **single multi-toy file** (e.g., NetCDF with a `toy_index` dimension or an HDF5 group hierarchy) rather than hundreds of individual `.nc` files.
- This preserves the ability to query individual toy/systematic results while keeping physical file counts small, eliminating inode pressure, and dramatically speeding up notebook loading.

### Keep IDs machine-stable and provide filesystem-readable views

The opaque ID should not be the only way to identify an artifact for a human. It remains the unambiguous key, but a person should be able to browse a study using familiar labels, filenames, and standard shell tools without first running an analysis Python session.

Keep one canonical payload directory per artifact, and generate a lightweight readable view that references it:

```text
<store>/
  artifacts/fit/<artifact-id>/
    fit.h5 -> blob <content-id>
    fit.nc -> blob <content-id>
    manifest.json
  views/<campaign-name>--<campaign-id>/
    model=tech_note_200/
      nominal/fit--<short-id> -> artifact <artifact-id>
      systematics/co_shift/0p5/fit--<short-id> -> artifact <artifact-id>
      toys/bias/toy_00042/fit--<short-id> -> artifact <artifact-id>
    artifacts.tsv
```

`artifacts/` contains the derivation manifest and ordinary filenames that point to immutable blobs; `views/` contains model/campaign aliases to artifact records. Both levels can use relative symlinks on the shared POSIX filesystem (or generated small reference files for views if links are unsuitable); neither copies payload bytes. A root-level TSV/CSV index and each artifact's `manifest.json` make the contents legible in a terminal/editor. The SQLite catalog can also be inspected with the `sqlite3` command-line tool. The readable labels are convenience aliases; the derivation manifest and content IDs remain authoritative. This combines stable identity with familiar `ls`, `less`, `jq`, and direct `h5py`/ROOT access.

This is a common separation of concerns in practice: MLflow groups runs and logs artifacts, DVC gives readable workspace paths that point into a content-addressed cache, and HEP tools make files/datasets addressable by names while retaining physical-file identities. Do not force a choice between human-readable labels and stable IDs; expose both.

### Use Python campaign definitions for dynamic study logic; pick orchestration separately

Campaign specification, workflow execution, and scientific configuration are three different decisions. Since model/campaign construction here includes non-Cartesian points and method-specific branches, define models and studies as typed Python objects/functions rather than trying to make YAML express dynamic logic. For example, a `ModelSpec` registry and a `Study` builder can resolve into explicit `ModelPoint`/`FitSpec` task records. Serialize the resolved request and task list to a versioned JSON campaign manifest before execution; that snapshot is the reviewable, replayable record, even if the source study was Python.

Keep general environment/runtime settings in the existing environment configuration. Keep physics defaults/configuration in the agreed scientific config layer, but avoid a second set of implicit model choices in path names. For production, require the campaign to name the model keys and fit/toy/systematic profiles it intends to use.

**Snakemake is a better fit than Hamilton if the choice is strictly between those two for the campaign/workflow layer.** Snakemake describes file-producing rules and dependencies, expands model/point wildcards, displays the DAG, and can resume/invalidate jobs. Its Python-backed Snakefile can consume a Python-generated list of model points and express fan-out/fan-in, while the model/study definitions stay in ordinary typed Python modules. A September 2026 ATLAS analysis report describes Snakemake used for a ROOT-ntuple-to-results chain and scans over seeds, regions, signal models, and injected strengths.

There is an important NERSC qualification. NERSC recommends Snakemake mainly for single-node execution and discourages its usual Slurm executor pattern because it can submit many individual jobs and burden the scheduler. That is a recommendation about Snakemake's execution mode; it does **not** mean the analysis must run as one multi-node MPI job. QPyM2's desired parallelism is many independent single-node fit tasks. A Slurm job array can submit these together, with each element requesting one node; Slurm runs as many concurrently as account/QOS availability permits, and a `%N` array limit can cap concurrency. NERSC also notes that array elements count as individual jobs for scheduling, so arrays reduce submission overhead but do not make the work invisible to scheduler limits. For long fits, the current array/packed-job approaches are a sound starting point; measure their throughput before changing them. The Snakemake Slurm plugin has a `--slurm-array-jobs` option, but array submission adds synchronization and still needs validation against NERSC guidance and the actual task fan-out. For now, have the QPyM2 planner emit the existing TSV/array jobs; later, Snakemake can provide DAG and resume behavior with either a local executor for small studies or an array-capable Slurm adapter proven in a representative pilot. NERSC's TaskFarmer is another pattern for a set of independent tasks inside one multi-node allocation, but its manager reserves a node, it is aimed at longer tasks, and NERSC says it is not suited for complex dependency graphs; it is not necessary for the first conversion.

Hamilton is a Python-native **micro dataflow** library: functions declare dependencies, and the driver computes requested values, with materializers to read/write external files and adapters for local/thread/Dask/Ray-type parallelism. Its own guidance says it does not provision compute or schedule long-running jobs. It can make internal transformations and histogram-building functions more composable/testable, but does not replace the campaign planner, artifact store, or SLURM task manager. A good layering is to call Hamilton (if useful) from inside a QPyM2 task, not to use it as the top-level campaign engine.

If a Python-native task system with direct HEP batch-system support is a stronger requirement than Snakemake's file-rule workflow, also evaluate **law** (Luigi Analysis Workflows). Law is explicitly built for HEP analysis, separates workflow code from run/storage/software locations, supports branches/targets and batch submission including Slurm, and has published HEP analysis examples. It is more infrastructure/convention than a small custom campaign planner, so pilot it against an actual QPyM2 study before adoption. Regardless of a later engine choice, the typed QPyM2 campaign definition and resolved task manifest should remain the scientific source of truth; workflow engines should consume or translate that plan rather than reimplement its model/dependency semantics.

### Execution topology for QPyM2

The first execution target should preserve the existing operating model: a fit task receives an explicit `FitSpec`, requests **one node** when submitted to SLURM, and can run concurrently with independent fit tasks on other nodes. Slurm decides how many are active based on available resources and project limits; generated arrays can optionally cap concurrency. The same resolved task records must also run locally without a scheduler—serially or through the existing local pool—even if a study takes much longer. A campaign may contain stages with dependencies (stage → inputs → fits → toy/systematic aggregation), but the planner, task functions, artifact IDs, and output formats should not depend on which executor runs them. The planner records dependencies and task IDs; execution adapters handle local runs or SLURM arrays. This lets the project first correct identity, storage, and path semantics without also changing scheduling behavior.

Later, Snakemake can encode those dependencies and generate the dry-run DAG. Snakemake's Slurm executor supports arrays (`--slurm-array-jobs`), but NERSC's guidance is cautious about its cluster executor and scheduler load, so validate job grouping, array size, submission rate, array concurrency, failure recovery, and one-node resource requests using a representative QPyM2 campaign before relying on it. MLflow can be added separately as a parent campaign run with child fit runs; it should log the QPyM2 artifact IDs and metadata rather than become a prerequisite for publishing a fit file. Neither future tool is needed for the initial clean conversion.

### MLflow is a tracking layer, not automatically the artifact store or campaign engine

MLflow maps naturally to *one tracking run per actual QPyM2 fit*: use an experiment for an analysis/campaign family, a parent run for one study execution, and child runs for nominal/systematic/toy fits. Record model/config/point IDs, input artifact IDs, fit role/variation, seed, source revision, status, runtime, and scalar fit diagnostics as params/tags/metrics. Its dataset-tracking API can record dataset source/digest/metadata. This gives the familiar ML experiment comparison UI and is a credible way to align with ML tooling.

Do not interpret MLflow's “model registry” as QPyM2's physics-model registry. In MLflow, a model generally means a trained/deployable model artifact; QPyM2's model is a statistical definition/configuration, and an individual Bayesian posterior fit is a run/result. Use QPyM2 `model_config_id` and `model_point_id` as explicit tags/parameters rather than forcing them into MLflow's registered-model abstraction.

MLflow's metadata backend and artifact store are separate. If the canonical QPyM2 artifact store already holds large ROOT/HDF5/NetCDF files, log their artifact IDs and paths/URIs in MLflow and log small diagnostics/reports there; do not upload a second copy of every large file per child run. Let MLflow own those payloads only if QPyM2 deliberately chooses MLflow's artifact store as the single canonical payload store and verifies the needed layout/access patterns. For a team tracking server, use its supported server/database setup; a shared SQLite file concurrently written directly by array tasks has the same NERSC filesystem/concurrency issue described above. Thus MLflow is worth piloting for a UI and experiment comparisons, but it should complement the QPyM2 artifact identity/reuse layer rather than replace it by default.

### ROOT/HEP has complementary tools, not one universal MLflow equivalent

ROOT is the analysis and I/O framework: ROOT files can hold self-describing objects and a RooWorkspace can persist a statistical model and its data. ROOT is not itself an experiment tracker, workflow scheduler, or general artifact catalog. In large experiments, those responsibilities are split across systems:

- **Rucio** is distributed scientific data management: names/scopes, files/datasets, metadata, replicas, transfers, and storage policy at ATLAS scale. It solves a much larger multi-site data-placement problem than QPyM2's local/shared-filesystem fit artifact index.
- **PanDA** is a workload/task management system used by ATLAS to distribute and monitor large-scale production and analysis jobs; it works alongside Rucio. It is not a per-analysis fit-result browser. CMS has experiment-specific workload tools such as CRAB. These are experiment infrastructure, not drop-in Python libraries for a small analysis repository.
- **law/Luigi and REANA** are closer analogues for analyst-authored workflows and reproducibility. Law is a Python HEP workflow library with batch/storage abstractions; REANA is a CERN analysis-preservation platform that packages inputs, code/environment, workflow, and outputs, and can run several workflow-description systems. REANA is a deployable platform, not just a local package.
- **Snakemake** is also now a demonstrated HEP option: the 2026 ATLAS anomaly-search workflow uses it across ROOT input processing and systematic/parameter scans. That is evidence of fit to particle-physics analysis structure, but does not override NERSC's local Slurm execution guidance.

For QPyM2, the practical comparison is therefore: an optional MLflow tracking server for searchable fit runs; QPyM2's own small artifact catalog for canonical objects/reuse; and a campaign/execution layer chosen for NERSC (custom planner + current arrays first, with Snakemake or law piloted). Rucio/PanDA/REANA become relevant only if the requirement expands to experiment-scale distributed data or an institutional reproducibility service.

## Patterns in neutrino, dark-matter, radio, and multi-messenger projects

These projects separate **scientific data identity**, **physical storage/replicas**, **processing lineage**, and **workflow execution** in different ways. The table focuses on patterns that transfer to QPyM2, not on adopting the full experiment infrastructure.

| Area / example | Data and metadata pattern | Workflow / analysis tools | What transfers to QPyM2 |
|---|---|---|---|
| **DUNE (neutrinos)** | DUNE separates MetaCat (logical files/datasets, searchable file and provenance metadata) from Rucio (logical identifiers, physical replicas, and placement). MetaCat queries can select files by metadata and parent/child relationships. | `justIN` workflows use MetaCat queries to select work and match jobs to data locations; it runs over the collaboration's HTCondor/glideinWMS infrastructure. | A data catalog and an artifact-location/replica manager are distinct concepts. For QPyM2, a small catalog can map an analysis artifact ID to one local path; Rucio-like federation is only justified if files need distributed placement. The task planner should carry logical inputs and dependencies instead of embedding every location in model names. |
| **IceCube (neutrinos)** | The I3 file format stores event frames and per-file/run detector metadata; file paths/names have historically encoded useful fields, while a separate file catalog was identified as necessary for checksums, status, provenance, and time/event ranges. | IceTray is the event-processing framework. IceProd handles organized production; user analyses have used HTCondor DAGMan. IceCube has evaluated Rucio for replication rather than treating ROOT/I3 files as the catalog. | Human-readable paths are valuable for browsing but do not replace a queryable catalog. Keep a readable view and a structured manifest/catalog together. Note that public IceCube Rucio material includes a 2018 evaluation, so it is evidence for the separation pattern, not a claim about today's complete IceCube stack. |
| **XENON1T/XENONnT (direct-detection dark matter)** | `strax` identifies a processed data product by run ID, data type, and its complete lineage: producer plugin and version, configuration values, and upstream lineage. Its storage frontend locates products while backends read/write chunks; it can cache and reuse compatible processed data. `straxen` supplies detector-specific plugins/contexts; `XeDocs` separately versions shared detector corrections and metadata. | Python plugin/context dataflow for event processing. XENON's `alea` is a separate Python inference package for statistical models, fitting, confidence intervals, sensitivity, and toy-MC studies. | This is the closest design analogue for QPyM2's staged/input products: artifact identity should include the exact upstream data and processing lineage, but not unrelated fit settings. Keep event processing/storage, physics/inference code, and shared detector/analysis metadata behind separate APIs. The fit output itself still needs a separate fit-run ID and role. |
| **LUX-ZEPLIN (dark matter)** | Raw data files are paired with metadata; the NERSC production workflow catalogs raw and reconstructed products in a MongoDB-backed Data Catalog accessed via REST. Successful processing validates ROOT output, registers it, and moves it to shared storage. | SPADE coordinates movement/warehousing; `spade-arrival` triggers processing/cataloging; PSQUARED plus RabbitMQ manages production job state; LZap performs reconstruction. | This is direct evidence from a NERSC-hosted dark-matter experiment for keeping files on shared storage and maintaining a separate catalog/API and job-state service. Its MongoDB, replicated data centers, and continuous production load are much larger than QPyM2 needs; borrow the boundary and post-write validation, not its service footprint. |
| **Radio astronomy (Stimela/CARACal, SKA, ASKAP/CASDA)** | Stimela recipes are human-readable YAML workflow definitions. SKA's data-product metadata package writes an adjacent `ska-data-product.yaml` with execution/processing block IDs, processing configuration, files, descriptions, status, size, CRC, and ObsCore fields. CASDA offers searchable science-ready FITS images/cubes/catalogues and CASA Measurement Sets, with descriptive product types/names and validation/processing records. | Stimela2 wraps radio tools as typed “cabs,” composes recipes, uses containers, and has Slurm/Kubernetes backends. SKA SDP uses processing scripts/workflows and product metadata. CASDA is a data discovery/archive service, not an individual analysis run manager. | Sidecar manifests are a proven, shell-readable way to describe complex scientific files and their processing. Keep descriptive names for product type/variation, add a stable ID and complete manifest, and keep the numerical payload in its native format. YAML is common for static workflow recipes in astronomy, but QPyM2 can use typed Python for dynamic model/campaign construction and emit a frozen JSON manifest. |
| **Gravitational-wave and multi-messenger astronomy** | IGWN separates primary-data discovery (`datafind`), event records (GraceDB), archival/distributed data (including Rucio and OSDF/Pelican), and pipeline-specific derived products. Bilby writes a self-contained result file (JSON by default or HDF5) with posterior and auxiliary fit information. SkyPortal stores GCN events/localizations and coordinates cross-messenger follow-up; it is an event/follow-up platform, not a batch-fit artifact store. | PyCBC/IGWN offline search workflows commonly use Pegasus and HTCondor DAGMan; parameter estimation tools include Bilby, LALInference, and others. Online alert ingestion uses streaming/event services such as GCN Kafka and SkyPortal. | Distinguish a **campaign/run** from an **event**, and both from the persistent fit result. Put event/model/fit/campaign IDs in metadata and link them explicitly. Use workflow/task engines for execution and catalogs for discovery; do not use an alert portal as a fit catalog. |

### Cross-domain implications

1. **The closest scientific-data precedent is `strax`.** Its lineage key makes “same run” insufficient for reuse: the producer version, config, and dependency chain determine whether a cached data type is compatible. QPyM2 should apply this at stage/input/fit boundaries and maintain different identities for each layer.
2. **DUNE and LZ reinforce catalog/file separation.** DUNE's MetaCat and Rucio split file metadata/querying from replica placement; LZ's Data Catalog indexes raw and reconstructed files while production jobs run elsewhere. QPyM2 can implement this boundary with JSON manifests plus SQLite because it has one configured filesystem root and no requirement to replicate data across institutions.
3. **Keep a manifest readable without the service.** SKA writes data-product metadata alongside the file bundle; IceCube and CASDA also use recognizable filenames/products to help browsing. QPyM2 should keep `manifest.json` with each object and generated readable views/TSV indices, while never treating the readable label as its identity.
4. **Workflow and data catalog are orthogonal.** `justIN`, Stimela, Pegasus, IceProd, and PSQUARED manage task execution or job state; MetaCat, Rucio, LZ's Data Catalog, and SKA/CASDA products handle data discovery, location, or metadata. QPyM2 should not expect adopting one workflow engine to solve artifact identification or deduplication.
5. **Use tools at their intended scale.** Rucio is now used or evaluated across DUNE, SBN/ICARUS, gravitational-wave collaborations, SKA, Rubin, and other communities, but its main value is multi-facility data placement and replication. For QPyM2, borrow logical IDs and metadata queries, but keep a local SQLite catalog until cross-site replication is a real requirement. Likewise, do not deploy production job services such as `justIN`, Pegasus, or PSQUARED just to organize one analysis's fit runs.

These examples strengthen the current proposal rather than changing its basic shape: use a filesystem and portable per-artifact manifests, keep human-readable names/views, and have a small catalog resolve logical artifacts to physical paths. The most useful patterns to borrow are strax-style lineage-sensitive cache keys, SKA-style sidecar metadata, and DUNE/LZ-style separation of data catalog from workflow control.

## Migration and implementation sequence

1. **Measure before moving.** Inventory artifact counts/sizes by stage and existing model folder. Sample repeated staging/input files by size and checksum to estimate true duplication and identify widely shared products. Do not run a full checksum sweep over all production data without measuring its I/O cost.
2. **Define the store contract.** Settle on model config ID, the canonical derivation-key descriptor, payload content hashing, manifest schema, completeness/verification rules, and one environment-configured store root. Specify which settings affect each stage so human model keys/full configs do not invalidate reusable upstream products. Reuse the dependency guidance in the workflow review rather than creating a second keying system.
3. **Implement the minimal artifact library.** Use `pathlib`, JSON manifests, and SQLite in the standard Python stack. Provide artifact publish/get/query/verify/ensure and an index-rebuild command. Keep the catalog rebuildable and payloads in their existing formats.
4. **Convert the current pipeline cleanly, without changing its execution options.** Keep the current command-line runners and local serial/pool execution as well as the TSV/SLURM submission path. Have the Python campaign planner produce explicit task records and arrays where each fit task requests one node; let Slurm schedule the allowed number concurrently. First route staging and input creation through the store and prove that differently named model points with identical effective input dependencies resolve to the same derivation IDs and shared blobs. In a small controlled test, remove one regenerable derived object and verify that `ensure` rebuilds it from its recorded parents and recipe.
5. **Make reusable histograms/templates addressable, then move fits.** Give arrays/histograms stable derivation keys and content IDs. Make fit-input manifests refer to shared observed-data/template products; where current fit readers require HDF5, create one bundle per unique composition. Write each actual fit to a unique fit artifact bundle with its role and parent/input IDs. Make toy/systematic summaries and batches reference the individual fits they consume. Preserve current numerical behavior while replacing reconstructed model-folder paths with catalog/artifact references.
6. **Provide human-readable access and index legacy data read-only.** Generate a campaign/model/fit-role view and TSV index. Create catalog records from existing folders and filenames, retaining original paths and treating inferred metadata as legacy/uncertain. Do not rename, delete, or relocate existing files in the first migration.
7. **Add future tool adapters only after the native contract works.** Add Snakemake to visualize/execute the resolved dependency graph and test its Slurm-array mode against a representative fan-out. Add MLflow parent/child run logging for campaign/fit metadata and metrics; avoid copying canonical large artifacts into a second store. Keep both optional and ensure the current runner can still use the core library directly.
8. **Benchmark and tune storage.** Compare ROOT and Parquet staging on event counts, histograms/systematics, file size, read/filter time, and actual shared-filesystem metadata behavior. Add checksums/stronger source identity according to measured cost and scientific requirements.
9. **Revisit database backend only with evidence.** If catalog indexing is small and performed by one collector, SQLite remains the minimal solution. If users require simultaneous cross-node updates or a live shared service, move the catalog backend behind the same small API to PostgreSQL or an approved service.

## Acceptance criteria

- A model is identified by an explicit key and resolved config identity, never by its output folder name alone.
- Artifact lookup returns the exact manifest and a direct file path without requiring the caller to parse a directory/filename convention.
- Repeated campaigns referencing an identical compatible stage/input reuse one artifact rather than produce copied payloads.
- Different model labels with identical data-producing dependencies resolve to the same derivation artifact and shared content blobs; model labels never determine payload paths or keys.
- Independent models can reference the same stored event arrays/histograms/templates, while a changed selection, binning, component template, or fit setting creates only the affected new product.
- Derivation IDs include only effective operation dependencies; output files are stored once per exact content hash, with provenance kept in manifests and association records.
- Changed dependencies/settings create distinct IDs; irrelevant downstream settings do not invalidate upstream data.
- An artifact can be validated using its manifest, schema, and file checksums; an incomplete or mismatched output is not silently treated as a cache hit.
- The catalog can be recreated from artifact manifests, and ordinary tools can read the underlying files without the QPyM2 catalog API.
- A missing regenerable blob can be materialized from recorded dependencies and an identified code/environment; missing source or environment dependencies are reported explicitly.
- Each independent fit task has a resource profile for one node; generated arrays can cap concurrency while preserving a stable mapping from task index to artifact ID.
- The same resolved campaign can run locally (serially or with a local pool) without requiring SLURM, a tracking server, or another cluster service.
- Many SLURM tasks can publish artifacts without concurrent writes to one shared SQLite database.
- Legacy results remain discoverable during gradual migration, and cleanup does not remove artifacts still referenced by retained results.

## Conclusion

The immediate recommendation is a **clean native conversion first**: a small derivation-key cache to reuse work before execution, a content-addressed blob store to keep one copy of each exact payload, reusable array/histogram products, fit-input manifests that compose references to those products, and `ensure()`-style rematerialization from recorded lineage when a derived blob is absent. Keep readable model/campaign views, a rebuildable SQLite index, typed Python model/campaign definitions, and executor-independent task functions. The same resolved tasks should run locally (serially or with a local pool) or through the existing TSV/SLURM array executor. Model names are aliases that point through many-to-many relationships to artifacts; they never appear in artifact or blob paths/keys. Each independent fit task requests one node when submitted, and multiple such tasks run in parallel as Slurm capacity and project limits allow. Add Snakemake later for DAG construction/resume and MLflow later for searchable campaign/fit tracking, after both can consume the same artifact/task contracts. The adjacent experiments reinforce the separation of data catalog, workflow control, and scientific analysis; they do not require QPyM2 to take on their distributed infrastructure.

## Sources

External documentation checked 2026-10-17. Links point to the project documentation as checked; details may evolve.

### ML/AI data and experiment management

- [DVC internal files and content-addressed cache](https://dvc.org/doc/user-guide/project-structure/internal-files)
- [DVC run cache](https://dvc.org/doc/user-guide/pipelines/run-cache)
- [DVC pipeline dependencies and parameter invalidation](https://dvc.org/doc/user-guide/pipelines/defining-pipelines)
- [MLflow backend stores](https://mlflow.org/docs/latest/tracking/backend-stores.html)
- [MLflow artifact stores](https://mlflow.org/docs/latest/tracking/artifact-stores.html)
- [MLflow tracking APIs and nested runs](https://mlflow.org/docs/latest/ml/tracking/tracking-api)
- [MLflow dataset tracking](https://mlflow.org/docs/latest/ml/dataset/)

### Workflow engines and HEP examples

- [Snakemake rules, wildcards, and dynamic input/output declarations](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
- [NERSC workflow tools overview](https://docs.nersc.gov/jobs/workflow/)
- [NERSC guidance for Snakemake](https://docs.nersc.gov/jobs/workflow/snakemake/)
- [NERSC TaskFarmer: task-list execution and restart behavior](https://docs.nersc.gov/jobs/workflow/taskfarmer/)
- [NERSC Slurm job-array example and guidance](https://docs.nersc.gov/jobs/examples/#job-arrays)
- [Slurm job arrays: one-node tasks and concurrency limits](https://slurm.schedmd.com/job_array.html)
- [Snakemake Slurm executor: job-array options](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [NERSC guidance for Parsl](https://docs.nersc.gov/jobs/workflow/parsl/)
- [Apache Hamilton dynamic DAG/parallel execution](https://hamilton.apache.org/concepts/parallel-task/)
- [Apache Hamilton: why use Hamilton and its orchestration scope](https://hamilton.apache.org/get-started/why-hamilton/)
- [law: Luigi Analysis Workflows](https://github.com/riga/law)
- [End-to-End Analysis Automation over Distributed Resources with Luigi Analysis Workflows (2024)](https://arxiv.org/abs/2402.17949)
- [ATLAS 2026: A Snakemake-Based Workflow for Reproducible Particle Physics Analyses](https://cds.cern.ch/record/2971347)
- [REANA workflow specification](https://docs.reana.io/reference/reana-yaml/)
- [REANA ROOT/RooFit analysis example](https://github.com/reanahub/reana-demo-root6-roofit)

### Neutrino and dark-matter experiments

- [DUNE Data Catalog / MetaCat documentation](https://dune.github.io/DataCatalogDocs/)
- [MetaCat documentation](https://fermitools.github.io/metacat)
- [MetaCat: metadata catalog for data management systems (Fermilab, 2021)](https://lss.fnal.gov/archive/2021/conf/fermilab-conf-21-061-scd.pdf)
- [justIN workflow command documentation](https://dunejustin.fnal.gov/docs/justin_command.man_page.md)
- [DUNE 2026 data-management status presentation](https://indico.fnal.gov/event/72820/contributions/341338/attachments/199647/278244/2026_data_handling_dpf_sct.pdf)
- [IceCube IceTray public software](https://github.com/icecube/icetray-public)
- [IceCube data-management / Rucio evaluation presentation (2018)](https://indico.cern.ch/event/676472/contributions/2905770/attachments/1609072/2554431/20180301_IceCube_Data_Management_-_Rucio_Workshop.pdf)
- [Strax storage API and lineage key](https://strax.readthedocs.io/en/latest/developer/storage.html)
- [Strax lineage compatibility and cache behavior](https://strax.readthedocs.io/en/latest/advanced/fuzzy_for.html)
- [Straxen XENONnT analysis framework](https://straxen.readthedocs.io/)
- [XeDocs XENON metadata management](https://github.com/XENONnT/xedocs)
- [Alea XENONnT toy-MC inference framework](https://github.com/XENONnT/alea)
- [LZ real-time data workflow at NERSC (2025)](https://doi.org/10.1051/epjconf/202533701232)
- [LZ Data Management Principles](https://lz.lbl.gov/wp-content/uploads/sites/6/2021/09/LZ-DMP-Operations.pdf)

### Radio and multi-messenger astronomy

- [Stimela2 documentation](https://stimela.readthedocs.io/)
- [Stimela2 framework paper](https://arxiv.org/abs/2412.10080)
- [SKA SDP Data Product Metadata](https://developer.skao.int/projects/ska-sdp-dataproduct-metadata/en/latest/)
- [SKA SDP data processing/workflow model](https://developer.skao.int/projects/ska-sdp-integration/en/latest/design/data-processing.html)
- [CASDA product types, filenames, and ancillary metadata](https://research.csiro.au/casda/data-products/)
- [Pegasus LIGO/PyCBC workflow example](https://pegasus.isi.edu/application-showcase/ligo)
- [IGWN data discovery, distribution, and processed data documentation](https://computing.docs.ligo.org/guide/data)
- [Bilby sampler output formats](https://bilby-dev.github.io/bilby/api/bilby.core.sampler.run_sampler.html)
- [SkyPortal multi-messenger event and localization model](https://skyportal.io/docs/mma.html)
- [Rucio cross-experiment community (2026)](https://indico.fnal.gov/event/73933/contributions/340014/attachments/197728/275065/rucio_for_wfms_workshop.pdf)

### ROOT and experiment-scale HEP data/workload systems

- [ROOT overview: storage, access, analysis, and Python bindings](https://root.cern.ch/about)
- [Rucio: scientific data management](https://rucio.cern.ch/)
- [Rucio documentation](https://rucio.cern.ch/documentation)
- [PanDA workload management system](https://pandawms.org/)

### Data formats and catalogs

- [Apache Parquet overview](https://parquet.apache.org/docs/overview/)
- [Zarr Python user guide](https://zarr.readthedocs.io/en/stable/user-guide/index.html)
- [PyIceberg documentation](https://py.iceberg.apache.org/)
- [Delta Lake overview](https://docs.databricks.com/aws/en/delta)

### SQLite and NERSC filesystems

- [SQLite: appropriate uses](https://sqlite.org/whentouse.html)
- [SQLite: WAL mode and network filesystem limitation](https://sqlite.org/wal.html)
- [NERSC filesystems overview](https://docs.nersc.gov/filesystems/)

### QPyM2 sources

- [`qpym2/cfg/model_utils.py`](../../../qpym2/cfg/model_utils.py): current model directory naming and resolution.
- [`qpym2/cfg/config_builders.py`](../../../qpym2/cfg/config_builders.py): current output path helpers.
- [`config/environments/tech_note.yaml`](../../../config/environments/tech_note.yaml): example split paths and model directory convention.
- [`docs/todos/analysis-workflow-and-config-reorganization.md`](../../todos/analysis-workflow-and-config-reorganization.md): campaign and artifact identity plan.
- [`docs/todos/staging-parquet.md`](../../todos/staging-parquet.md): existing proposal for backend-neutral staged-event storage.
- [`2026-10-15-workflow-and-model-testing-review.md`](2026-10-15-workflow-and-model-testing-review.md): related task/campaign and provenance recommendations.
