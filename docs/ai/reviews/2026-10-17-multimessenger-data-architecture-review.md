# Multimessenger astronomy data architecture and distributed systems review

- **Review date:** 2026-10-17
- **Scope:** Multimessenger astronomy data lifecycles, event-alert infrastructure, large-scale distributed processing, heterogeneous data formats and standards, and tools in use across major astronomy systems.
- **QPyM2 context:** Relates field-level architecture patterns to the existing CUORE M2 analysis framework. This review complements [`2026-10-17-data-management-review.md`](2026-10-17-data-management-review.md), which focuses on QPyM2 artifact identity, storage, and discovery.
- **Research note:** External service and project status is time-sensitive; links point to documentation checked on the review date.

## Executive summary

Multimessenger astronomy is a **federated system**, not a single data lake. Observatories operate their own detectors, reduction pipelines, archives, schemas, access policies, and analysis facilities. Shared alert services and community brokers provide rapid coordination; standards-based APIs help researchers discover and retrieve data from separately operated archives.

The dominant architectural separation is between:

1. **The alert/event plane:** low-latency, compact, machine-readable notices distributed to many consumers.
2. **The science-data plane:** large images, spectra, strain series, event tables, and analysis products retained in observatory archives and data facilities.
3. **The processing plane:** local or distributed batch and streaming compute, workflows, data replication, and user analysis environments.
4. **The interoperability plane:** schemas, identifiers, semantics, and APIs that let heterogeneous services work together without forcing identical internal storage.

Kafka is increasingly used for alerts, but it is not a substitute for long-term product archives. Brokers turn very high-volume survey streams into filtered, cross-matched, classified, and prioritized products. Large-scale processing uses separate tools for workflow scheduling, data movement/replication, file abstraction, distributed catalog queries, and chunked-array analysis.

QPyM2 currently follows a file-oriented HPC batch pattern—ROOT staging, HDF5 fit inputs, PyMC fits, NetCDF/HDF5 outputs, and SLURM task arrays. It is not an alert-stream consumer or a cross-observatory data federation. The most transferable ideas are explicit data contracts, revision-aware event identity if event feeds are ever added, durable provenance, separation of metadata from payloads, and execution near large datasets. Kafka, Rucio, Qserv, or a distributed catalog are not implied requirements for the current QPyM2 workload.

## 1. Multimessenger data and why integration is difficult

The term *multimessenger* covers distinct measurement systems, not merely different file formats:

| Messenger / domain | Typical products | Important context |
|---|---|---|
| Electromagnetic: gamma ray through radio | Images, spectra, light curves, source catalogs, exposure and quality maps | Instrument response, bandpass, calibration, detection limits, and observing cadence |
| Gravitational waves | Detector strain time series, candidate parameters, sky-localization probability maps, data-quality segments | Large sampled time series and compact candidate notices; parameter estimates can be revised as analyses improve |
| Neutrinos | Event direction and energy estimates, event quality, astrophysical signalness, skymaps | Directional uncertainty and sub-threshold events can be central to a joint search |
| Cosmic rays | Arrival direction, energy, shower/reconstruction products, exposure | Reconstruction and exposure vary by detector and event class |

A cross-messenger association is usually **probabilistic**, not a direct join on a shared event ID. A search asks whether signals overlap in time and sky position, given each instrument's uncertainty, background, exposure, and selection function. AMON is a representative system: it searches event streams for spatial and temporal coincidences, including combinations of sub-threshold events, and issues alerts when a correlation is significant. [AMON overview](https://www.amon.psu.edu/) · [AMON system](https://www.amon.psu.edu/amon-system)

At minimum, a useful event record needs more than a coordinate and timestamp. It needs a time scale, coordinate frame, uncertainty or probability representation, measurement units, data-quality context, event version/state, and references to detailed products.

## 2. Reference architecture: alert plane and science-data plane

```text
OBSERVATORY / DETECTOR
  raw measurements
       │
       ├── low-latency processing ──► ALERT / EVENT PLANE
       │                                GCN or SCiMMA Kafka topics
       │                                          │
       │                                          ▼
       │                              brokers / coincidence services
       │                           validate · cross-match · classify
       │                                          │
       │                             event record / archive / API
       │                                          │
       │                              telescope and follow-up systems
       │
       └── calibration and reprocessing ──► SCIENCE-DATA PLANE
                                             files / object storage
                                             catalogs / archive databases
                                                       │
                                         discovery and access APIs
                                                       │
                                          analysis near the data
```

### Alert/event plane

This path prioritizes latency, fan-out, message size, and machine readability. A message typically communicates what was detected, when and where, the current significance or confidence, the notice's relationship to previous notices, and how to locate supporting products.

- **NASA GCN** now recommends Kafka for automated alert consumption. Its documentation lists the legacy VOEvent Transport Protocol brokers as retired on **April 6, 2026**. New GCN notice types are available through its new Kafka service using JSON schemas; legacy formats remain available on specified channels. [GCN Notices](https://gcn.nasa.gov/docs/notices)
- **LVK public alerts** are published through GCN on the Kafka topic `igwn.gwalert`. Alert states include early warning, preliminary, initial, update, and retraction. Sky-localization and source estimates can change over hours or days, so consumers must handle revisions to an event rather than treating every notice as an unrelated event. [GCN LVK mission page](https://gcn.nasa.gov/missions/lvk)
- **SCiMMA Hopskotch** provides Kafka-based publish/subscribe, public and access-controlled topics, an archive REST API, and large-message offload. Its documented Kafka message retention is about 28 days; archived messages are retained separately. Oversized payloads are placed in the archive and represented on the topic by a small pointer message. [SCiMMA Hopskotch](https://www.scimma.org/hopskotch)
- **Event services** add durable event records, notification, authorization, and coincidence logic. GraceDB provides programmatic access and notifications for gravitational-wave candidates; AMON performs cross-stream correlation. [GraceDB](https://gracedb.ligo.org/) · [AMON](https://www.amon.psu.edu/)

Kafka provides scalable stream delivery and consumer decoupling. It should not be conflated with a mission archive, full science-product catalog, or the durable system of record for every downstream scientific result.

### Science-data plane

Large science products remain in observatory-controlled archives or data facilities. Alerts point to these products rather than carrying all raw measurements. Archive APIs typically expose product metadata and permit selective retrieval.

GWOSC, for example, exposes observation runs, event catalogs, parameter estimates, data-quality information, and strain files via a read-only REST API. Its API v2 includes pagination and request throttling; downloadable strain products include formats such as HDF5 and GWF. [GWOSC API](https://gwosc.org/api)

Rubin demonstrates the value of an intermediate broker tier for a survey-scale alert stream. Its alert packets are distributed to brokers rather than directly to individual users, and the same information is queryable through the Prompt Products Database and Rubin Science Platform. Rubin's broker ecosystem includes full-stream brokers and specialized downstream brokers. [Rubin alerts and brokers](https://rubinobservatory.org/for-scientists/data-products/alerts-and-brokers) · [Recent Rubin data releases](https://rubinobservatory.org/for-scientists/data-products/recent-data-releases)

## 3. Major systems and tools

| Layer | Examples | Function |
|---|---|---|
| Alert transport | Apache Kafka, NASA GCN, SCiMMA Hopskotch | Publish and fan out live notices to independent consumers |
| Alert brokers | Rubin community brokers, Fink, ALeRCE, ANTARES, Lasair, Babamul | Filter, enrich, cross-match, classify, prioritize, and expose tailored alert streams and APIs |
| Coincidence / event state | AMON, GraceDB | Correlate events; maintain candidate records, revisions, and access controls |
| Follow-up coordination | HERMES, Target-and-Observation Managers, telescope schedulers | Turn a notice into an observing request; coordinate observations and reporting |
| Science archives and user platforms | GWOSC, Rubin Science Platform, mission archives | Discover and retrieve data products, catalogs, and analysis context |
| Distributed catalog query | Rubin Qserv | Execute spatially partitioned SQL queries over very large catalogs |
| Workflow / job management | PanDA/iDDS, SLURM, Pegasus/HTCondor, Nextflow | Schedule dependency graphs or high-throughput tasks, manage retries, and distribute compute |
| Data management / movement | Rucio, FTS, XRootD, WebDAV, S3-compatible storage | Track files and replicas, enforce placement policies, and move data between facilities |
| Chunked array processing | Dask + xarray; Zarr, NetCDF, HDF5 | Work lazily and in parallel with large arrays and time-dependent data |

### Rubin as an end-to-end example

Rubin is a useful case study because it combines low-latency alert production with large-scale annual processing and distributed data access. Its architecture uses Kafka in prompt processing, container orchestration including Kubernetes/Knative, object storage, the Butler abstraction for locating and interpreting datasets, an alert/prompt-products database, and Qserv for catalog queries. PanDA manages workloads across data facilities; Rucio and transfer services coordinate file replication. [Rubin software architecture paper](https://arxiv.org/abs/2211.13611) · [Butler and pipeline execution](https://dmtn-229.lsst.io/) · [PanDA integration paper](https://arxiv.org/abs/2312.04921) · [Qserv documentation](https://qserv.lsst.io/) · [Rucio](https://rucio.cern.ch/)

The tools are complementary rather than interchangeable:

- **Workflow manager:** chooses which computation runs where and manages dependencies.
- **Data management:** tracks file identity, location, replicas, and transfer policy.
- **Data abstraction:** shields pipeline code from physical paths and storage details.
- **Catalog database:** supports low-latency structured queries over catalogs.
- **Alert bus:** transports time-sensitive event messages.

Rucio is appropriate when an experiment must manage large numbers of files and replicas across facilities. Qserv is specialized for petascale catalog queries. Neither is a general-purpose replacement for a filesystem, a workflow engine, or an alert service.

## 4. Formats and interoperability standards

There is no single format suitable for all astronomical products. Interoperability relies on descriptive metadata, schemas, and access protocols as much as on file formats.

| Product | Common representation / access pattern | Notes |
|---|---|---|
| Live notices | JSON for new GCN notices; Avro for high-volume survey alert packets; VOEvent XML for legacy and supported pathways | Keep schema/version and event revision explicit. Avro efficiently carries nested alert records; GCN's newer unified JSON schemas define shared notice fields. |
| Images and spectra | FITS; VOTable for tabular exchange | Alerts generally link to full products; they may embed or refer to cutouts and smaller supporting payloads. |
| Gravitational-wave strain | HDF5 and GWF, with JSON metadata and APIs | Keep detector, sample rate, time interval, quality flags, and product release/version attached. |
| Sky localizations | HEALPix probability maps, commonly FITS | Represent a sky probability distribution, not just a central coordinate and radius. |
| Sky/time coverage | IVOA MOC | MOC is for fast coverage representation/comparison; it is distinct from a probabilistic event localization map. |
| Catalogs and event tables | Relational databases, TAP/ADQL, VOTable | Large catalogs can be spatially partitioned internally while retaining standard query interfaces. |
| Large N-D arrays / time series | HDF5, NetCDF, and Zarr | Chunking and layout should match expected access patterns; xarray with Dask enables lazy, parallel computation. |

IVOA standards make independently operated data services more interoperable without requiring identical databases or internal file layouts:

- **TAP + ADQL:** service-based access to tables and catalog metadata, including spatial queries and asynchronous jobs. TAP describes multi-service query orchestration as a way to perform distributed cross-matches. [IVOA TAP](https://www.ivoa.net/documents/TAP)
- **ObsCore and DataLink:** common observation-product metadata and links from discovery records to data products or related resources. [IVOA DataLink](https://www.ivoa.net/documents/DataLink)
- **SIA / SSA and VOTable:** discovery/access patterns for image and spectral products, and a standard tabular response format.
- **MOC and HiPS:** interoperable coverage maps and progressive sky-survey visualization. [IVOA MOC](https://www.ivoa.net/documents/MOC) · [IVOA HiPS](https://www.ivoa.net/documents/HiPS)
- **VOEvent:** an IVOA model for event content, separate from any transport protocol. It remains relevant to existing notices and event descriptions, even though the legacy GCN VOEvent Transport Protocol is retired. [VOEvent 2.1 document](https://www.ivoa.net/documents/VOEvent/20250703/PR-VOEvent-2.1-20250703.html)

These standards are best understood as **adapters and contracts** over local systems, not a mandate that every observatory store all data in the same way. For large numerical arrays, Dask/xarray offers a computation model: data can be chunked, processed lazily, and executed across a cluster. It does not replace astronomy-specific coordinate, unit, or archive metadata. [xarray and Dask](https://docs.xarray.dev/en/stable/user-guide/dask.html)

## 5. Heterogeneity and event-data contract

The recurring integration challenges are:

1. **Time semantics:** timestamps require a defined scale/reference (for example UTC or GPS), precision, and conversion rules.
2. **Coordinate semantics:** a direction needs a coordinate frame and uncertainty. Broad or multimodal localizations should be represented as probability maps, not reduced to an unjustified point.
3. **Measurement semantics:** flux, energy, significance, and signalness depend on detector calibration and analysis definitions. Units, uncertainties, quality flags, and definitions matter.
4. **Revision semantics:** notices can improve or retract earlier candidates. Consumers must preserve event identity and state/version history.
5. **Schema evolution:** common fields should be stable while mission-specific data remains extensible and typed. Schema version compatibility must be explicit.
6. **Access policy:** some streams/products are public, private, embargoed, or collaboration-restricted. Authorization and provenance are part of the architecture.
7. **Selection effects:** joint statistical interpretation requires each instrument's background, exposure, and event-selection properties, not merely an apparent temporal/spatial match.

A compact event record should therefore include:

- Publisher, instrument, stable event identifier, and notice revision/state.
- Event time, time scale, precision, and any relevant interval.
- Coordinates/frame and a localization product reference when applicable.
- Significance, false-alarm rate, signalness, or probability with a defined meaning.
- Measurements, units, uncertainties, and quality flags.
- Schema/version and links to authoritative science products.
- Access policy and provenance: which pipeline/version created or revised the notice.

The message should carry what consumers need for triage and routing. Large data products should generally remain retrievable from their archive or be offloaded behind a resolvable reference.

## 6. Distributed architecture patterns and trade-offs

### Move compute toward data when data movement dominates

Moving petabyte-scale data to each user's laptop is neither practical nor economical. Data facilities and science platforms place compute near archives; remote users submit queries or run notebooks in managed environments. For large transfers between processing centers, explicit replica tracking and policy-driven transfer systems are used.

### Keep low-latency processing and reprocessing distinct

Prompt pipelines optimize latency and bounded processing of new measurements. Periodic data-release pipelines can consume much more compute, apply improved calibrations, and produce internally consistent reprocessed products. An alert should retain a link to its source/revision context even as later data releases supersede or refine the underlying products.

### Partition according to query and compute shape

- Catalog queries are partitioned spatially (Qserv is an example) to distribute work across database nodes.
- Time series and N-D arrays are chunked along dimensions used by computation (Dask/xarray/Zarr patterns).
- Independent parameter studies are split into tasks and scheduled using job arrays or workflow engines.
- Data replication is handled separately from task scheduling.

Poor partitioning, skewed work, very small files, expensive cross-node shuffles, and mismatched array chunking can dominate performance even when the compute framework scales horizontally.

### Maintain provenance through asynchronous processing

Event messages, follow-up observations, revised localizations, calibration products, and offline data releases are created by separate systems over different timelines. Stable IDs and links between records are necessary to reconstruct the chain later. Durable metadata should record which product revision was used in a scientific result.

## 7. Application to QPyM2

The repository's current pipeline is:

```text
raw ROOT
  → staged ROOT (filtered events and transformed u/v coordinates)
  → HDF5 fit inputs (data, MC histograms, model tables, bin edges)
  → PyMC Bayesian fits
  → ArviZ NetCDF trace + HDF5 components
  → systematics / toys / sensitivity / bias results
```

The project uses ROOT `RDataFrame` for event processing, pandas/h5py for HDF5, PyMC/ArviZ for inference, Hydra/OmegaConf for configuration, and TSV parameter grids executed via SLURM arrays on Perlmutter. This is documented in [`README.md`](../../../README.md) and implemented in [`qpym2/io/rdf.py`](../../../qpym2/io/rdf.py), [`qpym2/io/hdf.py`](../../../qpym2/io/hdf.py), [`qpym2/run/run_single_fit.py`](../../../qpym2/run/run_single_fit.py), and [`scripts/run_job.sh`](../../../scripts/run_job.sh).

This is **high-throughput batch analysis**, not a distributed alert-processing or archive-federation system. Independent fit points run as separate tasks; they consume file-based stage products and produce file-based results. That is a sound pattern for the current analysis and is compatible with the existing artifact-store direction.

Transferable improvements and boundaries:

- **Use explicit manifests for reproducibility:** record raw input identity, relevant resolved config, code revision, software environment, grid/task identity, random seed, outputs/checksums, and scheduler status.
- **Make stage contracts explicit:** describe units, transformations, histogram edges, and schema version alongside HDF5 content; use file metadata and sidecars to preserve interpretation.
- **Keep retries safe:** distinguish complete outputs from partial task writes and preserve enough task status to retry failed SLURM array elements deterministically.
- **Retain current formats unless a measured access pattern motivates a change:** ROOT remains suitable for current event processing; HDF5/NetCDF fit products remain readable by existing tools.
- **If event feeds are ever integrated:** retain notice ID/revision and product references; don't place full science products in a message queue.
- **Adopt IVOA services if external discovery/access becomes a requirement:** TAP/ADQL and ObsCore can make public observational metadata discoverable. They do not need to be introduced for an internal model-fit pipeline.
- **Defer distributed infrastructure until its requirement exists:** Kafka addresses live fan-out, Rucio multi-site placement, and Qserv petascale catalog queries. They solve problems QPyM2 does not currently show.

For the detailed QPyM2 artifact store, campaign, provenance, concurrency, and workflow recommendations, see [`2026-10-17-data-management-review.md`](2026-10-17-data-management-review.md) and [`2026-10-15-workflow-and-model-testing-review.md`](2026-10-15-workflow-and-model-testing-review.md).

## Sources

### Alerts, brokers, and multimessenger services

- [NASA GCN Notices: protocols, formats, schemas, and VOEvent Transport Protocol retirement](https://gcn.nasa.gov/docs/notices)
- [GCN LIGO/Virgo/KAGRA mission page: Kafka topic, notice types, and latencies](https://gcn.nasa.gov/missions/lvk)
- [GCN IceCube mission page and alert archive formats](https://gcn.nasa.gov/missions/icecube)
- [SCiMMA Hopskotch: Kafka topics, retention, archive, large-message offload](https://www.scimma.org/hopskotch)
- [SCiMMA: services, HERMES, and Rubin cross-match service](https://scimma.org/)
- [GraceDB candidate event database](https://gracedb.ligo.org/)
- [AMON overview](https://www.amon.psu.edu/)
- [AMON system architecture](https://www.amon.psu.edu/amon-system)
- [Rubin alerts and brokers](https://rubinobservatory.org/for-scientists/data-products/alerts-and-brokers)
- [Rubin recent data releases / Early Science alert status](https://rubinobservatory.org/for-scientists/data-products/recent-data-releases)
- [The Zwicky Transient Facility Alert Distribution System](https://arxiv.org/abs/1902.02227)
- [Fink broker architecture paper](https://arxiv.org/abs/2009.10185)

### Archives, distributed processing, and storage

- [GWOSC API](https://gwosc.org/api)
- [Rubin Observatory software architecture and system design](https://arxiv.org/abs/2211.13611)
- [Rubin Data Butler and Pipeline Execution System](https://dmtn-229.lsst.io/)
- [PanDA integration with Rubin Observatory](https://arxiv.org/abs/2312.04921)
- [Qserv: distributed petascale catalog database documentation](https://qserv.lsst.io/)
- [Rubin data movement model](https://arxiv.org/abs/2507.12955)
- [Rucio scientific data management](https://rucio.cern.ch/)
- [xarray: parallel computing with Dask](https://docs.xarray.dev/en/stable/user-guide/dask.html)
- [Nextflow workflow documentation](https://nextflow.io/docs/latest/)
- [Pegasus workflow documentation](https://pegasus.isi.edu/documentation)

### IVOA standards

- [IVOA standards and recommendations index](https://www.ivoa.net/documents)
- [Table Access Protocol (TAP)](https://www.ivoa.net/documents/TAP)
- [DataLink](https://www.ivoa.net/documents/DataLink)
- [Multi-Order Coverage map (MOC)](https://www.ivoa.net/documents/MOC)
- [Hierarchical Progressive Survey (HiPS)](https://www.ivoa.net/documents/HiPS)
- [Sky Event Reporting Metadata (VOEvent) 2.1](https://www.ivoa.net/documents/VOEvent/20250703/PR-VOEvent-2.1-20250703.html)
- [FAIR standards for astronomical data](https://arxiv.org/abs/2203.10710)

### QPyM2 sources

- [`README.md`](../../../README.md): current pipeline and execution overview.
- [`qpym2/io/rdf.py`](../../../qpym2/io/rdf.py): ROOT RDataFrame input and staging output.
- [`qpym2/io/hdf.py`](../../../qpym2/io/hdf.py): HDF5 fit input/output structure.
- [`qpym2/run/run_single_fit.py`](../../../qpym2/run/run_single_fit.py): fit-input loading, PyMC fit, and output writing.
- [`scripts/run_job.sh`](../../../scripts/run_job.sh): TSV task selection and SLURM task execution.
- [`2026-10-17-data-management-review.md`](2026-10-17-data-management-review.md): QPyM2 artifact storage and discovery review.
- [`2026-10-15-workflow-and-model-testing-review.md`](2026-10-15-workflow-and-model-testing-review.md): QPyM2 campaign, workflow, provenance, and testing review.
