#import "common.typ": *
#import "@preview/akatable:0.1.0": academic-table

= Data <data>

This analysis utilizes our observed data and the corresponding Monte Carlo simulations for the expected background and #ndbd signal equivelant to the CUORE 1 ton-yr Background Model (BM) data release. This includes 15 datasets and 84 background source simulations in addition to the simulated signal. 

To optimize computationally heavy procedures within the Bayesian framework, intermediate staging files are generated containing only the relevant M2 events within the region of interest, and actual fit inputs are all converted to numpy compatible formats.

== Experimental Data <sec:experimental_data>

The analysis incorporates 15 datasets (ds3601-ds3615), constituting an exposure of $1038.4 "kg" dot "yr"$ of $"TeO"_2$ from the first four years of CUORE operations (out of $1122.9 "kg" dot "yr"$ exposure)#footnote[Reading and summing exposures from PeakShape files. Check if 1038.4 or 1038.7 matches the final official value.]. Data production parallels the standard single-site #ndbd analysis framework. The raw data is pre-processed to extract pulse parameters, followed by filtering, thermal gain stabilization, and energy calibration.

Consistent with the 1 ton-yr BM and low-energy analyses, a 30 ms temporal coincidence window is applied. Events are classified as Multiplicity 2 (M2) if they trigger exactly two crystals with a center-to-center distance of less than 150 mm. Please refer to the Internal Note 137D for a detailed review of the Background Model and the data processing steps#footnote[which technote details the data processing steps for 1ton-yr?].

Standard data quality criteria are enforced, including basic pulse quality validation, pileup rejection, Principal Component Analysis (PCA) profiling and lineshape quality. Both events in an M2 pair must independently pass all selection criteria; otherwise, the entire multiplet is rejected#footnote[Check if this is the case? I forgot if we decided to use filtering only on the main event]

#align(center)[
#box(width: 90%, [
#academic-table(
  [Analysis cuts applied to the M2 dataset.],
  header: ([*Cut*], [*Variable*], [*Description*]),
  (
    [Signal Events],        [`IsSignal`],                                       [Removes noise and pulser triggers],
    [Pileup rejection],     [`NumberOfPulses==1`, `NoHeaterInWindow`],          [Rejects closely spaced pulses],
    [Reject bad intervals], [`RejectBadIntervals`],                             [Excludes periods of known detector instability],
    [Bad for Analysis],     [`BadForAnalysis`],                                 [Removes poorly performing channels],
    [PCA],                  [`PCA`],                                            [Pulse shape discrimination],
    [Included],             [`Included`],                                       [Standard analysis inclusion flag (lineshape)],
    [Multiplet validation], [`FullMultiplet && Validation`],                    [Rejects coincident pairs if either event fails a preceding cut],
  ),
  columns: (auto, 1fr, 1.5fr),
  label: <tbl-analysis-cuts>,
)
])
]

After these cuts, the overall selection efficiency is 2.51% (M2 selection 2.57%, PCA 97.7%).

== Simulations <sec:simulations>

*Background*: The expected background spectra are derived from Monte Carlo simulations using the `geant4`-based `qshields` framework, which precisely reconstructs the experiment's geometry and material compositions. Energy deposits are processed using the `Ares` software to simulate the detector response, incorporating specific channel livetimes. These simulated spectra mirror the standard BM inputs and include contaminant sources from the $"TeO"_2$ crystals, the copper cryostat inner components, the shields, the tower assembly, cosmic muons, and the metastable $""^(110m)$ Ag state, totalling 84 independent physical components.

*Signal*: The nominal ndbd signal template is generated from $10^8$ simulated decays using `qshields`, with the detector response applied identically to the background components.

#figure(
  image("../images/02-expected_spectra.png", width: 80%),
  caption: [Expected Spectra (TODO: Add Z-axis units/label)]
)

== Data Pipeline <sec:data_pipeline>

Because this analysis utilizes only a small subset of the total CUORE exposure, we reduce the full experimental data into ROOT files containing a minimal set of branches. These *staging files* apply the aforementioned selection criteria for both data and simulations and serve as the foundational inputs. 

Given that the Bayesian fit evaluates binned background spectra alongside event-level data, the framework generates optimized binned histograms from the staging files depending on the selected model configuration. These are saved in HDF5 format as *input files*, while the event lists are retained. Since the analytical framework (`qpym2`) is implemented in Python, saving these structures natively enables rapid memory loading as `numpy` arrays for iterative model testing and evaluations.

