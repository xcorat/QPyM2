#import "common.typ": *
#import "@preview/akatable:0.1.0": academic-table

= Data <data>

The core inputs for this fit are simulated spectra for the expected background and #ndbd signal, alongside observed data. These are equivalent to the BM input files for the 1 ton-yr data release.

To improve performance, we generate "staging files" containing only relevant M2 events in the fit region. These staging files serve as the primary inputs for the analysis.

== Dataset Details

The fit uses 15 datasets (ds3601-ds3615), totaling $1038.7 "kg.yr"$#footnote[Reading and summing exposures from PeakShape files, not 1038.4? check!] of $"TeO"_2$ exposure from the first 4 years of CUORE operation. Data production mirrors the main #ndbd analysis. As detailed in 1-ton-yr BM #footnote[Section 2, Internal Note 137D], data is pre-processed to calculate pulse parameters, then we apply filtering, gain-stabilizion #footnote[heater or calibration, add details?], and finally calculate the calibrated energies.

Consistent with the BM and low-energy analyses, we use a 30 ms coincidence window, and events must occur in crystals less than 150 mm apart 
(center-to-center) to be classified as M2.

We apply the standard main analysis cuts, including basic pulse quality, pileup rejection, PCA, and line-shape validity checks.

#academic-table(
  [Analysis cuts applied to the M2 dataset.],
  header: ([*Cut*], [*Variable*], [*Description*]),
  (
    [Signal Events],        [`IsSignal`],                                       [Removes noise and pulser triggers],
    [Pileup rejection],     [`NumberOfPulses==1`, `NoHeaterInWindow`],          [],
    [Reject bad intervals], [`RejectBadIntervals`],                             [],
    [Bad for Analysis],     [`BadForAnalysis`],                                 [],
    [PCA],                  [`PCA`],                                            [],
    [Included],             [`Included`],                                       [Further analysis cuts],
    [Multiplet validation], [`FullMultiplet && Validation`],                    [Rejects coincident pairs if either event fails a preceding cut],
  ),
  columns: (auto, 1fr, 1fr),
) #footnote[Update this]

After these cuts, the overall selection efficiency is 2.51% (M2 selection 2.57%, PCA 97.7%) #footnote[Present this better? different numbers?].

== Simulations

*Background*: Background contaminant sources are simulated using `geant4` based `qshields` which includes the experiment's geometry and materials to accurately model the interactions. The deposited energies from decays are then processed using `Ares` to apply the detector response (including the livetimes) to obtain the expected background spectra for each of the sources. These files are copied directly from the BM model input, and include the contaminants in the TeO2 crystals, the copper parts from the cryostat, shields and tower assembly, as well as cosmic muons and a contribution from Silver metastable state #Ag110m. Altogether, there are 84 components in the
background model #footnote[any components that are used only for systematic testing?]. A detailed description of the BM components can be found in the internal note 137D.

*Signal*: We generate $10^8$ #ndbd events using qshields and apply detector response using the same approach as the background components.

#figure(
  image("../images/02-expected_spectra.png", width: 80%),
  caption: "Expected Spectra"
) #footnote[Add z axis units/label]

== Data

The observed data is processed with same parameters as above, which follows the CUORE low energy analysis#footnote[is this correct?], and the selected values from the processing tree are saved in a ROOT file (`ReducedData.root`) and copied directly from the BM fit input. The events from this file is then filtered according to our selection cuts before being fed into the scripts that create staging files for the fit.

== Data pipeline

Since only few percent of all the CUORE exposure is used in this analysis, we create input files with only the events needed for the M2 analysis. We save these files as ROOT files with a minimum set of branches, which acts as the basic input for the analyses. These ROOT are referred to as *Staging files*#footnote[staging file generation at:], and the selection for data and simulations indicated above are applied while creating these files. 

While looking for the optimal analysis model, we tested multiple event selection criteria with different fit models. Since the fit uses binned background spectra (instead of events) in addition to observed data, we create histograms using the staging files depending on the fit model and save them as HDF5 *input files*#footnote[input generation] along with the data as an event list. As the fit and the rest of the analysis is based on python/numpy, these histograms and data as saved as numpy arrays.

