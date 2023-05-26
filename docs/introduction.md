# Introduction

In this analysis we aim to use Multiiplicity Two events to recover some of CUORE exposure to $0\nu\beta\beta$ in TeO2. 



Multiplicity Two (M2) events register in two different bolometer channels (Ch1, Ch2) as energies (E1, E2) and taken as in-coincidence with 150mm and 30ms cuts.

This analysis aim to provide a pathway into searching for $0vbb$ decay in CUORE multi-site events. The main CUORE analysis focuses on multiplicity 1 events which account for ~88% of the expected exposure, while this  analysis uses events that deposit energy in two bolometers (multiplicity 2 or M2, events) to recover sensitivity from about 4-5% of events that are expected to be M2-like. The region of interest (ROI) for this analysis spans multiple features in the expected background, and we use the CUORE background model to inform us about the expected background. We use simulated data for reconstruct our expected backgorund and signal (0vbb) event distributions, and extract the maximum likelihood rate of 0v for our data by varying the strengths of each background component. The minimization is done using the Baysian MC tool PyMC3. 

## Motivation

While the share of exposure that could be recovered is around 5%, M2 analysis benefits from much lower background levels (~2e-3 ckky) than the M1 analysis, thus the expected sensitivity is ~6-8% of the main analysis.

## Setup

+ ROI: 
+ Rotation:
+ Binning:
+ Smoothing:

+ Fit: 4 bkg components and 0v
+ Data from Nt2021
+ Background


This analysis is a search for $0\nu\beta\beta$ decay in ${}^{130}\mathbf{Te}$ with a selection of CUORE data where two coincident energy depositions occur near the region of interest of 0vbb.