#import "common.typ": *
#import "@preview/akatable:0.1.0": academic-table

#let btable(..args) = block(stroke: (bottom: 1.5pt))[
  #set table(
    stroke: (x, y) => (
      left: none,
      right: none,
      top: if y == 0 { 1.5pt } else if y == 1 { 0.75pt } else { none },
      bottom: none,
    ),
    inset: (x: 8pt, y: 5pt),
    align: left,
  )
  #table(..args)
]

= Systematic Uncertainties

This section details the systematic uncertainties in the M2 fit, including contributions from the CUORE background model (BM), data production (lineshape errors), and model selection. For most systematic effects, the shape of the component spectra changes when merging the individual MC spectra. Since recreating these components and fitting at each and every grid point is computationally prohibitive, we adopt an approximate approach of sampling the systematic variations around the reference fit for each possible contribution independently, and estimate the individual and combined effect of these variations on the edges of the final posterior distribution#footnote[margins of? whats a good wording without mentioning 68% HDI?] of the #ndbd event count.

The main sources of systematic uncertainty are:

1. BM fit statistical errors.
2. PSA efficiency.
3. Line-shape modeling in BM simulations.
4. Cobalt shift selection.
5. #ndbd spectrum shift.  

The model shown in section 5 with $Delta = 1.0$ keV/MeV Cobalt shift is used as the reference.

== Estimating Systematic Uncertainties

For each source of systematic uncertainty, we estimate its impact by refitting the data with a model with the varied parameters, and obtaining the posterior for bias corrected #ndbd number of events to compare with the reference fit. All studies other than background model systematics use one or two points in the systematic parameter space. For these posteriors, we record the shift in the left and right margins of the 68% HDI around mode from the reference fit (for LSS dQ we sample two points, and take the maximum shift of the two), and combine them in quadrature separately for left and right shifts, to estimate the widening of our reference posterior due to systematic effects. Then we simulate left-right smearing separately using a split-Gaussian to sample the original reference posterior and obtain a simulated marginalized posterior on #ndbd events. We use this posterior to obtain the posterior on #ndbd rate, from which we cite the left/right errors, as well as the limits.

We tested multiple interval estimators for quantifying the systematic shift, including the shift in mode, but chose to report the left and right shifts of the 68% high density interval as asymmetric systematic uncertainties.

== BM fit errors 

#figure(
  placement: auto,
  image("../images/06-bmmc_trace_combined_posterior_biascorr_cs1p0.png", width: 100%),
  caption: [Bias-corrected combined posterior vs. nominal ($"cs" = 1.0$). Dashed lines: mode; dotted: p16/p84; dash-dotted: p90.],
)<fig-bmmc>

We assess the impact of BM statistical fluctuations by repeating the M2 fit with background spectra generated from 500 random samples of the BM posterior (JAGS output). For each sample, we perform the M2 fit, and combine the posteriors across all samples, (TODO: weighted by their maximum likelihoods) to obtain a distribution of the $0 nu beta beta$ rate. We record a leftward shift of the $n_(#ndbd)$ of around 1 event (@fig-bmmc).
// The resulting systematic uncertainty is derived from the spread of this distribution compared to the referece posterior.


#academic-table(
  [BM MC systematic: mode and p90 shifts (raw and bias-corrected). All values in events.],
  header: ([*BMMC syst*], [*Nominal*], [*Syst (BMMC)*], [*$Delta$ mode*], [*$Delta$ 68% HDI*]),
  (
    [mode], [15.25], [14.47], [-0.79], [-1.044, +0],
    // [p90],  [24.56], [23.18], [-1.38], [-1.41],
  ),
  columns: (1.5fr, 1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center, center),
)

== Lineshape Systematics

The data and the BM MC spectra used here have been processed with the same lineshape scaling (LSS) corrections as the main 1 tonne-yr analysis. To account for the variation of LSS parameters, we rerun the fits with varied samples of these parameters according to their statistical fluctuations in the BM fit.

=== Lineshape functions and parameters
The LSS functions and parameters are defined as follows:

$ E_a = 1/(1 + A_L + A_R) ( "gaus"(E; mu, sigma) + A_L "gaus"(E; mu, a_L, sigma) + A_R "gaus"(E; mu, a_R, sigma)) $

and the scaling functions:

$ Delta  Q = mu' - mu = (p_0 + p_1E + p_2 E^2) $
$ sigma(E) = q_0 + q_1E $

The parameters $A_L$, $A_R$, $a_L$, $a_R$, and $sigma$ are determined from the peakshape fit#footnote[to Tl-208 2615 peak and saved into `PeakShape_dsNNNN.root files`], and the scaling functions and parameters ($arrow(p), arrow(q)$) are saved as polynomials in energy#footnote[Scaling functions are calculated using the fits to other peaks and the parameters are taken from `residual_and_width_vs_energy_dsNNNN.root` ]. However, trying to apply variations for all events in all the MCs is computationally expensive. Instead, we simulate the errors by applying only a $1 sigma$ variation to the parameters from ares output spectra (staging files).

=== Systematic estimation: 

For $Delta Q$ systematics, we apply a shift equivalent to $plus.minus 1 sigma$ variation in the $p_0$ and $p_1$ parameters to all the events (per dataset) by the error of the residual fit. The resulting change in the $0 nu beta beta$ mode defines the systematic error for that dataset; we run the systematic on each dataset while keeping the events from the rest of the datasets untouched.

For the resolution ($sigma$-scaling) systematic, we smear each event by a Gaussian draw from $delta sigma(E) = delta q_1 E dot sigma("ds", "ch")$ at $+1 sigma$.


Rather than combining perturbed posteriors into a mixture, each perturbed fit is compared directly to the nominal. Here too, the systematic for each dataset is taken as the shift in the left and right edges of the 68% HDI around mode of the $0 nu beta beta$ posterior (M2 estimator).

#figure(
  placement: auto,
  grid(
    columns: 1,
    row-gutter: 0.5em,
    image("../images/lss_dq_perds_barplot.png", width: 90%),
    image("../images/lss_sscl_perds_barplot.png", width: 90%),
  ),
  caption: [LSS systematic per dataset (HDI interval shift, $"cs" = 1.0$). Each bar shows $Delta n_i$ for the $Delta Q$ (left/right) and $sigma$-scaling perturbations.],
)

// We also evaluated two alternative interval estimators — M1 (Only mode shift) and M3 (mode-anchored half-interval) — and found that while M1 gives the most conservative (largest) error estimates, it oversamples the statistical uncertainties; results for M1 and M3 are reported in the appendix.

We record the maximum shift left of the left edge and right of the right edge, (or zero if an  edge of the systematic posterior is shifted towards the mode compared to the reference posterior) and add these per-dataset left/right shifts of the HDI in quadrature to obtain the combined LSS uncertainty.

  $ n_+ = sqrt(sum_("ds", plus.minus) (Delta^+)_i^2); quad n_- = sqrt(sum_("ds", plus.minus) (Delta^-_i)^2); quad Delta^plus.minus_i < 0$
  
where the indices run over the datasets and samples; $plus.minus sigma$ for $Delta Q$ and just $+1 sigma$ for $sigma$-scaling.

// #align(center,
  #table(
    columns: (auto, auto, auto),
    align: (left, right, right),
    table.header([*Source*], [*$n_+$*], [*$n_-$*]),
    [$Delta Q: plus.minus 1 sigma$], [+4.70], [$-$3.59],
    [$sigma$-scaling: $+1 sigma$],   [+1.17], [$-$1.26],
  )
// )

=== Final LSS-marginalized posterior

The bias-corrected reference posterior is smeared with a split-Normal kernel whose widths $sigma^+_"LSS"$ and $sigma^-_"LSS"$ are derived from the 45 per-dataset mode shifts (15 datasets $times$ 3 perturbations: $sigma$-scaling, dQ $+1sigma$, dQ $-1sigma$), combined in quadrature per side. This gives the final LSS-marginalized posterior shown in @fig-lss-splitgaus-posterior.

#figure(
  image("../images/lss_syst_posteriors.png", width: 100%),
  caption: [LSS split-Gaussian marginalized posterior vs.\ bias-corrected reference (cs=1.0). Dashed: mode; dotted: 68% HDI;],
) <fig-lss-splitgaus-posterior>

// #academic-table(
//   [LSS systematic: mode and p90 shifts. All values in events (bias-corrected).],
//   header: ([*LSS*], [*Reference (BC)*], [*LSS-marginalized (BC)*], [*$Delta$*]),
//   (
//     [mode], [15.55], [13.63], [-1.92],
//     [p90],  [27.49], [27.76], [+0.27],
//   ),
//   columns: (1.5fr, 1fr, 1fr, 1fr),
//   label: <tbl-lss-splitgaus-stats>,
// )

== Systematics due to possible shift in $0 nu beta beta$ decay

We also evaluate potential systematic bias from a shift in the $0 nu beta beta$ spectrum, analogous to the Cobalt shift. Using a signal spectrum shifted by $Delta = 1.0$ keV/MeV, we repeat the fitting procedure, and record the shift in the posterior HDI.
#figure(
  image("../images/06-ndbd_shift_posterior_biascorr_cs1p0.png", width: 100%),
  caption: [Bias-corrected posterior with shifted $0 nu beta beta$ signal (blue) vs.\ nominal (red).]
) <fig-ndbd-shift>

#academic-table(
  [Summary of $0 nu beta beta$ shift systematic: mode and 90% upper limit for nominal and shifted fits, raw and bias-corrected.],
  header: ([*NDBD shift*], [*Nominal*], [*Syst (shift)*], [*$Delta$ mode*], [*$Delta$ 68% HDI*]),
  (
    [mode], [15.25], [14.19], [-1.07], [-.42, +0],
    // [p90],  [24.55], [22.53], [-2.03], [-1.97],
  ),
  columns: (1.5fr, 1fr, 1fr, 1fr, 1fr),
  label: <tbl-ndbd-shift-stats>,
)

// NOTE: we only need this if we are using ndbd uncorrected counts for systematics. using bias-corrected 
// errors standardize all the errors prior to PSA.
// == Marginalized Posterior for $n_(0 nu beta beta)$

// Studying the systematic effets from the Bias correction and PSA efficiency is done using the posterior distribution marginalized across the three studies above. We combine the posteriors from these three studies assuming our posteriors are close to a split-gaussian distribution, and shift in mode for each study are combined in quadrature (seperately for left and right shifts) to obtain the left and right widths of the gaussian broadening of the reference posterior. 


== Cobalt shift selection 

// The combined posterior is wider than the reference, with a mode shift of *XX* and a 90% upper limit of *XX*.
#figure(
  placement: auto,
  grid(
    columns: 1,
    row-gutter: 1em,
    align(center)[
      #btable(
        columns: (auto, auto, auto, auto, auto, auto, auto, auto),
        inset: 6pt,
        align: center,
        table.header(
          [*$Delta$ 0.0*], [*0.3*], [*0.5*], [*0.8*], [*1.0*], [*1.2*], [*1.5*], [*(keV/MeV)*],
        ),
        [0.000], [0.003], [0.024], [0.530], [0.309], [0.044], [0.090], [*weight*],
      )
    ],
    image("../images/comarg_ref_vs_comarg_bc.png"),
  ),
  caption: [Nominal ($Delta = 1.0$ keV/MeV) vs.\ co-shift combined bias-corrected posterior. Left: per-shift likelihood weights used for the combination.]
) <fig-comarg>

We consider the marginalization over the Cobalt shift as a source of systematic uncertainty by comparing the width of the combined posterior over all shifts to that of the reference fit (with $Delta = 1.0$ keV/MeV). The combined posterior is obtained via adding the bias corrected posterior for each of the shifts weighted by their mean likelihood, which is then compared to the bias corrected reference posterior (@fig-comarg). 

The combined posterior is shifted to the left by approximately half an event compared to the reference as evident by the shifts in left and right 68% HDI edges, and we cite the average of the two ($-0.51$ events) as the systematic bias due to the Cobalt shift selection.

== Bias Correction Systematics

#figure(
  placement: auto,
  grid(
    columns: 1,
    row-gutter: 1em,
    align(center)[
      #btable(
        columns: (auto, auto, auto, auto, auto),
        inset: 6pt,
        align: center,
        table.header(
          [*Quantity*], [*Reference*], [*Systematic*], [*$Delta$ Mode*], [*$Delta$ 68% HDI*],
        ),
        [mode],    [15.10], [15.97], [+0.87], [-0.11, +0.17],
        // [p16],     [ 9.74], [ 9.61], [$-$0.13],
        // [p84],     [25.17], [25.26], [+0.09],
        // [p90 UL],  [27.79], [27.92], [+0.13],
        // [$Delta sigma_+$], [10.08], [9.29], [$-$0.78],
        // [$Delta sigma_-$], [ 5.36], [6.36], [+1.00],
      )
    ],
    image("../images/bias_corr_comparison.png"),
  ),
  caption: [Reference (central bias-correction parameters) vs.\ systematic (bias-correction parameters sampled from their covariance) bias-corrected posterior. All values in bias-corrected event counts. The systematic uncertainty is taken as $[+Delta sigma_+, -Delta sigma_-] = [-0.78, -1.00]$ events.]
) <fig-biascorrsyst>

To evaluate how the results change when applying the bias correction, we generate a posterior distribution by applying the bias correction while sampling its parameters with their covariant uncertainties (@fig-biascorrsyst). The mode of the systematically varied posterior shifts by $+0.87$ events, while the 68% interval edges move by $-0.13$ (lower) and $+0.09$ (upper), indicating a small asymmetric broadening. We take the change in the split-Gaussian widths, $[+Delta sigma_+, -Delta sigma_-] = [-0.78, -1.00]$ events, as the systematic uncertainty due to the bias correction procedure.

== PSA Systematics

Uncertainties in PSA efficiency curves can affect the fit. Ares simulations account for efficiency by discarding events probabilistically. The exposure-averaged signal efficiency is 98%. The error on efficiency is 0.12%, and this is propagated directly to the rate calculation from the marginalized bias corrected posterior that include all the other systematics.

== Summary of Results

@tbl-syst-summary lists all systematic contributions to the $0 nu beta beta$ event count, each evaluated as a signed shift in the bias-corrected posterior mode relative to the reference fit ($Delta = 1.0$ keV/MeV, `gamma10_v150` cut, `u1740` signal cutoff). All values are in bias-corrected event counts. Right ($+$) and left ($-$) uncertainties are combined in quadrature separately to give the total split-Gaussian smearing widths.

#academic-table(
  [Summary of systematic uncertainties on the $0 nu beta beta$ event count (bias-corrected). All values in events. The TOTAL row gives the quadrature sum of each side independently.],
  header: ([*Source*], [*$+sigma$ (right)*], [*$-sigma$ (left)*]),
  (
    [LSS — $Delta Q$ ($plus.minus 1 sigma$)],              [$+4.74$], [$-$3.59],
    [LSS — $sigma$-scaling ($+1 sigma$)],                  [$+1.17$], [$-1.24$],
    [BM MC statistical],                                   [$+0.00$], [$-1.04$],
    [$0 nu beta beta$ spectrum shift ($1.0$ keV/MeV)],     [$+0.00$], [$-0.42$],
    [Cobalt shift selection],                              [$+0.00$], [$-0.28$],
    [Bias correction],                                     [$+0.17$], [$-0.10$],
    table.hline(),
    [*TOTAL (quad-sum)*],                                  [*$+4.89$*], [*$-3.97$*],
  ),
  columns: (2.2fr, 1fr, 1fr),
  align: (left, center, center),
  label: <tbl-syst-summary>,
)

PSA efficiency uncertainty (0.12% on $epsilon_"PSA" = 0.98$) is propagated independently at the rate-conversion stage and contributes a negligible $tilde.op 0.1%$ relative shift to the rate, subdominant to the statistical uncertainty.

The bias-corrected reference posterior has mode $= 15.3$ and $p_{90} = 27.9$ events. After marginalising over all systematic contributions via split-Gaussian KDE smearing (widths $sigma^+ = 6.9$, $sigma^- = 5.3$ events), the marginalised posterior has mode $= 17.0$ and $p_{90} = 31.3$ events (@fig-final-rate).

Following rate conversion with exposure $= 1038.7$ kg$dot$yr, isotopic fraction $eta = 0.342$, $N_("TeO"_2) = 3.773 times 10^{24} "mol"^{-1}$, signal efficiency $epsilon_"MC"$, and PSA efficiency $epsilon_"PSA"$, the final marginalised rate posterior gives:

$ Gamma_(0 nu) = 3.5 times 10^(-25) "yr"^(-1) quad ("mode") $
$ Gamma_(0 nu) < 6.33 times 10^(-25) "yr"^(-1) quad (90% "C.L.") $

which translates to a lower limit on the half-life of:

$ T_(1/2)^(0 nu) > 1.09 times 10^24 "yr" quad (90% "C.L.") $

#figure(
  placement: auto,
  image("../images/06_syst_rate_new_final_rate_posterior.png", width: 70%),
  caption: [Final $0 nu beta beta$ decay rate posterior, marginalised over all systematic uncertainties and PSA efficiency. Red: mode; green dashed: 90% upper limit; shaded band: 68% HDI.],
) <fig-final-rate>
