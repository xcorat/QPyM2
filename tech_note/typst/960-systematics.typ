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

// ── Systematic fit values ──────────────────────────────────────────────────
// Source of truth: output/syst_rate_results.md  (2026-05-31)
// Placeholder "?" marks values not yet extracted from fit output files.

// Nominal fit (reference: cs=1.0, gamma10_v150, u1740)
#let nom_raw_mode  = "11.70"
#let nom_raw_hdi_l = "4.95"
#let nom_raw_hdi_r = "20.38"
#let nom_bc_mode   = "15.26"
#let nom_bc_hdi_l  = "8.59"
#let nom_bc_hdi_r  = "23.81"

// BM MC statistical
#let bmmc_raw_mode   = "10.90"
#let bmmc_raw_hdi_l  = "3.89"
#let bmmc_raw_hdi_r  = "18.91"
#let bmmc_raw_dhdi_l = "-1.06"
#let bmmc_raw_dhdi_r = "-1.47"
#let bmmc_bc_mode    = "14.47"
#let bmmc_bc_hdi_l   = "7.55"
#let bmmc_bc_hdi_r   = "22.36"
#let bmmc_bc_dhdi_l  = "-1.04"
#let bmmc_bc_dhdi_r  = "-1.45"

// LSS ΔQ (±1σ)
#let lss_dq_raw_mode   = "12.13"
#let lss_dq_raw_hdi_l  = "3.08"
#let lss_dq_raw_hdi_r  = "23.37"
#let lss_dq_raw_dhdi_l = "-1.86"
#let lss_dq_raw_dhdi_r = "+3.00"
#let lss_dq_bc_mode    = "15.68"
#let lss_dq_bc_hdi_l   = "6.76"
#let lss_dq_bc_hdi_r   = "26.77"
#let lss_dq_bc_dhdi_l  = "-1.84"
#let lss_dq_bc_dhdi_r  = "+2.96"

// LSS σ-scaling (+1σ)
#let lss_sscl_raw_mode   = "12.41"
#let lss_sscl_raw_hdi_l  = "4.16"
#let lss_sscl_raw_hdi_r  = "21.09"
#let lss_sscl_raw_dhdi_l = "-0.78"
#let lss_sscl_raw_dhdi_r = "+0.71"
#let lss_sscl_bc_mode    = "15.96"
#let lss_sscl_bc_hdi_l   = "7.82"
#let lss_sscl_bc_hdi_r   = "24.51"
#let lss_sscl_bc_dhdi_l  = "-0.77"
#let lss_sscl_bc_dhdi_r  = "+0.70"

// Signal energy shift (Δ = 1.0 keV/MeV)
#let ndbdshift_raw_mode   = "10.62"
#let ndbdshift_raw_hdi_l  = "4.53"
#let ndbdshift_raw_hdi_r  = "18.86"
#let ndbdshift_raw_dhdi_l = "-0.42"
#let ndbdshift_raw_dhdi_r = "-1.51"
#let ndbdshift_bc_mode    = "14.19"
#let ndbdshift_bc_hdi_l   = "8.18"
#let ndbdshift_bc_hdi_r   = "22.32"
#let ndbdshift_bc_dhdi_l  = "-0.42"
#let ndbdshift_bc_dhdi_r  = "-1.49"

// Cobalt shift selection (likelihood-weighted marginalisation) — bc only
#let cobalt_bc_mode    = "15.83"
#let cobalt_bc_hdi_l   = "8.24"
#let cobalt_bc_hdi_r   = "23.54"
#let cobalt_bc_dhdi_l  = "-0.35"
#let cobalt_bc_dhdi_r  = "-0.27"

// Bias correction (central vs. covariance-sampled BC parameters)
// raw row = nominal raw (BC does not affect pre-correction posterior)
#let biascorr_raw_mode   = nom_raw_mode
#let biascorr_raw_hdi_l  = nom_raw_hdi_l
#let biascorr_raw_hdi_r  = nom_raw_hdi_r
#let biascorr_bc_mode    = "15.98"
#let biascorr_bc_hdi_l   = "8.43"
#let biascorr_bc_hdi_r   = "24.01"
#let biascorr_bc_dhdi_l  = "-0.16"
#let biascorr_bc_dhdi_r  = "+0.20"

// Marginalised posterior and split-Gaussian smearing widths
#let marg_mode    = "16.60"
#let marg_hdi_l   = "8.04"
#let marg_hdi_r   = "25.70"
#let marg_p90     = "29.81"
#let marg_sigma_p = "4.89"
#let marg_sigma_m = "3.98"

// ── syst_table helper ─────────────────────────────────────────────────────
// Two-row (non-corrected / bias-corr) table.
// Columns: Fit | Mode | Δ Mode | 68% HDI | Δ HDI
// d_mode:  shift in mode vs. nominal (shown on non-corrected row)
// d_hdi_l/r: shift in 68% HDI edges vs. nominal (shown on bias-corr row)
#let syst_table(
  raw_mode:    "?",
  raw_hdi_l:   "?",
  raw_hdi_r:   "?",
  d_mode:      "?",
  raw_d_hdi:   [—],
  d_bc_mode:   "?",
  bc_mode:     "?",
  bc_hdi_l:    "?",
  bc_hdi_r:    "?",
  d_hdi_l:     "?",
  d_hdi_r:     "?",
) = btable(
  columns: (auto, 1fr, 1fr, 2fr, 1.5fr),
  inset: 6pt,
  align: (left, center, center, center, center),
  table.header([*Fit*], [*Mode*], [*$Delta$ Mode*], [*68% HDI*], [*$Delta$ HDI*]),
  [non-corrected], [#raw_mode], [#d_mode],    [#raw_hdi_l, #raw_hdi_r], raw_d_hdi,
  [bias-corr],     [#bc_mode],  [#d_bc_mode], [#bc_hdi_l, #bc_hdi_r],  [#d_hdi_l, #d_hdi_r],
)
// ─────────────────────────────────────────────────────────────────────────

= Systematic Uncertainties <systematics>

This section evaluates the impact of systematic uncertainties on the extracted #ndbd event count. The evaluation metric, HDI-edge shift method, and marginalization procedure are described in @methods. The following sources are considered:

1. BM fit statistical errors.
2. Lineshape modeling (LSS).
3. Signal energy shift.
4. Cobalt shift selection.
5. Bias correction parameters.
6. PSA efficiency.

The reference model uses the fit from @results with cobalt shift $Delta = 1.0$ keV/MeV, `gamma10_v150` cobalt cut, and `u1740` signal cutoff.

== BM fit errors <sec:bm_fit_errors>

#figure(
  placement: auto,
  image("../images/06-bmmc_trace_combined_posterior_biascorr_cs1p0.png", width: 100%),
  caption: [Bias-corrected combined posterior vs. nominal ($Delta = 1.0 times 10^(-3)$).\ Dashed lines: mode; dotted: 68% HDI],
)<fig-bmmc>

We assess the impact of BM statistical fluctuations by repeating the M2 fit with 500 background spectra, each drawn from a random sample of the BM posterior (JAGS output). The per-sample posteriors are combined via likelihood-weighted mixture (@methods), yielding the combined posterior in @fig-bmmc. The combined posterior shifts left by approximately 1 event relative to the reference.


#align(center)[
  #syst_table(
    raw_mode: bmmc_raw_mode, raw_hdi_l: bmmc_raw_hdi_l, raw_hdi_r: bmmc_raw_hdi_r,
    d_mode: "-0.80", raw_d_hdi: [#bmmc_raw_dhdi_l, #bmmc_raw_dhdi_r],
    d_bc_mode: "-0.79",
    bc_mode: bmmc_bc_mode, bc_hdi_l: bmmc_bc_hdi_l, bc_hdi_r: bmmc_bc_hdi_r,
    d_hdi_l: bmmc_bc_dhdi_l, d_hdi_r: bmmc_bc_dhdi_r,
  )
]

== Lineshape Systematics <sec:lineshape_systematics>

The data and BM MC spectra are processed with lineshape scaling (LSS) corrections, which adjust the energy scale ($Delta Q$) and resolution ($sigma(E)$) per dataset. Miscalibration near the ROI can bias signal extraction. We evaluate this systematic by generating spectra with the LSS parameters perturbed by their $plus.minus 1 sigma$ uncertainties and refitting.

=== Lineshape functions and parameters <sec:lineshape_functions>
The lineshape model uses an asymmetric Gaussian core:

$ E_a = 1/(1 + A_L + A_R) ( G(E; mu, sigma) + A_L  G(E; mu, a_L, sigma) + A_R G(E; mu, a_R, sigma)) $

with per-dataset energy scale and resolution corrections:

$ Delta Q = mu' - mu = p_0 + p_1 E + p_2 E^2, quad sigma(E) = q_0 + q_1 E $

The asymmetry parameters ($A_L$, $A_R$, $a_L$, $a_R$) are fixed from the peakshape fit to the $ ""^(208)"Tl"$ 2615 keV peak; the scaling parameters ($arrow(p)$, $arrow(q)$) are derived from fits to other calibration peaks#footnote[TODO: add reference to the 1-tonne analysis for LSS parameter derivation]. Applying full-MC variations require recreating spectra for each sampled LSS point per input source #footnote[recreating spectral histograms is the most time consuming part of the fit, and done only once per fit in current setup]. Instead, the systematic is evaluated by applying only the $plus.minus 1 sigma$ ($+1 sigma$) shifts to LSS parameters, and re-calculating the staged event energies by applying a corrective $Delta Q$ ($sigma$ scaling) #footnote[When creating staging files, a 10 keV padding is added to account for these LSS variations] before creating the input spectra.

=== Systematic estimation <sec:lss_estimation>

For $Delta Q$, a $plus.minus 1 sigma$ shifts in covariant $p_0$ and $p_1$ (from the residual fit uncertainties) are applied to events in one dataset at a time, leaving other datasets unchanged. For $sigma$-scaling, each event is smeared by $delta sigma(E) = delta q_1 E dot sigma("ds", "ch")$ at only $+1 sigma $ (since we cannot sample a narrower gaussian).

Each perturbed fit is compared directly to the nominal; the 68% HDI edge shifts are measured and extremum shifts $Delta_plus.minus$ for 45 independent pertubations (15 datasets $times$ 3 perturbations: $sigma$-scaling, $2 times "dQ" "at" plus.minus 1sigma$) are recorded.
//where the sum runs over datasets and perturbation samples ($plus.minus sigma$ for $Delta Q$; $+1 sigma$ for $sigma$-scaling

// We also evaluated two alternative interval estimators — M1 (Only mode shift) and M3 (mode-anchored half-interval) — and found that while M1 gives the most conservative (largest) error estimates, it oversamples the statistical uncertainties; results for M1 and M3 are reported in the appendix.

#figure(
  placement: auto,
  grid(
    columns: 1,
    row-gutter: 0.5em,
    image("../images/lss_dq_perds_barplot.png", width: 90%),
    image("../images/lss_sscl_perds_barplot.png", width: 90%),
  ),
  caption: [LSS systematic per dataset (HDI interval shift, $Delta = 1.0 times 10^(-3)$). Each bar shows $Delta n_i$ for the $Delta Q$ (left/right) and $sigma$-scaling perturbations.],
  )

=== Final LSS-marginalized posterior <sec:lss_marginalized>

The combined effect of left-right bradening of the HDI interval due to LSS ($n_plus.minus$ per dQ and $sigma$) is calculated by combining the above $Delta_plus.minus$ values in quadrature per side per dataset as described in @methods.


  $ n_+ = sqrt(sum_("ds", plus.minus) (Delta^+_i)^2), quad n_- = sqrt(sum_("ds", plus.minus) (Delta^-_i)^2) $

The reference posterior is smeared with a split-Gaussian kernel of these widths, giving the LSS-marginalized posterior in @fig-lss-splitgaus-posterior.

#figure(
  image("../images/lss_syst_posteriors.png", width: 100%),
  caption: [LSS split-Gaussian marginalized posterior vs.\ bias-corrected reference (cs=1.0). Dashed: mode; dotted: 68% HDI;],
) <fig-lss-splitgaus-posterior>

#align(center)[
  #btable(
    columns: (auto, auto, 1fr, 1fr, 2fr, 1.5fr),
    inset: 6pt,
    align: (left, left, center, center, center, center),
    table.header([*Systematic*], [*Fit*], [*Mode*], [*$Delta$ Mode*], [*68% HDI*], [*$Delta$ HDI*]),
    table.cell(rowspan: 2)[$Delta Q$ ($plus.minus 1 sigma$)],
    [non-corrected], [#lss_dq_raw_mode],  [+0.43], [#lss_dq_raw_hdi_l, #lss_dq_raw_hdi_r], [#lss_dq_raw_dhdi_l, #lss_dq_raw_dhdi_r],
    [bias-corr],     [#lss_dq_bc_mode],   [+0.42], [#lss_dq_bc_hdi_l, #lss_dq_bc_hdi_r],  [#lss_dq_bc_dhdi_l, #lss_dq_bc_dhdi_r],
    table.cell(rowspan: 2)[$sigma$-scaling ($+1 sigma$)],
    [non-corrected], [#lss_sscl_raw_mode], [+0.71], [#lss_sscl_raw_hdi_l, #lss_sscl_raw_hdi_r], [#lss_sscl_raw_dhdi_l, #lss_sscl_raw_dhdi_r],
    [bias-corr],     [#lss_sscl_bc_mode],  [+0.70], [#lss_sscl_bc_hdi_l, #lss_sscl_bc_hdi_r],  [#lss_sscl_bc_dhdi_l, #lss_sscl_bc_dhdi_r],
  )
]

== Signal energy shift systematic <sec:signal_shift_syst>

A mismatch between the true $0 nu beta beta$ peak position and the signal template similar to the cobalt shift could bias the extracted event count. To account for the possibility that the signal might also be shifted due to the same mechanism as the #Co60 peaks, we refit with a signal spectrum shifted by $Delta = 1.0$ keV/MeV and record the 68% HDI edge displacement (@fig-ndbd-shift).
#figure(
  image("../images/06-ndbd_shift_posterior_biascorr_cs1p0.png", width: 100%),
  caption: [Bias-corrected posterior with shifted $0 nu beta beta$ signal (blue) vs.\ nominal (red).]
) <fig-ndbd-shift>

#align(center)[
  #syst_table(
    raw_mode: ndbdshift_raw_mode, raw_hdi_l: ndbdshift_raw_hdi_l, raw_hdi_r: ndbdshift_raw_hdi_r,
    d_mode: "-1.09", raw_d_hdi: [#ndbdshift_raw_dhdi_l, #ndbdshift_raw_dhdi_r],
    d_bc_mode: "-1.07",
    bc_mode: ndbdshift_bc_mode, bc_hdi_l: ndbdshift_bc_hdi_l, bc_hdi_r: ndbdshift_bc_hdi_r,
    d_hdi_l: ndbdshift_bc_dhdi_l, d_hdi_r: ndbdshift_bc_dhdi_r,
  )
] <tbl-ndbd-shift-stats>

// NOTE: we only need this if we are using ndbd uncorrected counts for systematics. using bias-corrected 
// errors standardize all the errors prior to PSA.
// == Marginalized Posterior for $n_(0 nu beta beta)$

// Studying the systematic effets from the Bias correction and PSA efficiency is done using the posterior distribution marginalized across the three studies above. We combine the posteriors from these three studies assuming our posteriors are close to a split-gaussian distribution, and shift in mode for each study are combined in quadrature (seperately for left and right shifts) to obtain the left and right widths of the gaussian broadening of the reference posterior. 


== Cobalt shift selection <sec:systematics_co60>

The cobalt shift $Delta$ is fixed to $1.0$ keV/MeV in the reference model, but its optimal value is uncertain. We treat this choice as a systematic by forming a likelihood-weighted mixture posterior over all seven tested shift values (weights shown in @fig-comarg), following the procedure in @methods, and comparing it to the reference.

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
    grid(
      columns: (1fr, 1fr),
      column-gutter: 0.5em,
      image("../images/comarg_ref_vs_comarg_bc.png", width: 100%),
      image("../images/comarg_per_shift_bc_posteriors.png", width: 85%),
    ),
  ),
  caption: [Cobalt shift systematic. Left: likelihood-weighted combined bc posterior vs.\ nominal ($Delta = 1.0$ keV/MeV). Right: bc posteriors for each tested shift value.]
) <fig-comarg>

The combined posterior shifts left by approximately $0.5$ events. We take the average of the left and right 68% HDI displacements ($-0.51$ events) as the systematic bias from cobalt shift selection.

// No raw posterior available for cobalt shift (likelihood-weighted marginalisation only).
// bc: mode=#cobalt_bc_mode, HDI=[#cobalt_bc_hdi_l, #cobalt_bc_hdi_r], Δ=[#cobalt_bc_dhdi_l, #cobalt_bc_dhdi_r]

== Bias Correction Systematics <sec:bias_correction_syst>

Low-statistics posteriors carry a reconstruction bias that is corrected before quoting results (TODO: add cross-reference to bias correction derivation). The correction parameters carry statistical uncertainties from calibration fits; we propagate these by sampling them from their covariance matrix and reapplying the correction.

#figure(
  placement: auto,
  grid(
    columns: 1,
    row-gutter: 1em,
    align(center)[
      #syst_table(
        raw_mode: biascorr_raw_mode, raw_hdi_l: biascorr_raw_hdi_l, raw_hdi_r: biascorr_raw_hdi_r,
        d_mode: "+0.00",
        d_bc_mode: "+0.73",
        bc_mode: biascorr_bc_mode, bc_hdi_l: biascorr_bc_hdi_l, bc_hdi_r: biascorr_bc_hdi_r,
        d_hdi_l: biascorr_bc_dhdi_l, d_hdi_r: biascorr_bc_dhdi_r,
      )
    ],
    image("../images/bias_corr_comparison.png"),
  ),
  caption: [Reference (central bias-correction parameters) vs.\ systematic (bias-correction parameters sampled from their covariance) bias-corrected posterior. All values in bias-corrected event counts. The systematic uncertainty is taken as $[+Delta sigma_+, -Delta sigma_-] = [-0.78, -1.00]$ events.]
) <fig-biascorrsyst>

The varied posterior shifts by $+0.87$ events in mode, with 68% HDI edges moving by $-0.13$ (lower) and $+0.09$ (upper). We take $[+Delta sigma_+, -Delta sigma_-] = [-0.78, -1.00]$ events as the systematic uncertainty (@fig-biascorrsyst).

== PSA Systematics <sec:psa_systematics>

Uncertainties in PSA efficiency curves can affect the fit. Ares simulations account for efficiency by discarding events probabilistically. The exposure-averaged signal efficiency is 98%. The error on efficiency is 0.12%, and this is propagated directly to the rate calculation from the marginalized bias corrected posterior that include all the other systematics.

== Summary of Results <sec:systematics_summary>

@tbl-syst-summary lists all systematic contributions to the $0 nu beta beta$ event count, evaluated as signed 68% HDI edge shifts relative to the reference fit ($Delta = 1.0$ keV/MeV, `gamma10_v150` cut, `u1740` signal cutoff). All values are in bias-corrected event counts; right ($+$) and left ($-$) contributions are combined in quadrature independently.

#academic-table(
  [Summary of systematic uncertainties on the $0 nu beta beta$ event count (bias-corrected). All values in events. The TOTAL row gives the quadrature sum of each side independently.],
  header: ([*Source*], [*$+sigma$ (right)*], [*$-sigma$ (left)*]),
  (
    [LSS — $Delta Q$ ($plus.minus 1 sigma$)],              [$+4.74$], [$-$3.59],
    [LSS — $sigma$-scaling ($+1 sigma$)],                  [$+1.17$], [$-1.24$],
    [BM MC statistical],                                   [$+0.00$], [$-1.04$],
    [$0 nu beta beta$ spectrum shift ($1.0$ keV/MeV)],     [$+0.00$], [$-0.42$],
    [Cobalt shift selection],                              [$+0.00$], [$-0.35$],
    [Bias correction],                                     [$+0.20$], [$-0.16$],
    table.hline(),
    [*TOTAL (quad-sum)*],                                  [*$+4.89$*], [*$-3.98$*],
  ),
  columns: (2.2fr, 1fr, 1fr),
  align: (left, center, center),
  label: <tbl-syst-summary>,
)

PSA efficiency uncertainty (0.12% on $epsilon_"PSA" = 0.98$) is propagated at the rate-conversion stage and contributes a negligible $tilde.op 0.1%$ relative shift, subdominant to the statistical uncertainty.

The total uncertainty is dominated by lineshape modeling: LSS $Delta Q$ alone contributes $+4.7$ and $-3.6$ events. Bias correction and BM MC contributions are subdominant; cobalt and signal shift terms are negligible. No single systematic variation qualitatively reshapes the posterior.

The bias-corrected reference posterior has mode $= #nom_bc_mode$ events. After marginalising over all sources via split-Gaussian smearing (widths $sigma^+ = #marg_sigma_p$, $sigma^- = #marg_sigma_m$ events), the marginalised posterior has mode $= #marg_mode$ and $p_90 = #marg_p90$ events. Rate conversion and the final result are presented in @conclusion.

