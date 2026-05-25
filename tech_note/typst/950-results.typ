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

// --- Fit results: reference fit (co_shift=1.0, gamma10_v150, u1740) ---
// Source: 05-results-plots.ipynb, cell 8 (summary_tbl)
// Columns: mode | prior_mu (expected) | prior_sigma | p90
#let fit_co60   = (mode: "318.0", exp: "328.5", std: "77.4",  p90: "343.1")
#let fit_u238   = (mode: "313.5", exp: "332.6", std: "234.8", p90: "348.2")
#let fit_th232  = (mode: "343.5", exp: "631.4", std: "332.6", p90: "392.7")
#let fit_others = (mode: "435.6", exp: "119.6", std: "142.3", p90: "486.9")
#let fit_ndbd   = (mode: "12.8",  exp: "--",    std: "--",     p90: "25.2")

= Fit Results <results>

We construct our final posterior by selecting the model and combining the posteriors from each cobalt shift fit. The output of each fit is independently bias corrected and then combined with weights proportional to their maximum likelihoods. In this section, we discuss one reference fit in detail, followed by bias correction procedure and the final combined posterior after marginalization over the cobalt shift.

== Reference Fit <results.reference_fit>

For the selected model:

  - maximal cobalt cuts: $Delta E equiv |E_1-E_2| > 212$ keV and exclusion of 20 keV regions around $gamma$ peaks,
  - moderate $0 nu beta beta$ cuts: $E_"sum" > 1740 sqrt(2) = 2460.4$ keV,

the optimal cobalt shift correction is $Delta = 0.8$ keV/MeV. The following subsections describe the fit parameter values, spectral reconstruction, quality diagnostics, and component correlations for this reference fit.

=== Fit Values <results.fit_values>

The maximum likelihood fit yields #fit_ndbd.mode events in the $0 nu beta beta$ component. As shown in @tbl-fit-values, the fitted component strengths are consistent with background model expectations within uncertainties.

#figure(
  grid(
    columns: 2,
    column-gutter: 1em,
    align: horizon,
    btable(
      columns: 3,
      align: (left, right, right, right),

      [*Comp.*], [*Mode*], [*Expected ($mu plus.minus sigma$)*],

      [co60],   [#fit_co60.mode],   [$#fit_co60.exp  plus.minus #fit_co60.std$],
      [u238],   [#fit_u238.mode],   [$#fit_u238.exp  plus.minus #fit_u238.std$],
      [th232],  [#fit_th232.mode],  [$#fit_th232.exp plus.minus #fit_th232.std$],
      [others], [#fit_others.mode], [$#fit_others.exp plus.minus #fit_others.std$],
      [ndbd],   [#fit_ndbd.mode],   [#fit_ndbd.exp],
    ),
    image("../images/05_prior_post_mpl.png"),
  ),
  caption: [Left: fit results for the reference model. Mode (number of events) given in events; expected values and standard deviations are from the JAGS background model. Right: prior and posterior distributions for the background components.]
) <tbl-fit-values>

=== Spectral Reconstruction <results.spectral_reconstruction>

@fig-esum-ediff shows the fit in the rotated $E_"sum"$ and $E_"diff"$ axes, overlaid on the observed data. The residuals show no significant deviations.

#figure(
  grid(
    columns: 2,
    column-gutter: 1em,
    align: horizon,
    image("../images/05-01-esum_projection.png"),
    image("../images/05-01-ediff_projection.png"),
  ),
  caption: [Spectral reconstruction of the fit projected onto the $E_"sum"$ (left) and $E_"diff"$ (right) axes, with pull residuals.]
) <fig-esum-ediff>


=== Fit Quality and Convergence <results.fit_quality>

All $hat(R)$ convergence diagnostics are less than 1.001, indicating good Markov chain mixing. @fig-trace shows the MCMC trace for each component parameter across the four chains.

#figure(
  grid(
    columns: 2,
    column-gutter: 1em,
    align: horizon,
    image("../images/05-fit_trace.png"),
    image("../images/05-fit_energy.png"),
  ),
  caption: [Left: MCMC trace plots for all fit components; stable horizontal bands across chains confirm convergence. Right: fit energy spectrum.]
) <fig-trace>

=== Component Correlations <results.correlations>

@fig-correlations shows the pairwise posterior correlations between fit components. Some correlation is present between the $""^(232)$Th and "others" components, as expected from the similarity of their spectral shapes in this energy region. Crucially, no significant correlation is observed between the $0 nu beta beta$ and cobalt components, validating the shift correction procedure. #footnote[discuss 0v and other correlation, show some numbers for expected components 2nu, muon etc here or section above.]

#figure(
  image("../images/05-01-posterior_pair_kde.png"),
  caption: [Pairwise posterior correlation matrix for all fit components.]
) <fig-correlations>


== Bias correction

To correct the best-fit posterior for scaling and bias using toy MC studies, we evaluate the bias and scaling of each model. We run 200 toy fits on data generated from the best-fit (mode) spectrum with a known number of $0 nu beta beta$ events ($n_(0nu"-inj")$) and record the fitted number of events $n_(0nu"-fit")$.

- $n_(0nu"-inj") \~ text("Uniform")(5, 200)$
- fit components $f_i$ are the same as the original fit.
- Best-fit modes $c_i$ are used to generate the toy data, with the injected $0 nu beta beta$ events added on top: 
    
    $lambda_"tot" = sum_i c_i + n_(0nu"-inj")$.
// - Toy component strengths $n_i$ are the best fit values from the original posterior.

We plot the best-fit values and 68% confidence intervals (as errors) against the true $n_(0nu)$#footnote[discuss the lower cutoff for the fits]. Linear fits to these points yield the bias and scaling parameters. For each model (varying Cobalt shift), for all models scaling and bias is consistent with unity ($1-2 sigma$ variation).#footnote[add table]


#figure(
  image("../images/05-bias_focus_gamma10_v150_u1740.png"),
  caption: "Marginalization over Co shifts"
)

== Marginalizing over the Cobalt shift <results.marg_co60>

As discussed in Section 4.5, a Cobalt shift of $Delta = 0.8$ keV/MeV yields the highest likelihood. However, to utilize the full information, we manually marginalize over all fitted shifts by fitting a 2nd order polynomial to the maximum likelihood values of the each shift fit. From the polynomial fit, we estimate the best fit cobalt shift to be $Delta = 1.0$ keV/MeV, which we select as the final best fit.

We then apply the bias correction to each sample in the refence posterior distribution to obtain the best fit posterior on the $n_(0nu beta beta)$
#grid(
  columns: 2,
  column-gutter: 1em,
  align: horizon,
  image("../images/05-bias_marg.png"),
  image("../images/05_ref_cs10_raw_vs_corrected.png"),
  )





== Best fit 0vbb rate

After applying bias correction to the reference posterior ($"cs" = 1.0$), we find the bias-corrected $0 nu beta beta$ event count and convert it to a decay rate using the formula:

$ Gamma = n_(0 nu) / (f_(130) dot N_A dot (1000 / M_("TeO"_2)) dot "MT" dot epsilon_("MC") dot epsilon_("PSA")) $

The bias correction parameters for the reference fit (500 toy fits, $n_"inj" in [5, 200]$, excluding $n_"inj" < 15$) are $a = (1.010 plus.minus 0.016)$ and $b = -(2.46 plus.minus 1.55)$, consistent with unity slope and negligible offset.

#figure(
  grid(
    columns: 2,
    column-gutter: 1em,
    align: horizon,
    btable(
      columns: 2,
      align: (left, right),
      table.header([*Quantity*], [*Event Count*]),
      [Mode (measured)],     [$12.8$ ],
      [Mode (corrected)],    [$14.4$ ],
      [68% HDI],             [$[6.45, 26.44]$],
      [90% C.I.],            [$< 26.9$],
      table.header([*Quantity*], [*Value*]),
      [$f_(130)$],           [$0.34167$],
      [$"MT"$],            [$1038.4 "kg" dot "yr"$],
      [$epsilon_("MC")$],    [$0.03547$],
      [$epsilon_("PSA")$],   [$0.98$],
    ),
    image("../images/05_ref_cs10_raw_vs_corrected.png"),
  ),
  caption: [Left: bias-corrected $0 nu beta beta$ event count and rate conversion parameters for the reference fit ($"cs" = 1.0$, `gamma10_v150`, `u1740`). Right: raw (blue) vs bias-corrected (orange) posterior distributions.]
)

The resulting $0 nu beta beta$ decay rate and half-life limit are:

#figure(
  grid(
    columns: 2,
    column-gutter: 1em,
    align: horizon,
    btable(
      columns: 2,
      align: (left, right),
      table.header([*Quantity*], [*Value*]),
      [$Gamma$ (mode)],        [$3.3 times 10^(-25) "yr"^(-1)$],
      [68% HDI ($Gamma$)],     [$[1.5,\ 4.9] times 10^(-25) "yr"^(-1)$],
      [90% C.I. ($Gamma$)],    [$< 5.8 times 10^(-25) "yr"^(-1)$],
      [$T_(1/2)$ (90% C.L.)],  [$> 1.2 times 10^(24) "yr"$],
    ),
    image("../images/05-rate_final_marg.png"),
  ),
  caption: [Left: bias-corrected $0 nu beta beta$ rate and half-life results for the reference fit. Right: rate posterior distribution.]
)
