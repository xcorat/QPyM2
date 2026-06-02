#import "common.typ": *

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
#let fit_ndbd   = (mode: "11.7",  exp: "--",    std: "--",     p90: "24.6")

= Fit Results <results>

The final posterior is constructed by bias-correcting the selected model fit and converting the observed number of #ndbd events into a decay rate. This section details the reference fit, the bias correction procedure, and the selection of the best fit Cobalt shift.

== Reference Fit <results.reference_fit>

As detailed in @sec:fit_models The selected reference model utilizes:
- maximal #Co60 removal criteria ($|E_1-E_2| > 212$ keV and the exclusion of 20 keV regions surrounding the gamma peaks),
- $0 nu beta beta$ continuum threshold of $u < 1740$ keV ($E_"sum" > 2460.4$ keV). 

//The best fit cobalt shift parameter for this configuration is $Delta = 0.8 times 10^(-3)$, but the marginalization as described in <sec:marg_co60> yields $Delta = 1.0 times 10^(-3)$.

The maximum likelihood fit extracts #fit_ndbd.mode events for the $0 nu beta beta$ component. As detailed in @tbl-fit-values, the fitted background normalizations are close to expectations, even though we see a possible movement of events from the #Th232 to the 'others' component in the final numbers. Spectral reconstructions projected onto the rotated energy axes (@fig-esum-ediff) are well in agreement with the observed data, demonstrating structureless flat pull residuals.

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
  caption: [Left: Fit results for the reference model. Mode given in events; expected priors derived from Monte Carlo variants. Right: Prior densities and scaled posterior distributions for all background components.]
) <tbl-fit-values>

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

All $hat(R)$ convergence diagnostics #footnote[Describe what these are here or in an appendix] are less than 1.001, indicating good Markov chain mixing. @fig-trace shows the MCMC trace for each component parameter across the chains, confirming stability. 

@fig-correlations shows the pairwise posterior correlations between the fit components. Some correlation is present between the $""^(232)$Th and "others" components, as expected from the similarity of their spectral shapes. Crucially, no significant correlation is observed between the $0 nu beta beta$ signal and the cobalt components, validating that the #Co60 shift correction does not bias the signal rate#footnote[Discuss additional component correlations involving the continuum or muons here].

#figure(
  grid(
    columns: 2,
    column-gutter: 1em,
    align: horizon,
    image("../images/05-fit_trace.png"),
    image("../images/05-fit_energy.png"),
  ),
  caption: [Left: MCMC trace plots for all fit components; stable horizontal bands verify sampler convergence. Right: Fit energy spectrum.]
) <fig-trace>

#figure(
  image("../images/05-01-posterior_pair_kde.png", width: 90%),
  caption: [Pairwise posterior correlation matrix for the tested components in the reference model.]
) <fig-correlations>


== Bias Correction <results.bias_correction>

We run 200 toy fits with injected signal events and evaluate the bias and scaling of our fitted #ndbd posterior as detailed in @sec:methods_bias_correction. Across all tested models with varying cobalt shifts, the bias and scaling strictly conforms to unity within a $1-2 sigma$ statistical variation for all models. The bias correctd posterior, where we apply the correction $n_"corrected" = ( n_("sampled") - b ) / a $, yields $14.88$ events.

#figure(
  image("../images/05-bias_focus_gamma10_v150_u1740.png", width: 90%),
  caption: [Signal extraction bias and scaling calibration for the reference configuration.]
)

== Marginalizing over the Cobalt Shift <results.marg_co60>

As discussed in @sec:model_selection_shift, the cobalt shift of $Delta = 0.8$ keV/MeV yielded the highest likelihood. However, we manually marginalize over all the shifts to leverage full information from all 7 fits by fitting a second-order polynomial to the log-likelihood maxima vs the shifts.

This polynomial fit gives a continuous global best-fit estimate of $Delta = 1.0$ keV/MeV, which we select as our nominal reference fit. The full marginalized raw posterior due to cobalt shift is constructed by combining all these posteriors weighted by their maximum likelihood, which we consider as a systematic, discussed in @sec:systematics_co60.

#figure(
  grid(
    columns: 2,
    column-gutter: 1em,
    align: horizon,
    image("../images/05-bias_marg.png", width: 90%),
    image("../images/05_ref_cs10_raw_vs_corrected.png", width: 90%),
  ),
  caption: [Left: Polynomial marginalization over discrete cobalt shift fits. Right: Reconstructed raw vs bias-corrected posterior distributions for the optimal configuration.]
)

== Extracted 0vbb Rate <sec:extracted_rate>

Following the application of the bias correction to the optimal reference ($Delta = 1.0 times 10^(-3)$) model, the finalized absolute $0 nu beta beta$ event count is utilized to compute the decay rate via the expression:

$ Gamma = n_(0 nu) / (f_(130) dot N_A dot (1000 / M_("TeO"_2)) dot "MT" dot epsilon_("MC") dot epsilon_("PSA")) $

The bias formulation determined for the reference model configuration utilizes 500 toy simulations specifically within the injected scope of $n_"inj" in [20, 200]$. We remove the lower left quarter of the toy resutls ($n_"inj" < 20$) from the fit to remove possible artefacts due to the zero bound priors. This establishes the intercept at $b = -(3.46 plus.minus 2.55)$ and the gradient stringently at $a = (1.010 plus.minus 0.016)$. These represent near unity parameter scaling and a completely statistically negligible offset.

#figure(
  grid(
    columns: 2,
    column-gutter: 3em,
    align: horizon,
    btable(
      columns: 2,
      align: (left, right),
      table.header([*Quantity*], [*Event Count*]),
      [Mode (measured)],     [#fit_ndbd.mode ],
      [Mode (corrected)],    [$14.4$ ],
      [68% HDI],             [$[6.45, 26.44]$],
      [90% C.I.],            [$< 26.9$]
    ),
    btable(
      columns: 2,
      align: (left, right),
      table.header([*Quantity*], [*Value*]),
      [$f_(130)$],           [$0.34167$],
      [$"MT"$],            [$1038.4 "kg" dot "yr"$],
      [$epsilon_("MC")$],    [$0.03547$],
      [$epsilon_("PSA")$],   [$0.98$],
    ) 
  ),
  caption: [Left: Corrected $0 nu beta beta$ event count \ Right: Rate conversion parameters for the reference fit.]
)

The resulting $0 nu beta beta$ decay rate and half-life limit is:

$ T^(0 nu)_(1/2) < 1.2 times 10^(-24) "yr"^(-1) quad 90%"C.L" $

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
  // caption: [Left: bias-corrected $0 nu beta beta$ rate and half-life results for the reference fit. Right: rate posterior distribution.]
)
