#import "common.typ": *
#import "@preview/akatable:0.1.0": academic-table

= Analysis Methodology <methods>

This section formulates the Bayesian Markov Chain Monte Carlo (MCMC) framework applied to extract the #ndbd rate. It outlines the specific likelihood definitions, prior selections, validation checks against toy Monte Carlo datasets and the marginalization logic used to mitigate the #Co60 energy shift as well as methods used for bias correction and systematic studies.

== Bayesian Framework <sec:bayesian_framework>

The optimal #ndbd decay rate is extracted by comparing the observed spectrum against scaled background and signal templates constructed according to BM output (@sec:merged_spectra). This is achieved by maximizing an extended unbinned likelihood function over the component normalization factors with posterior obtained using the MCMC method.

Using Bayes' theorem, the posterior distribution of the fit parameters $arrow(theta)$ (i.e the source strengths for background and the signal components) is calculated as:

$ p(arrow(theta)|cal(D)) = ( p(cal(D)|arrow(theta)) p(arrow(theta)) ) / p(cal(D)) prop p(cal(D)|arrow(theta)) p(arrow(theta)) $

Here, $p(cal(D)|arrow(theta))$ represents the likelihood function for the parameters $arrow(theta)$ with respect to observed data $cal(D)$ and $p(arrow(theta))$ defines the prior knowledge. The best-fit values we cite correspond to the mode of the resulting posterior distribution#footnote[Mode, as opposed to mean or median is selected as the best estimator for this analysis]. We approximate this distribution utilizing the No-U-Turn Sampler (NUTS, an extension to the Hamiltonian Monte-Carlo sampler) implemented via the PyMC (v5) probabilistic programming framework. 

== Probability Model <sec:probability_model>

In our fit construction, the priors $p(arrow(theta))$ reflect the expected event counts from the BM fit and its systematic variations#footnote[See Appendix for details on the BM fit and systematic variations: TODO]. The likelihood $p(cal(D)|arrow(theta))$ evaluates the unbinned events within the ROI against the binned positional probabilities defined by the merged component configurations:

$
 p(cal(D)|arrow(theta)) = ( lambda^n e^(-lambda) ) / ( n! ) product_(i=1)^n f_i (arrow(E)_i|arrow(theta))
prop  e^(-lambda ) product_i ( sum_(k=0)^(m) lambda_k f_(i k) )
$

where $lambda = sum_(k=0)^(m) lambda_k$ represents the total expected event count across $m$ background components plus the signal. The terms $f_(i k)$ denote the normalized probability density value for component $k$ at the specific energy coordinates of data point $i$. $f_(i k)$ is calculated once and used for likelihood calculations while sampling, which reduces the computation heavy calls to bin-search algorithm.

== Prior Distributions <sec:prior_distributions>

Background component priors are defined as truncated normal distributions bounded below at zero#footnote[Implemented as `TruncatedNormal`; assess if additional details are required], while the signal prior acts as a strictly positive uniform distribution also bounded above at the total number of observed events. The mean ($mu$) and width ($sigma$) of each background prior are determined by parametrically propagating uncertainties from the BM fits. The generative procedure operates as follows:

1. Select a random vector from the BM posterior.
2. Sample from the maximal-correlation matrix#footnote[Check Appendix xxx for discussion of what this correlation matrix is] to integrate systematic variations.
3. Aggregate the relevant sub-components to predict the expected M2 event rate per sum category.
4. Iterate 100 times to construct an ensemble distribution, record $mu$ and $sigma$,
5. Assign $mu_k = mu$ and $sigma_k = 2 sigma$ for truncated Gaussian priors per component $k$.

== Fit Validation <sec:fit_validation>

Fit quality is evaluated by comparing empirical fit results against that of toy Monte Carlo datasets generated from the best-fit model. The posterior mode of the data fit is benchmarked against the distribution of modes from the toy ensemble for both the log-likelihood and the signal rate. To generate toy data ensambles, we take the fitted mode for each component (backgrounds and signal,) and sample a Poisson smeared number of events from the relevant histogram to create a list of toy events. Across all tested models, nominal fits match well with toy expectations; routinely, the best log-likelihood estimate for actual data falls within the 68% highest density interval (HDI) of the toy ensemble, and all empirical fits are within the expected distribution of signal rates derived from the toy datasets.


== Sensitivity Comparison <sec:sensitivity_comparison>

Model sensitivity is defined as the 68% C.L. upper limit on the signal rate, calculated using ensembles of background-only toy Monte Carlo datasets. The toy generation procedure is identical to the above only without any signal injection. Evaluation across all hyper-parameter configurations produces comparable sensitivity limits#footnote[Sensitivity was initially evaluated prior to unblinding in May 2024 to compare the models, whose plots are not reproduced here.].

== Fit Models <sec:fit_models>

The analysis tests three main categorical variations to the fit configuration:

+ #ndbd Cutoff: Cut the fitted #ndbd signal on total energy.
  - `none`: None,
  - `u1740`: $u<1740$ keV
  - `u1775`: $u<1775$ keV

+ #Co60 Removal: Remove regions from the fit (remove the events from both data and simulations)
  - `none`: None
  - `v_150`: $|v| < 150$ keV cut,
  - `gamma10`: $|E_(1|2) - E_gamma| < 10$ keV, $E_gamma in (1173, 1334)$
  - `gamma10_v150`: both `v150` and `gamma10`

+ #Co60 Peak Shift $Delta$: Scale the $E_1$ energy in Cobalt component, $E_1^' = (1+Delta)E_1$
  - $Delta in {0.0, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5}$

The fits from each point in the grid of above paraters were tested for fit validity and sensitivity (see <sec:loglike-validation>). During blinded analysis, we decided to keep the more loose cut on #ndbd region (`u1740`) and strongest cut on the #Co60 peaks (`gamma10_v150`). However, looking at unblinded either of these measures can be used to rank models authoritatively, and we resort to a more heuristic selection procedure as discussed below.

#figure(
  image("../images/04-sens_comparison.png", width: 100%),
  caption: "Sensitivity (90% C.L. upper limit on signal rate) comparison"
) 
// #align(center)[
// #box(width: 70%, [
// #academic-table(
//   [Hyper-parameters tested across fit models.],
//   header: ([*Parameter*], [*Values*]),
//   (
//     [$0 nu beta beta$ Cutoff],  [None; $u < 1740$ keV; $u < 1770$ keV],
//     [$"Co"^60$ Removal],        [None; $|v| < 150$ keV; Co $gamma$ peaks; both],
//     [Co Peak Shift $Delta$],    [$\{0, 3, 5, 8, 10, 12, 15\} times 10^(-3)$],
//   ),
//   columns: (auto, 1fr),
//   label: <tbl-hyperparams>,
// )
// ])
// ]

=== Selecting optimal cobalt shifts <sec:model_selection_shift>

To validate the cobalt shift procedure without signal bias, an initial quasi-blinding protocol excluded the $0 nu beta beta$ peak region (2515--2540 keV, maintaining $tilde$95% $0 nu beta beta$ coverage). This confirmed that the background normalization is independent of the cobalt shift parameters#footnote[Refer to pre-unblinding analysis (May 2024) for detailed pseudo-blinding metrics]. Following this verification, the analysis proceeds unblinded. The methodology for identifying the optimal cobalt shift via marginalization is detailed in @results.marg_co60.

=== Selecting Cobalt and signal cuts <sec:selecting_cuts>

Because the tested models have similar sensitivities and toy-fits agree with in each other, the reference model can be selected heuristically to maximize robustness. To minimize the influence of the cobalt shift, the maximal removal criteria are applied, excluding events within 20 keV of the cobalt gammas or sharing nearly symmetric energies ($|v| < 150$ keV). For the $2nu beta beta$ continuum cutoff, $u < 1740$ keV is chosen to effectively discard the flat continuum while preserving the regions most adjacent to the ROI. 

// == Systematics studies and marginalization

// - most systematics cannot be included in MCMC
// - perform descreet marginalization, 
//   - how to combine posteriors with many samples, 
//   - how to combine when +-1 sigma?
//   - why use HDI? what estimator to use given multiple posteriors?
// - how to combine all systematics studies,
//   - how to get the combined-marginalized posterior, left-right split-gaus

== Bias correction <sec:methods_bias_correction>

The unbinned fits are are evaluated and corrected for biases using toy fits with injected signal events. We run toy fits on data sampled using the best-fit (mode) spectrum and plot the fitted vs injected number of signal events.

+ Inject a uniform distribution of signal events: $n_(0nu"-inj") tilde text("Uniform")(5, 200)$.
+ Load the fit components $f_i$ equal to the original fit.
+ Run fit and record the posterior modes of #ndbd ($n_(0nu"-fitted")$), and standard deviations $sigma_(0nu)$ for each toy.
+ Fit a 1st degree polynomial ($n_(0nu"-fitted") = a n_(0nu"-inj") + b$) to the fitted vs injected ( $n_(0nu"-inj")$) number of signal events, with error bars representing the recorded $sigma_(0nu)$.

We record these parameters and their uncertainties to correct for bias in our final posteriors.

== Systematic treatment framework <sec:systematic_framework>

Incorporating all systematic effects directly into the MCMC sampler would require regenerating the background spectra for every point in each systematic-parameter space, which is computationally prohibitive. Instead, we run a small number of targeted refits around the reference model and compare the resulting posteriors to it.

=== Measuring individual systematics <sec:measuring_systematics>

For each perturbed refit, we compare the 68% highest-density interval (HDI) of the $n_(#ndbd)$ posterior to the reference. The HDI is the shortest interval containing 68% of the probability mass. This is analogous to a $1 sigma$ confidence interval, but better captures both displacement and asymmetric widening of the posterior in a way that a simple mode shift does not.

We record the outward shift of each edge separately:
- *Left shift* $Delta^-$: how far the lower HDI edge moves to the left relative to the reference.
- *Right shift* $Delta^+$: how far the upper HDI edge moves to the right.

Inward shifts (i.e. an edge moving toward the mode) are set to zero, since those indicate a narrowing of the posterior rather than a broadening.

When the systematic is sampled at two points (e.g. $plus.minus 1 sigma$ for the lineshape corrections), the maximum outward shift on each side across all samples and datasets is recorded. Contributions from systematic effects that are considered independent (Ex: per dataset posteriors for $Delta Q$, or completely different parameters,) are summed in quadrature:

$ sigma^+ = sqrt(sum_i (Delta^+_i)^2), quad sigma^- = sqrt(sum_i (Delta^-_i)^2) $

When many refits are available (e.g. the 500 background-model Monte Carlo samples), a single mixture posterior is generated by weighted addition where the weights $w_i$ are given by:

$ w_i = exp("LL"_i - "LL"_"max") / (sum_j exp("LL"_j - "LL"_"max") ) $

=== Producing the marginalized posterior <sec:marginalized_posterior>

Once the left and right systematic widths $sigma^(plus.minus)_k$ from all sources $k$ are obtained, they are added in quadrature independently on each side:

$ sigma^+_"tot" = sqrt(sum_k (sigma^+_k)^2), quad sigma^-_"tot" = sqrt(sum_k (sigma^-_k)^2) $

The reference posterior is then convolved with a split-Gaussian kernel of widths $sigma^(plus.minus)_"tot"$, producing the final systematic-marginalized posterior. The mode, 68% HDI, and 90% upper limit quoted in the results all derive from this final posterior. The numerical values of all individual contributions and posterior comparison plots are discussed in @systematics. //sec:systematics.
