#import "common.typ": *
#import "@preview/akatable:0.1.0": academic-table

= Fit Model <fit_model>

The M2 Region of Interest (ROI) contains a limited number of spectral features. To accommodate the sparse event statistics, the 84 background components are consolidated into four primary source categories. To optimize convergence, the phase space is rotated, the background distributions are smoothed, and localized cuts and parametric shifts are introduced to mitigate discrepancies in spectra perticularly near the #Co60 sum peak lineshape.

== Merged Input Spectra <sec:merged_spectra>

The background model is simplified into four source-specific composite spectra: #Co60, #U238, #Th232, and all remaining ('other') contaminants. The prominence of the 1173 keV and 1334 keV gammas, alongside specific #U238 and #Th232 chain peaks, dictates this grouping.

To construct these composites, simulated events passing the selection criteria (@sec:experimental_data) are binned and normalized by their respective Background Model (BM) posterior modes#footnote[See @sec:bm_mode for mode algorithm]. The individual scaled histograms are then summed within their respective categories. @tbl-merged-comps details the leading contributors to each merged component. Both the expected integral counts and normalized shape templates are archived for the fit.

#align(center)[
#box(width: 90%, [
#academic-table(
  [Merged component event counts derived from JAGS-scaled Monte Carlo simulations. The three most prominent sources within each category are listed.],
  header: ([*Component (Top Sources)*], [*Sources*], [*Expected Events*], [*% Total*]),
  (
    [*#Co60*], [4], [1904], [62.3%],
    [  CuNOSV], [], [1865], [],
    [  TeO2], [], [20], [],
    [  S600mK], [], [19], [],
    
    // U238 component
    [*#U238*], [28], [342], [11.2%],
    [  CuNOSV Surf], [], [105], [],
    [  TopPb], [], [89], [],
    [  S10mK Surf], [], [53], [],
    
    // Th232 component
    [*#Th232*], [18], [641], [21.0%],
    [  RomanPb], [], [231], [],
    [  S600mK], [], [205], [],
    [  CuNOSV Surf], [], [114], [],
    
    // Other component
    [*Other*], [34], [127], [4.1%],
    [  Muons], [], [88], [],
    [  TeO2 #super[110m]Ag], [], [39], [],
  ),
  columns: (1.5fr, auto, auto, auto),
  label: <tbl-merged-comps>,
) 
])
]

#figure(
  image("../images/03-ucomps_stacked_0v.png", width: 100%),
  caption: "Background and signal components projected on the u-axis."
) 
#figure(
  image("../images/03-vcomps_stacked_0v.png", width: 100%),
  caption: "Background and signal components projected on the v-axis."
) 

== Rotation and Smoothing <sec:rotation_smoothing>

The diagonal orientation of the M2 ROI introduces binning artifacts at the boundaries by introducing partially included bins. To reduce this effect, the primary crystal energies ($E_1$, $E_2$) are mapped via a $45^(degree)$ rotation into orthogonal $u$ and $v$ coordinates:
$ u equiv (E_1+E_2)/sqrt(2), quad v equiv (E_1-E_2)/sqrt(2) $

To prevent sampling failures driven by zero-density bins, the merged background histograms undergo a single pass of a $5 times 5$ smoothing kernel ('k5b') via the ROOT `TH2D.Smooth()` algorithm prior to integration. The kernel weight matrix is detailed in @fig-rot-bak-kernel.



#let k5b = (
  (0, 1, 2, 1, 0),
  (1, 2, 4, 2, 1),
  (2, 4, 8, 4, 2),
  (1, 2, 4, 2, 1),
  (0, 1, 2, 1, 0),
)

#figure(
  grid(
    columns: 2,
    column-gutter: 2em,
    align: horizon,
    image("../images/03-rotated2D_expected_spectra.png", height: 160pt),
    table(
      columns: 5,
      align: center,
      inset: (x: 10pt, y: 7pt),
      stroke: 0.5pt,
      fill: (x, y) => {
        let v = k5b.at(y).at(x)
        color.mix((luma(245), 8 - v), (orange, v))
      },
      ..k5b.flatten().map(v => [#v]),
    ),
  ),
  caption: [
    Left: rotated background and signal components in the ROI. \
         Right: 5x5 smoothing kernel weights ('k5b').
  ]
)<fig-rot-bak-kernel>

== #ndbd Cutoff <sec:ndbd_cutoff>

Approximately 10% of dual-site #ndbd events do not deposit their full energy within the adjacent detectors, generating a low-energy tail ($<2510$ keV) that can bias the likelihood estimation if modeled imperfectly. To evaluate and negate this bias, progressive energy thresholds are tested to exclude the low-energy continuum:

1. No threshold (100% boundary efficiency).
2. $u < 1740 "keV"$: Excludes the generic flat continuum region (98.6% efficiency).
3. $u < 1775 "keV"$: Excludes the continuum up to the X-ray escape peak (95.8% efficiency).

// TODO:? Section 4.2 quantifies the bias reduction to justify the nominal threshold selection.

== Modeling #Co60 sum peak <sec:co60_modeling>

A systematic difference exists between the #Co60 sum peak (1173 keV and 1334 keV coincidences) reconstructed in data compared to simulations. Minor lineshape deviations in a volume sparse of other features heavily distort fit convergences. Specifically, observed #Co60 sum events reconstruct at marginally higher energies than theoretical predictions, an effect that becomes more pronounced as the energy asymmetry between two depositions increases (see @co-shift-check).

To correct for this, events within 20 keV of the individual $gamma$ lines ($|E_(1|2) - E_gamma| < 10$)#footnote[abs|E-E_co| > 10, validate] and events presenting nearly symmetric energy distriibution ($Delta E = |E_1 - E_2| < 212 "keV"$, or $v < 150$) are rejected. For remaining events, we model an energy-dependent shift applied to the primary energy deposit, $E_1$#footnote[Add analysis efficiencies as a table].

$ E_1 ' = E_1(1 + Delta), quad Delta in {0, 3, 5, 8, 10, 12, 15} times 10^(-3) $

The analysis manually marginalizes over this variable shift to extract the unaffected #ndbd rate (see @methods).

#figure(
  grid(
    rows: 2,
    column-gutter: 1em,
    align: left,  
    image("../images/co_shift_check_g10v150_vlo.png", width: 100%),
    image("../images/co_shift_check_g10v150_vmid.png", width: 100%)
  ),
  caption: "Co60 peak modeling",
)<co-shift-check>
