#import "common.typ": *
#import "@preview/akatable:0.1.0": academic-table

= Fit Model <fit_model>

The M2 analysis fit region contains few background features. Therefore, we merge the individual BM Monte Carlo spectra into four main components representing these features. The fit is optimized by rotating the axes and selecting cuts that minimize discrepancies between simulated and observed spectra, particularly near peaks. After evaluating the fits and different models, we also introduce cuts on the fitted #ndbd region as well as cuts and a shift to the #Co60 component to better reconstruct data.

== Merged Input Spectra

This anlysis models the expected background as just four components. This simplification is fitting since we only see few observable features in the background spectra with our limited ROI and sparse data. Aside from the #Co60  peaks (1173 keV and 1334 keV), we primarily observe peaks from the #U238 and #Th232 decay chains. Consequently, we combine the 84 BM components into four source-specific groups: #Co60, #U238, #Th232, and 'other' (containing all remaining contaminants).

To generate the merged background spectra, we create histograms using the cuts and selection criteria from Section 2.1, scaling each by its best-fit value from the background model. These normalization factors are the modes of the BM fit posteriors #footnote[See Appendix for mode algorithm]. We then sum these scaled histograms to form the M2 background components. The sources are filtered by name to create the merged components, and the table below shows the major contributions for each component. We save both the integral (representing the expected event count) and the normalized spectrum for each component.

#academic-table(
  [Merged component event counts from JAGS-scaled MC, showing the three most prominent sources within each component.],
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
) 

#figure(
  image("../images/03-ucomps_stacked_0v.png", width: 100%),
  caption: "Background and signal components"
) 
#figure(
  image("../images/03-vcomps_stacked_0v.png", width: 100%),
  caption: "Background and signal components"
) 

== Rotation and Smoothing

Since the M2 ROI is a thin diagonal band, we rotate the fit axes by $45 degree$ to minimize binning artifacts at the boundaries. Additionally, we smooth the input spectra after merging to reduce empty bins that makes it harder for the fit to converge.

The rotation transforms the crystal energies, $E_1$ and $E_2$, into $u$ and $v$ axes:
$ u equiv (E_1+E_2)/sqrt(2), quad v equiv (E_1-E_2)/sqrt(2) $

Smoothing employs a $5 times 5$ kernel ('k5b') applied once to the merged histograms using the ROOT `TH2D.Smooth()` function.#footnote[Should we add the smoothed vs non plot comparison?] The kernel weights are shown in <fig-rot-bak-kernel>.



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
    "Left: rotated background and signal components in the ROI. \\\
         Right: smoothing kernel ('k5b')."
  ]
)<fig-rot-bak-kernel>

== #ndbd Cutoff

Not all dual-site #ndbd events deposit their full energy in the bolometers; approximately 10% record energies below 2510 keV, creating a low-energy tail. This tail can bias the fit if background features are imperfectly modeled. To mitigate this, we test ignoring the #ndbd tail below specific energy thresholds:

1. No cutoff
2. $u < 1740 "keV"$: excludes most of the flat region (98.6% efficiency).
3. $u < 1775 "keV"$: cuts near the X-ray escape peak (95.8% efficiency). #footnote[Validate]

Section 4.2 details the selection of the optimal cutoff.#footnote[Validate later]

== Modeling #Co60 sum peak

A key challenge in this analysis is the slight discrepancy between peak shapes in the background model simulations and data. Given the small ROI and few features, even minor lineshape deviations can bias the fit. We observe that the #Co60 sum peak (from 1334 keV and 1173 keV gammas) is reconstructed at slightly higher energies in data than in simulation. This shift is most pronounced when the energy distribution between crystals is asymmetric.#footnote[Update: what else did we check? 2615, DE, SE etc studies, did we check each part of Co too?]

To address this, we exclude regions near individual gamma peaks (minimal efficiency loss) and events with $Delta E = |E_1 - E_2| < 212 "keV"$ (v < 150). For remaining events, we model an energy-dependent shift applied to the primary energy deposit, $E_1$. We marginalize over these shifts manually to evaluate the best fit (see Section 4.2).  #footnote[Add analysis efficiencies as a table]

The shifted energy $E_1 '$ is defined as:

$ E_1 ' = E_1(1 + Delta), quad Delta in {0, 3, 5, 8, 10, 12, 15} times 10^(-3) $

#figure(
  grid(
    rows: 3,
    column-gutter: 1em,
    align: left,
  image("../images/co_shift_check_full.png", width: 100%),
  image("../images/co_shift_check_g10v150.png", width: 100%),
  image("../images/co_shift_check_g10v150_inv.png", width: 100%),
  ),

  caption: "Co60 peak modeling"
)

#figure(
  grid(
    rows: 2,
    column-gutter: 1em,
    align: left,  
    image("../images/co_shift_check_g10v150_vlo.png", width: 100%),
    image("../images/co_shift_check_g10v150_vmid.png", width: 100%)
  ),
  caption: "Co60 peak modeling"
)

===

==== TODO
- surface events? How much do we need to talk about the reasons for the shift?
- Select the plots above and also add comparison plots of ~2 shifts, none and 1.0?
- fit input generation appendix?
