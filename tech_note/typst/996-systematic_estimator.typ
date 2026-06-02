
== LSS Systematic: Alternative Methods (M2, M3) <sec:lss_alt_methods>

We evaluated two additional estimators for the LSS systematic uncertainty.

*Method M2 — 68% HDI edge shift*: the systematic for each dataset is taken as the shift of each HDI(68%) edge relative to the nominal:

$ Delta_"lo"^((i)) = "HDI_lo"_"syst"^((i)) - "HDI_lo"_"nom", quad Delta_"hi"^((i)) = "HDI_hi"_"syst"^((i)) - "HDI_hi"_"nom". $

#figure(
  image("../images/06-lss_m2_bars_cs1p0.png", width: 90%),
  caption: [LSS systematic per dataset — M2 (68% HDI edge shift, $"cs" = 1.0$).],
)

#table(
  columns: (auto, auto, auto),
  align: (left, right, right),
  table.header([*Source*], [*hi $n_+$*], [*lo $n_-$*]),
  [$Delta Q$],    [+4.508], [$-$2.732],
  [$sigma$-scl],  [+1.833], [$-$0.827],
)

*Method M3 — half-interval from mode*: the left (right) arm is the 31.73th (68.27th) percentile of posterior samples below (above) the mode.

#figure(
  image("../images/06-lss_m3_bars_cs1p0.png", width: 90%),
  caption: [LSS systematic per dataset — M3 (half-interval from mode, $"cs" = 1.0$).],
)

#table(
  columns: (auto, auto, auto),
  align: (left, right, right),
  table.header([*Source*], [*hi $n_+$*], [*lo $n_-$*]),
  [$Delta Q$],    [+4.503], [$-$2.671],
  [$sigma$-scl],  [+2.733], [$-$0.488],
)

M1 yields the largest (most conservative) uncertainties and is used as the primary result.
