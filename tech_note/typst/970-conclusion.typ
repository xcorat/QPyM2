#import "common.typ":*
#import "settings.typ": btable

= Conclusion <conclusion>

Combining the bias-corrected reference fit with the systematic marginalization described in @systematics, the final #ndbd decay rate and half-life from the mode of the posterior are //summarized in @fig-final-results.

#rect(
  stroke: black,
  // inset: 10pt,
  // radius: 4pt,
  fill: gray.lighten(90%),
  width: 100%
)[
  $ Gamma^(#ndbd) = 3.54  ""_(-1.54)^(+1.88) "stat" ""_(-0.85)^(+1.04) "syst" times 10^(-25)   thin thin "yr"^(-1) $
]

And we set a 90% confidence level Bayesian lower bounds for the #ndbd half life in the M2 analysis as:

#rect(
  stroke: black,
  // inset: 10pt,
  // radius: 4pt,
  fill: gray.lighten(90%),
  width: 100%
)[
  $ T^(#ndbd)_(1/2) > 1.09 times 10^(24) thin "yr" $
]

The median exclusion sensitivity (90% C.I) obtained by running 200 toy fits with background only sampled data is $1.77 times 10^(24) thin "yr"$. This result is consistant with the expected sensitivity even though given the higher bias correction for the specific Cobalt shift, we see a slight positive fluctuation of the signal. The probability of setting a weeker limit than one obtained is 6.5%.


#figure(
  grid(
    columns: (1fr, 1fr),
    // column-gutter: 1.5em,
    align: bottom,
    image("../images/06_syst_rate_new_final_rate_posterior.png"),
    image("../images/07_excl_sens.png"),
    

  ),
  caption: [Left: Systematic-marginalized #ndbd rate posterior. Red: mode; green dashed: 90% C.L.; shaded band: 68% HDI. Right: Final M2 #ndbd results.],
) <tbl-final-results>

#figure(
  grid(
    align: horizon,
    btable(
      columns: 2,
      align: (left, right),
      table.header([*Quantity*], [*Value*]),
      [$Gamma$ (mode)],
        [$3.54 times 10^(-25) thin "yr"^(-1)$],
      [$Gamma plus.minus$ stat],
        [$""^(+1.88)_(-1.54) times 10^(-25) thin "yr"^(-1)$],
      [$Gamma plus.minus$ syst],
        [$""^(+1.04)_(-0.85) times 10^(-25) thin "yr"^(-1)$],
      [$Gamma$ (90% C.L.)],
        [$< 6.36 times 10^(-25) thin "yr"^(-1)$],
      [$T^(0 nu)_(1\/2)$ (mode)],
        [$1.96 times 10^(24) thin "yr"$],
      [$T^(0 nu)_(1\/2) plus.minus$ stat],
        [$""^(+1.51)_(-0.68) times 10^(24) thin "yr"$],
      [$T^(0 nu)_(1\/2) plus.minus$ syst],
        [$""^(+0.62)_(-0.45) times 10^(24) thin "yr"$],
      [$T^(0 nu)_(1\/2)$ (90% C.L.)],
        [$> 1.09 times 10^(24) thin "yr"$],
    )
  ),
  caption: [Table of results]
) <tbl-final-results>

