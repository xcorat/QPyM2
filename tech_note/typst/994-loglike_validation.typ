#import "common.typ": *

== Log-likelihood Validation <sec:loglike_validation>

== Toy fits

=== Sampling algorithm:

```
sample_nevs = np.random.poisson(hist.sum())
samples = np.random.choice(np.arange(len(flat_hist)), size=sample_nevs, p=flat_hist/hist.sum())
```

#page(flipped: true)[
  #align(center + horizon)[
    #figure(
      image("../images/04-loglike_validation_grid.png", width: 90%),
      caption: "Log-likelihood validation and compatibility across tested model variations."
    )
  ]
]
