# Assess a fitted ecological taxonomy

Tests strict mean-step ordering and travelling concentration greater
than foraging (or encamped). For HDT, the unadjusted 95 percent Wald
interval of the target coefficient must be entirely positive. A
significantly negative target coefficient does not establish attraction.
Numerical failure and failed ecological rules remain distinct in the
returned checks.

## Usage

``` r
assess_taxonomy(
  fit,
  taxonomy,
  convention,
  target = "cos_diff_theta_to_target",
  min_occupancy = 0.05,
  persistence_bounds = c(0.05, 0.95),
  min_mean_gap = 0,
  score_tolerance = 0.01
)
```

## Arguments

- fit:

  An `ssf_fit` object.

- taxonomy:

  One of `HDT_F_R`, `T_F_R`, `T_E`, or `Other`, or a specification
  created with
  [`ssf_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_taxonomy.md).

- convention:

  Kernel convention passed to
  [`movement_metrics()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/movement_metrics.md).

- target:

  Name of a target-direction coefficient for HDT.

- min_occupancy:

  Minimum Viterbi occupancy, between zero and one.

- persistence_bounds:

  Open lower and upper self-transition bounds. For daily transitions,
  the bounds apply at each of the 24 hours.

- min_mean_gap:

  Nonnegative minimum absolute gap between ordered means, in step-length
  units. This is not a significance test.

- score_tolerance:

  Maximum absolute objective gradient for accepting a numerically
  stationary fit. Missing or nonfinite scores fail this check.

## Value

A list with `accepted`, named logical `checks`, state `labels`,
`metrics`, `occupancy`, and `reasons`. Wald intervals do not adjust for
model selection and numerical convergence alone does not establish
adequacy.
