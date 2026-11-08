# Fit and assess the descending taxonomies of the reference workflow

By default, tries HDT/F/R (if a target is supplied), then T/F/R, T/E,
and Other. Supply an ordered list to use custom taxonomies or mix them
with built-in candidates. The first candidate passing its checks is
selected. Returns every attempted candidate, including numerical
failures. If no candidate passes, status is `no_accepted_model`; it does
not reinterpret numerical failure as evidence of homogeneous behaviour.

## Usage

``` r
fit_guided_ssf(
  data,
  convention,
  target = NULL,
  n_starts = 5L,
  short_iter = 30L,
  long_iter = 300L,
  spread = 0.5,
  seed = NULL,
  min_occupancy = 0.05,
  persistence_bounds = c(0.05, 0.95),
  min_mean_gap = 0,
  engine = c("native", "R"),
  score_tolerance = 0.01,
  normalization = c("importance", "exact"),
  step_ar = FALSE,
  taxonomies = NULL,
  step_family = NULL
)
```

## Arguments

- data:

  Case-control step data.

- convention:

  Required kernel convention, see
  [`movement_metrics()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/movement_metrics.md).

- target:

  Optional name of a target-direction column for the built-in HDT
  candidate. Custom candidates supply covariates in their formula.

- n_starts, short_iter, long_iter, spread, seed:

  Fitting controls, see
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).

- min_occupancy, persistence_bounds, min_mean_gap, score_tolerance:

  Assessment controls, see
  [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md).

- engine:

  Fitting engine passed to
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).

- normalization:

  importance for case-control data, exact for observed ssf_steps and
  compatible analytically normalised kernels.

- step_ar:

  Logical. Include state-specific length autoregression in every
  candidate when using exact normalisation.

- taxonomies:

  NULL for the reference sequence, or a nonempty ordered list of
  [`ssf_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_taxonomy.md)
  specifications and/or built-in taxonomy names. Include `"Other"`
  explicitly to try the one-state fallback in a custom sequence.

- step_family:

  Optional length family passed to
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).

## Value

A `guided_ssf` object with selection status, selected fit and
assessment, and the complete candidate history. No model averaging,
post-selection confidence calibration or HSMM is performed.
