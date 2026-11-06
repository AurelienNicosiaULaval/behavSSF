# Evaluate importance-sampling or exact SSF emission densities

For J controls and the case in each choice set, returns eta_case -
log(sum(exp(eta) / w)) + log(J). The case is included in the denominator
and its numerator is unweighted, matching the published code. This is an
approximate density, not a probability over the finite choice set. The
measure underlying w must agree with the kernel being fitted. With exact
normalisation, controls and weights are absent and the returned
emissions are properly normalised canonical-kernel densities.

## Usage

``` r
ssf_log_emissions(
  data,
  formula,
  coefficients,
  normalization = c("importance", "exact")
)
```

## Arguments

- data:

  Case-control
  [`prepare_ssf_data()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/prepare_ssf_data.md)
  data, or observed
  [`prepare_ssf_steps()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/prepare_ssf_steps.md)
  data for exact normalisation.

- formula:

  One-sided SSF formula without an identifiable intercept.

- coefficients:

  Numeric matrix, with design columns in rows and states in columns. Row
  names must exactly match the design column names.

- normalization:

  importance uses the published approximation; exact uses observed
  ssf_steps and the analytic canonical-kernel normaliser.

## Value

A matrix of log emission densities, one row per stratum.
