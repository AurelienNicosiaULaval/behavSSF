# Recover movement summaries under an explicit kernel convention

For a kernel proportional to exp(beta_step \* r + beta_log \* log(r)),
the gamma shape is beta_log + 2 for a density with respect to planar
area, and beta_log + 1 for a density with respect to step length and
angle. The rate is -beta_step. On a finite disk this field can be
negative; the moments then describe the bounded radial kernel, not a
full-plane gamma distribution. These are intrinsic, habitat-free kernel
summaries; they are not the marginal distribution under spatially
varying habitat. With `lag_log_step`, mean and standard deviation
describe the baseline gamma kernel at the first step of a segment. The
mean also anchors the autoregression for subsequent steps. Hourly
step-length effects are separate from these baseline summaries. The
reported concentration is the intrinsic persistence proxy abs(beta_cos).
With a target-direction effect, use
[`angular_parameters()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/angular_parameters.md)
for the full conditional angular concentration.

## Usage

``` r
movement_metrics(coefficients, convention, step_family = NULL, step_max = Inf)
```

## Arguments

- coefficients:

  A numeric coefficient matrix with rows `step`, `log_step`, and
  `cos_angle` and states in columns.

- convention:

  Required: `"area"` for the article's +2 convention or `"step_angle"`
  for a step-length/angle density. The caller must verify consistency
  with the measure of the sampling density w.

- step_family:

  Optional length family supplied to
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).
  Use the named distribution parameters from the fitted coefficient
  matrix.

- step_max:

  Upper radial support, or Inf. With finite support the canonical mean
  and SD are conditional on (0, step_max\], and a positive linear step
  coefficient can be valid. Set this in
  [`prepare_ssf_data()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/prepare_ssf_data.md).

## Value

A data frame with state, shape/rate where applicable, mean/sd step
length, signed persistence, concentration, and kernel validity. Invalid
gamma kernels have NA moments, never fabricated finite values.
