# Define a configurable state taxonomy

Labels are assigned in increasing order of a fitted state signature. A
rule can then check separate movement, selection, or transition
contrasts. The same numerical and occupancy checks as in the guided
workflow apply.

## Usage

``` r
ssf_taxonomy(
  name,
  labels,
  formula = ~step + log_step + cos_angle,
  order_by = "mean_step",
  rules = NULL
)
```

## Arguments

- name:

  Unique identifier for the candidate. The four built-in names are
  reserved.

- labels:

  One to six unique state labels, in increasing signature order.

- formula:

  One-sided SSF formula. The default uses the canonical gamma terms, or
  only `cos_angle` if the fit specifies a length family. Supply angular
  and selection terms for a custom step-family formula. Assessment
  requires these terms in the fitted model; the fitted model may contain
  additional terms.

- order_by:

  Name of a numeric movement metric from
  [`movement_metrics()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/movement_metrics.md)
  (`mean_step`, `sd_step`, `persistence`, `concentration`, `shape`, or
  `rate`), or a function taking the unlabelled metrics and coefficient
  matrix and returning one score per state. Tied or non-finite scores
  fail the label assignment.

- rules:

  Optional function of a single list containing `metrics` (rows named by
  label), `coefficients` and `standard_errors` (columns named by label),
  `transition` (rows and columns named by label), and `occupancy` (named
  by label). Its return value is a nonempty logical vector with unique
  names. Each FALSE or NA rejects the candidate. Standard errors are NA
  if the fit has no covariance matrix. The list also contains the joint
  natural-scale coefficient `covariance` and a `parameter_index` table
  identifying each row by `state` and `term`. Use
  [`ssf_contrast()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_contrast.md)
  to form a contrast without numeric parameter indices.

## Value

An `ssf_taxonomy` specification for
[`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md)
and
[`fit_guided_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_guided_ssf.md).

## Examples

``` r
ssf_taxonomy("step_regimes", c("Short steps", "Long steps"),
  rules = function(s) {
    c(length_ratio = s$metrics["Long steps", "mean_step"] >
        2 * s$metrics["Short steps", "mean_step"])
  })
#> $name
#> [1] "step_regimes"
#>
#> $labels
#> [1] "Short steps" "Long steps"
#>
#> $formula
#> ~step + log_step + cos_angle
#> <environment: 0x10601e390>
#>
#> $default_formula
#> [1] TRUE
#>
#> $order_by
#> [1] "mean_step"
#>
#> $rules
#> function (s)
#> {
#>     c(length_ratio = s$metrics["Long steps", "mean_step"] > 2 *
#>         s$metrics["Short steps", "mean_step"])
#> }
#> <environment: 0x1060177e8>
#>
#> attr(,"class")
#> [1] "ssf_taxonomy"
```
