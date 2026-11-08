# Compare model effects between inferred states

Returns differences between state-specific coefficients using their
joint covariance. This applies to any estimable numeric, factor or
interaction term in the fitted model. The intervals condition on that
model and its assigned labels; they do not adjust for choosing states,
rules or contrasts.

## Usage

``` r
ssf_state_contrasts(
  x,
  terms = NULL,
  comparisons = NULL,
  assessment = NULL,
  level = 0.95
)
```

## Arguments

- x:

  A fitted `ssf_fit` or selected `guided_ssf` object.

- terms:

  Coefficient names. The default includes all non-radial terms.

- comparisons:

  Optional data frame with `from` and `to` state labels. Each contrast
  is `to - from`. The default gives every unordered pair once.

- assessment:

  Optional accepted assessment from
  [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md)
  for this fitted model. This makes a post-fit ecological interpretation
  visible without fitting the model again.

- level:

  Confidence level, strictly between zero and one.

## Value

A data frame with coefficient names, state pairs, differences, standard
errors and confidence limits. Missing covariance gives NA limits.
