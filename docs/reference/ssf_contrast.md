# Compute a coefficient contrast using labelled joint covariance

Contrasts condition on the fitted model and the selected labels. Wald
intervals do not adjust for selection of a model, solution or taxonomy.

## Usage

``` r
ssf_contrast(x, term, from, to = NULL, level = 0.95)
```

## Arguments

- x:

  A rule's input list, or the `context` returned by
  [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md).

- term:

  One or more distinct coefficient row names.

- from:

  Source state label. With `to`, the contrast is `to - from`.

- to:

  Destination label, or NULL to report the coefficient in `from`.

- level:

  Confidence level strictly between zero and one.

## Value

A data frame with `term`, `from`, `to`, `estimate`, `std_error`, `lower`
and `upper`. Uncertainty is NA if covariance is unavailable.
