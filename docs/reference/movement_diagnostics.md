# Compute observable movement summaries without crossing segment boundaries

The function does not compare decoded states with latent simulated
states. Each segment includes its supplied starting point for
displacement and MSD. Short-step runs are measured within the available
segment, including runs censored by segment boundaries, identically for
observed and simulated data.

## Usage

``` r
movement_diagnostics(
  track,
  segments,
  lags = c(1L, 2L, 4L, 8L, 16L),
  short_step = 0.05
)
```

## Arguments

- track:

  Data frame with ID, x, y, step and angle columns, in order.

- segments:

  Starting metadata as supplied to
  [`simulate_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/simulate_ssf.md).

- lags:

  Positive integer displacement lags in numbers of steps.

- short_step:

  Positive step threshold in the same units as coordinates.

## Value

Named numeric vector. Distances use input units, MSD squared units, and
short-run durations numbers of steps. Undefined correlations are NA.
