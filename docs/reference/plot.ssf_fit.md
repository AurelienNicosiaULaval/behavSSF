# Plot fitted behavioural states

Draw a track coloured by its Viterbi labels or show the smoothed state
probabilities in step order. Guided fits use the ecological labels of
the selected taxonomy. The probability curves are conditional on the
fitted model and stop at segment boundaries.

## Usage

``` r
# S3 method for class 'ssf_fit'
plot(
  x,
  ...,
  type = c("auto", "track", "posterior"),
  segment = NULL,
  step_range = NULL,
  coords = c("x0", "y0", "x", "y")
)

# S3 method for class 'guided_ssf'
plot(
  x,
  ...,
  type = c("auto", "track", "posterior"),
  segment = NULL,
  step_range = NULL,
  coords = c("x0", "y0", "x", "y")
)
```

## Arguments

- x:

  An object from
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md)
  or
  [`fit_guided_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_guided_ssf.md).

- ...:

  Reserved for the generic.

- type:

  Use `"track"` for a map, `"posterior"` for state probabilities, or
  `"auto"` for a map when the chosen coordinate columns are available
  and a probability plot otherwise.

- segment:

  Optional segment identifiers from the fitted `ID` column.

- step_range:

  Optional pair giving the first and last step to draw, counted in
  fitted order after selecting segments.

- coords:

  Four numeric columns in the fitted observations for start x, start y,
  end x and end y, respectively. Use projected coordinates in common
  units.

## Value

A ggplot object.
