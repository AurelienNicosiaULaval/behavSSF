# Show how inferred states acquire an ecological interpretation

Combines a labelled state table, a track, observed movement
distributions and coefficient intervals. The distributions group the
observations by their Viterbi labels; they are descriptive, not
independent validation. The same state colours are used throughout.
Coefficient intervals condition on the fitted model and labels and do
not adjust for model or taxonomy selection.

## Usage

``` r
plot_ssf_states(
  x,
  panel = c("overview", "table", "track", "distributions", "coefficients"),
  step_scale = 1,
  step_unit = "input units",
  segment = NULL,
  step_range = NULL,
  coords = c("x0", "y0", "x", "y"),
  terms = NULL,
  assessment = NULL,
  coord_unit = "input units",
  colours = NULL,
  term_labels = NULL,
  angle_bandwidth = 0.35
)
```

## Arguments

- x:

  A fitted `ssf_fit` or selected `guided_ssf` object.

- panel:

  `"overview"` combines the views. Use `"table"`, `"track"`,
  `"distributions"` or `"coefficients"` to draw them separately.

- step_scale:

  Positive multiplier applied to step lengths and their SDs. For input
  in kilometres, use 1000 to display metres.

- step_unit:

  Label for the resulting distance unit. Units are not inferred.

- segment, step_range, coords:

  Passed to
  [`plot.ssf_fit()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot.ssf_fit.md)
  for the map only. All other panels use the complete fitted
  observations.

- terms:

  Optional coefficient names for the coefficient panel. By default it
  shows angular, target, environmental and lag terms, excluding the
  distance-distribution parameters.

- assessment:

  Optional accepted assessment from
  [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md)
  for this fitted model. This makes a post-fit ecological interpretation
  visible without fitting the model again.

- coord_unit:

  Unit label for the original map coordinates. Coordinates are not
  rescaled by `step_scale`.

- colours:

  Optional colour vector, named by state label or in fitted state order.
  This can keep the same ecological colours across models.

- term_labels:

  Optional named character vector of readable coefficient names, indexed
  by the names in `terms`.

- angle_bandwidth:

  Circular smoothing bandwidth in radians, shared by all states. The
  default is 0.35 (about 20 degrees).

## Value

A ggplot object, or a patchwork for combined views. ggplot2 is required;
combined views also require patchwork. The turning-angle panel requires
ggcircular. Missing observations or covariance are identified in the
relevant panel.
