# Read the movement signature of each inferred state

Counts describe decoded steps. Length means and SDs describe the
intrinsic movement kernel, holding environmental terms constant. They
are baseline summaries when lengths have hourly or lag effects. Stay
probabilities are model estimates. A geometric mean run length is
returned only for homogeneous multi-state transitions; it is not an
observed bout length.

## Usage

``` r
ssf_state_summary(
  x,
  step_scale = 1,
  step_unit = "input units",
  assessment = NULL
)
```

## Arguments

- x:

  A fitted `ssf_fit` or selected `guided_ssf` object.

- step_scale:

  Positive multiplier applied to step lengths and their SDs. For input
  in kilometres, use 1000 to display metres.

- step_unit:

  Label for the resulting distance unit. Units are not inferred.

- assessment:

  Optional accepted assessment from
  [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md)
  for this fitted model. This makes a post-fit ecological interpretation
  visible without fitting the model again.

## Value

A data frame with state labels, decoded counts and fractions, fitted
movement summaries, stay-probability ranges and expected run lengths in
steps.
