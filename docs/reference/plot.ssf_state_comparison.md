# Plot validation scores or changes along the trajectory

Plot validation scores or changes along the trajectory

## Usage

``` r
# S3 method for class 'ssf_state_comparison'
plot(
  x,
  type = c("scores", "track"),
  reference = NULL,
  from = NULL,
  to = NULL,
  threshold = 0.7,
  coords = c("x0", "y0", "x", "y"),
  segment = NULL,
  step_index = NULL,
  ...
)
```

## Arguments

- x:

  An
  [`compare_ssf_states()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_states.md)
  result.

- type:

  Show held-out `scores` or a `track` comparison.

- reference:

  Optional candidate ID. Scores become paired gains per step against
  that candidate on each identical block, including its family.

- from, to, threshold:

  Passed to
  [`compare_ssf_state_labels()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_state_labels.md)
  for a map.

- coords:

  Names of origin x/y and endpoint x/y columns.

- segment:

  Optional segment IDs to display.

- step_index:

  Optional observed row indices to display on a map. These refer to the
  complete trajectory before selecting `segment`.

- ...:

  Reserved for future use.

## Value

A ggplot object. Failed or missing fold scores are not plotted. Inspect
`block_scores` and `candidates` for coverage and failure details.
