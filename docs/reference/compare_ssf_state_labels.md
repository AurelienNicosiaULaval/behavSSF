# Compare step labels between two state counts

Each label in `to` is mapped to the label in `from` sharing the most
decoded steps. Several fine labels may share a parent. Ties use the
label order in the source taxonomy. This descriptive mapping is not a
biological validation or a measure of held-out classification accuracy.

## Usage

``` r
compare_ssf_state_labels(comparison, from, to, threshold = 0.7)
```

## Arguments

- comparison:

  An
  [`compare_ssf_states()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_states.md)
  result.

- from, to:

  Candidate identifiers from `comparison$candidates$candidate`. Both
  need numerically eligible full-data fits with assignable labels.

- threshold:

  Flag steps whose largest smoothed probability in `to` is below this
  display threshold. It is not a calibrated error probability.

## Value

A list with `steps`, a `crosswalk` table, the majority `mapping`, and
full-data taxonomy status. `.changed` compares the source label with the
mapped parent. `.ambiguous` describes the destination fit only.
