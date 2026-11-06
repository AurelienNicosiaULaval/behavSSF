# Compare step labels between two restart solutions

Compare step labels between two restart solutions

## Usage

``` r
compare_ssf_solutions(diagnostic, from = 1L, to = 2L)
```

## Arguments

- diagnostic:

  Output from
  [`diagnose_ssf_restarts()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/diagnose_ssf_restarts.md).

- from, to:

  Distinct solution numbers from `diagnostic$solutions`.

## Value

An `ssf_solution_comparison` with a `steps` data frame in fitted
observation order, a `transitions` count table, and `posterior_from` and
`posterior_to` matrices with columns named by taxonomy label. The step
table preserves observed columns, including coordinates and timestamps
when supplied. Added columns are `.step_index`, `.segment_step`,
`.state_from`, `.state_to` and `.changed`. `.source_row` identifies the
row in the input as prepared by the fitting workflow. Keep original fix
IDs in separate input columns.

## See also

[`diagnose_ssf_restarts()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/diagnose_ssf_restarts.md),
[`plot_ssf_changes()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot_ssf_changes.md)
