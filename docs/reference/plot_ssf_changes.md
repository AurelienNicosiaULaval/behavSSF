# Map steps whose state labels change between restart solutions

Map steps whose state labels change between restart solutions

## Usage

``` r
plot_ssf_changes(comparison, coords = c("x0", "y0", "x", "y"), segment = NULL)
```

## Arguments

- comparison:

  Output from
  [`compare_ssf_solutions()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_solutions.md).

- coords:

  Four column names for start x, start y, end x and end y, respectively.
  Use projected coordinates in common units.

- segment:

  Optional segment identifiers from the `ID` column.

## Value

A ggplot object. Every observed step is drawn as its own line. A
coloured endpoint identifies each relabelled step, including steps whose
line is too short to see at the map scale.

## See also

[`compare_ssf_solutions()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_solutions.md)
