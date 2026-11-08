# Validate and organise case-control step data

Each ID must identify an independent, uninterrupted segment. Supply
separate IDs across observation gaps, even for the same animal. Strata
follow their first appearance in the input; rows are reordered within
strata to put the case first. If time is supplied, case times must
increase within each segment. No incomplete stratum or nonfinite value
is silently removed. If you rearrange rows after preparation, prepare
the original data again before fitting. Prepared choice sets must remain
contiguous and case-first.

## Usage

``` r
prepare_ssf_data(
  data,
  id = "ID",
  stratum = "stratum",
  observed = "obs",
  weight = "w",
  time = NULL,
  density_measure = c("unspecified", "area", "step_angle"),
  expected_interval = NULL,
  step_max = Inf
)
```

## Arguments

- data:

  A data frame of observed and control steps.

- id, stratum, observed, weight:

  Names of the segment, choice set, binary observation, and positive
  sampling density columns. Weights may be proportional to the
  conditioned proposal within a choice set when interpreting conditional
  choice probabilities. Absolute emission densities require the
  correctly normalised proposal on the availability domain.

- time:

  Optional name of a numeric or POSIXct observation time column.

- density_measure:

  Measure of the supplied sampling density: `"area"`, `"step_angle"`, or
  `"unspecified"`. This annotation records the convention; it does not
  alter or verify weights.

- expected_interval:

  Optional positive interval between consecutive cases within a segment.
  Requires time. Units are seconds for POSIXct, otherwise the units of
  the numeric time column. Gaps cause an error.

- step_max:

  Upper radial support in step-length units, or Inf for unbounded
  support. Declare the radius of a common availability disk here. Finite
  support currently uses the canonical importance kernel.

## Value

A validated data frame of class `ssf_data`, with standard columns `ID`,
`stratum`, `obs`, and `w`, and a `.source_row` mapping.

## Examples

``` r
d <- data.frame(ID = 1, stratum = rep(1:3, each = 2),
                obs = rep(c(1, 0), 3), w = 1, step = 1:6)
prepare_ssf_data(d)
#>   ID stratum obs w step .source_row  log_step
#> 1  1       1   1 1    1           1 0.0000000
#> 2  1       1   0 1    2           2 0.6931472
#> 3  1       2   1 1    3           3 1.0986123
#> 4  1       2   0 1    4           4 1.3862944
#> 5  1       3   1 1    5           5 1.6094379
#> 6  1       3   0 1    6           6 1.7917595
```
