# Validate and organise case-control step data

Each ID must identify an independent, uninterrupted segment. Supply
separate IDs across observation gaps, even for the same animal. Strata
follow their first appearance in the input; rows are reordered within
strata to put the case first. If time is supplied, case times must
increase within each segment. No incomplete stratum or nonfinite value
is silently removed.

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
  expected_interval = NULL
)
```

## Arguments

- data:

  A data frame of observed and control steps.

- id, stratum, observed, weight:

  Names of the segment, choice set, binary observation, and positive
  sampling density columns.

- time:

  Optional name of a numeric or POSIXct observation time column.

- density_measure:

  Measure of the supplied sampling density: `"area"`, `"step_angle"`, or
  `"unspecified"`. This annotation does not alter weights.

- expected_interval:

  Optional positive interval between consecutive cases within a segment.
  Requires time. Units are seconds for POSIXct, otherwise the units of
  the numeric time column. Gaps cause an error.

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
