# Prepare observed steps for analytically normalised movement kernels

No controls or proposal weights are required. Each row is one positive
observed step, and each ID is an uninterrupted independent segment.
Directional target effects require the target bearing relative to the
previous heading, chosen at the step origin before observing the
endpoint. The constructor can derive this from target_bearing and
previous_bearing. It cannot establish from numbers alone how a supplied
target was chosen.

## Usage

``` r
prepare_ssf_steps(
  data,
  density_measure,
  id = "ID",
  time = NULL,
  expected_interval = NULL
)
```

## Arguments

- data:

  Data frame with ID, step and angle columns. Angles are radians.
  Optional relative_target_bearing, or target_bearing and
  previous_bearing, describes the prospective target direction. Optional
  target_active is logical; FALSE removes its directional effect when
  direction is undefined.

- density_measure:

  Required: area or step_angle. This determines the reference measure of
  the exact likelihood, not an importance proposal.

- id:

  Name of the independent segment column.

- time:

  Optional numeric or POSIXct time column.

- expected_interval:

  Optional positive interval between consecutive observations in a
  segment, in seconds for POSIXct times.

## Value

An ssf_steps data frame with canonical covariates and metadata.
