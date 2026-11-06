# Generate gamma-von Mises control steps with explicit area densities

All observed rows and their segment identifiers are retained. The target
bearing is fixed at the step origin and shared by the case and controls.
Controls are not filtered by a raster or spatial observation window.
This implementation is for an untruncated gamma step proposal and a von
Mises angular proposal. A pool can be nested by retaining alternative
\<= J.

## Usage

``` r
generate_ssf_controls(observed, n_controls, proposal, seed = NULL)
```

## Arguments

- observed:

  One row per observed step, with ID, stratum, step, angle,
  previous_bearing and target_bearing. Bearings and angles are in
  radians.

- n_controls:

  Positive integer number of controls per observed step.

- proposal:

  Named numeric vector with gamma `shape`, gamma `rate`, angular
  `mean_angle` and angular `concentration`.

- seed:

  Optional integer RNG seed, restored on exit.

## Value

Validated `ssf_data` with area-density weights and alternative zero for
the case, 1:n_controls for controls. Source row order is preserved.
