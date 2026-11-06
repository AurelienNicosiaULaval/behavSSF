# Convert movement-kernel coefficients to another density measure

A radial-angle density equals the planar density multiplied by r, so its
log_step coefficient is larger by one. All other coefficients are
unchanged. This conversion preserves the underlying spatial model when
proposal densities are also converted with
[`convert_density_measure()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/convert_density_measure.md).

## Usage

``` r
convert_movement_coefficients(coefficients, from, to)
```

## Arguments

- coefficients:

  Numeric coefficient matrix with a unique log_step row.

- from, to:

  Source and target measures, `"area"` or `"step_angle"`.

## Value

A coefficient matrix in the target measure.
