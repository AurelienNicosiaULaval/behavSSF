# Convert a sampling density between planar area and polar step coordinates

Since dx dy = r dr dtheta, q_area = q_step_angle / r. This function
changes weights and their metadata together. It does not change
covariates or fit a model. Use
[`convert_movement_coefficients()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/convert_movement_coefficients.md)
for the corresponding change of kernel coefficients.

## Usage

``` r
convert_density_measure(data, to)
```

## Arguments

- data:

  Validated `ssf_data` with an explicit density measure and finite
  positive `step` values, or observed ssf_steps. For ssf_steps there are
  no proposal weights; only the reference-measure metadata changes.

- to:

  Target measure, `"area"` or `"step_angle"`.

## Value

An `ssf_data` object with converted w and updated density metadata.
