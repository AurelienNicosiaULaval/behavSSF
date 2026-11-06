# Exact log normalisers for canonical gamma-von Mises kernels

Integrates exp(beta_step*r + beta_log*log(r) + a*cos(phi) +
b*cos(phi-psi)) on the full unbounded plane or in step-angle
coordinates. The angular concentration is
sqrt((a+b*cos(psi))^2+(b*sin(psi))^2). No endpoint habitat covariates,
truncation or spatial barriers are supported.

## Usage

``` r
ssf_log_normalizer(
  coefficients,
  convention,
  relative_bearing = 0,
  target_active = TRUE
)
```

## Arguments

- coefficients:

  Coefficient matrix with canonical movement rows and an optional
  cos_diff_theta_to_target row.

- convention:

  Required density measure: area or step_angle.

- relative_bearing:

  Prospective target bearings relative to previous headings. Scalar or
  vector of length equal to target_active.

- target_active:

  Logical scalar or vector; FALSE removes the target term.

## Value

Matrix of log normalisers, rows for bearings and columns for states.
