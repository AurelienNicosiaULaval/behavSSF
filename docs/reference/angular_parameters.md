# Combine directional persistence and a target-direction effect

For turning angle phi and target bearing psi relative to the previous
heading, the angular predictor is a*cos(phi) + b*cos(phi-psi). Its
conditional von Mises resultant is (a+b*cos(psi), b*sin(psi)). Thus
abs(a) alone is the intrinsic persistence concentration, not the full
concentration when a target effect is present.

## Usage

``` r
angular_parameters(persistence, target_bias = 0, relative_bearing = 0)
```

## Arguments

- persistence:

  Signed coefficient a for cos(turning angle).

- target_bias:

  Signed coefficient b for direction toward the target.

- relative_bearing:

  Target bearing minus previous heading, in radians.

## Value

A data frame with `mean_angle` and nonnegative `concentration`. The mean
is NA for a uniform angular distribution. Arguments must have a common
length or length one.

## Examples

``` r
angular_parameters(2, 1, pi / 2)
#>   mean_angle concentration
#> 1  0.4636476      2.236068
```
