# Score and decode a fitted HMM-SSF with daily transitions

Applies a fitted hourly transition schedule to new observations. At the
start of each segment, the hidden state uses the corresponding hour of
the periodic stationary distribution. Step-length autoregression starts
afresh at each segment boundary. All daily lag terms use only
observations in that segment. Fitted hourly effects on mean step length
apply from the first step of each segment.

## Usage

``` r
decode_ssf_diel(data, fit, transition_hour)
```

## Arguments

- data:

  Observed hourly `ssf_steps`.

- fit:

  A fit from
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md)
  with daily transitions.

- transition_hour:

  Integer hour of the step origin, from 0 to 23. Supply one value per
  row in a consistent time zone.

## Value

The log likelihood, Viterbi states, and smoothed state probabilities.
