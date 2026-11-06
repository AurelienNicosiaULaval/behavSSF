# Decode independent hidden Markov sequences in log space

Uses a homogeneous transition matrix. Restarts both filtering and
Viterbi decoding at every segment boundary, including segments of length
one. Ties in Viterbi decoding are resolved by the first state index.

## Usage

``` r
decode_states(
  log_emissions,
  transition,
  segment = rep(1L, nrow(as.matrix(log_emissions))),
  initial = NULL
)
```

## Arguments

- log_emissions:

  Numeric matrix with one row per step and one column per state.

- transition:

  Row-stochastic transition probability matrix.

- segment:

  Independent sequence identifiers, one per step.

- initial:

  Initial state probabilities, used at each boundary. If NULL, the
  stationary distribution is used; it must be uniquely computable.

## Value

A list with `log_likelihood`, global Viterbi `state`, and smoothed
marginal state `posterior`. These probabilities condition on fixed
parameters.

## Examples

``` r
decode_states(log(matrix(c(.8, .2, .3, .7), 2, byrow = TRUE)),
              matrix(c(.9, .1, .1, .9), 2, byrow = TRUE))
#> $log_likelihood
#> [1] -1.599488
#> 
#> $state
#> [1] 1 1
#> 
#> $posterior
#>           [,1]      [,2]
#> [1,] 0.6732673 0.3267327
#> [2,] 0.5495050 0.4504950
#> 
```
