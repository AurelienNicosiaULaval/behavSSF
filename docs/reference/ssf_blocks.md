# Make consecutive validation blocks without splitting segments

Input order defines chronology. Cuts are placed near equally spaced
cumulative step counts, subject to leaving at least one segment per
block.

## Usage

``` r
ssf_blocks(data, n_blocks = 3L)
```

## Arguments

- data:

  Prepared `ssf_data` or `ssf_steps`.

- n_blocks:

  Number of validation blocks, at least two.

## Value

An integer vector named by segment ID, for
[`compare_ssf_states()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_states.md).
