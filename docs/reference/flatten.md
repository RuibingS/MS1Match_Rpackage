# Flatten matched PI result list into a data.frame

Filters out empty (unmatched) feature slots and row-binds the remainder
into a single flat `data.frame`. The result includes the new columns
added by the updated Rcpp function: `exactMass`, `Comp_mz`, `Comp_rt`,
`mz_diff`, and `ppm_error`.

## Usage

``` r
flatten(object)
```

## Arguments

- object:

  A `DirectTargetMatch` object, or a raw `list` as returned by
  `PImatch_target_rcpp()`

## Value

A `data.frame` of all matched feature–compound pairs, or an empty
`data.frame` if nothing matched
