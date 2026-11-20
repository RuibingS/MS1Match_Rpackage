# Define the M+H+ anchor feature for each compound instance

@description For each `Lib_group` (compound instance), selects the
single Tier 1a `[M+H]+` feature as the reference anchor. Where multiple
Tier 1a `[M+H]+` features exist for one instance, the one with the
highest intensity (or first by feature name) is kept.

## Usage

``` r
build_anchors(eval_full, full_anno = NULL)
```

## Arguments

- eval_full:

  `data.frame`. Output of
  [`build_eval_full`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_eval_full.md).
  Must contain: `Lib_group`, `reference_tier`, `adduct_match`,
  `pcgroup`, `feature_name`.

## Value

A `data.frame` with one row per `Lib_group`, containing: `Lib_group`,
`Compound_name`, `anchor_feature`, `anchor_pcgroup`, `anchor_mz`,
`anchor_rt`.
