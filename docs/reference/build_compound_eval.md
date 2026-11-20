# Build the per-compound tier-wise grouping recall table

For each compound instance × annotation tier, counts how many features
were co-grouped with the `[M+H]+` anchor (recall) versus split to a
different pcgroup.

## Usage

``` r
build_compound_eval(eval_full, anchors, purity)
```

## Arguments

- eval_full:

  `data.frame`. Output of
  [`build_eval_full`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_eval_full.md).

- anchors:

  `data.frame`. Output of
  [`build_anchors`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anchors.md).

- purity:

  `list`. Output of
  [`build_pcgroup_purity`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_pcgroup_purity.md).

## Value

A wide `data.frame` with one row per `Lib_group`, with columns
`n_MS1Match_<Tier>`, `n_CAMERA_grouped_<Tier>`, `n_CAMERA_split_<Tier>`
for each tier.
