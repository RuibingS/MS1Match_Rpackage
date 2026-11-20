# Build the tier-stratified benchmark summary table

Aggregates grouping recall (from `build_compound_eval`) and annotation
accuracy (from `build_anno_accuracy`) into a single table, one row per
annotation tier.

## Usage

``` r
build_benchmark_table(compound_eval_wide, anno_accuracy, purity)
```

## Arguments

- compound_eval_wide:

  `data.frame`. Output of
  [`build_compound_eval`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_eval.md).

- anno_accuracy:

  `list`. Output of
  [`build_anno_accuracy`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md).

- purity:

  `list`. Output of
  [`build_pcgroup_purity`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_pcgroup_purity.md).

## Value

A `data.frame` with columns: `Tier`, `n_MS1Match`, `n_CAMERA_grouped`,
`recall_pct`, `n_CAMERA_split`, `split_pct`,
`CAMERA_annotation_accuracy_pct`, `CAMERA_has_concept`.
