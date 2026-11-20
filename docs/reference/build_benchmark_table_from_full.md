# Recompute the tier-stratified benchmark table from build_full_feature_table()

Produces the identical column layout as
[`build_benchmark_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_benchmark_table.md),
but every count (`n_MS1Match`, `n_CAMERA_grouped`/recall,
`n_CAMERA_split`/split, `n_CAMERA_unassigned`) is aggregated directly
from `ref_status` on `full_table` — the same row-level object exported
for manual inspection — rather than from the separately computed
`compound_eval_wide` pathway. Running both
[`build_benchmark_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_benchmark_table.md)
and this function and comparing output is a useful consistency check:
any divergence indicates the two pipelines have drifted out of sync
somewhere upstream.

Annotation coverage/accuracy (`CAMERA_coverage_pct`,
`CAMERA_annotation_accuracy_pct`) are NOT re-derived here — they depend
on a ppm-tolerance match between CAMERA's own adduct call and the
expected mass, a different computation axis from `ref_status` entirely —
and are still sourced from `anno_accuracy$summary`, same as
[`build_benchmark_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_benchmark_table.md).

## Usage

``` r
build_benchmark_table_from_full(full_table, anno_accuracy)
```

## Arguments

- full_table:

  `data.frame`. Output of
  [`build_full_feature_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_full_feature_table.md).
  Must contain `ref_status` and `MS1Match_reference_tier`.

- anno_accuracy:

  `list`. Output of
  [`build_anno_accuracy`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md).

## Value

A `data.frame`, one row per tier plus TOTAL, with columns: `Tier`,
`n_MS1Match`, `n_CAMERA_grouped`, `recall_pct`, `n_CAMERA_split`,
`split_pct`, `CAMERA_coverage_pct`, `CAMERA_annotation_accuracy_pct`,
`CAMERA_has_concept`, `n_CAMERA_unassigned`, `unassigned_pct`.

## See also

[`build_full_feature_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_full_feature_table.md),
[`build_benchmark_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_benchmark_table.md)
