# Build the complete CAMERA feature table with reference-relative status

Extends
[`build_compound_instance_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_instance_table.md)
to cover every CAMERA feature in `CAMERA_annotated_pool_FT`, not only
reference-bearing ones. Reference-bearing rows are produced exactly as
`build_compound_instance_table` would; every CAMERA feature with no
reference entry at all is added as an additional row. Every row — both
kinds — carries a single `ref_status` column using the our definitions ,
computed identically to (and consistent with)
[`classify_ref_status`](https://ruibings.github.io/MS1Match_Rpackage/reference/classify_ref_status.md):

- **Ref-confirmed**:

  Co-grouped with the compound's anchor.

- **Ref-conflicting**:

  Reference entry exists, split from anchor.

- **Ref-unmatched, co-grouped**:

  No reference entry, but sits in a pcgroup that contains some OTHER
  compound's anchor — diagnostic of over-grouping.

- **Ref-unmatched, unrelated**:

  No reference entry, and its pcgroup contains no anchor at all — a
  reference-set coverage gap, not a CAMERA grouping outcome.

- **Ref-unassigned**:

  Reference entry exists, but the tool assigned the feature no group
  (zero for CAMERA).

Every row has `mz` and `rt`: reference rows from MS1Match, the other
rows from the CAMERA peak list (`mz`, `rt`).

For Ref-unmatched rows, `compound_instance`, `MS1Match_lib_group`,
`MS1Match_reference_tier`, `MS1Match_anno`, and `MS1Match_iso` are `NA`
(no compound owns the row). `CAMERA_Anchor_pcgroup_Refconfirmed_pct` and
`CAMERA_anno_correct` are also left `NA` for these rows by design: the
former is a specific compound's own share of its anchor pcgroup, the
latter requires an MS1Match-expected mass to check CAMERA's call against
— neither exists when no compound owns the row.

## Usage

``` r
build_full_feature_table(
  eval_full,
  anchors,
  purity,
  CAMERA_annotated_pool_FT,
  anno_accuracy
)
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

- CAMERA_annotated_pool_FT:

  `data.frame`. The FULL CAMERA-annotated feature table (every feature,
  not pre-filtered to reference matches).

- anno_accuracy:

  `list`. Output of
  [`build_anno_accuracy`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md).
  Passed through to
  [`build_compound_instance_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_instance_table.md)
  to populate `CAMERA_anno_correct`.

## Value

A `data.frame`: all `build_compound_instance_table` columns (including
`CAMERA_isotopes`, `CAMERA_adduct` as separate columns and
`CAMERA_anno_correct`), plus `ref_status` (ordered factor, 5 levels),
covering every feature in `CAMERA_annotated_pool_FT`.

## See also

[`build_compound_instance_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_instance_table.md),
[`classify_ref_status`](https://ruibings.github.io/MS1Match_Rpackage/reference/classify_ref_status.md),
[`build_benchmark_table_from_full`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_benchmark_table_from_full.md)
