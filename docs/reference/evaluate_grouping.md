# Evaluate a feature-grouping tool against the MS1Match reference set

For each library compound, the Tier 1a `[M+H]+` feature is the anchor
and the group the tool assigns to it defines the compound's expected
group. Each reference feature (Phase 1 and Phase 2) gets the first
status that applies: Ref-unassigned (no-anchor) if the tool did not
group the anchor, Ref-unassigned if the anchor is grouped but the
feature is not, otherwise Ref-confirmed or Ref-conflicting. Every other
feature of the tool is Ref-unmatched, co-grouped (in a group with a
compound anchor) or Ref-unmatched, unrelated (in a group without an
anchor, or in no group).

## Usage

``` r
evaluate_grouping(
  anno_result,
  grouping,
  tool_name = "tool",
  rt_window = NULL,
  feature_table,
  mz_tol = 0.001,
  rt_tol = 1
)
```

## Arguments

- anno_result:

  Output of
  [`join_annotation_results()`](https://ruibings.github.io/MS1Match_Rpackage/reference/join_annotation_results.md).
  Uses `export_feature_compound_summary` (reference rows) and
  `export_reduced_anno_table` (anchor search).

- grouping:

  Output of
  [`prepare_grouping_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/prepare_grouping_table.md),
  [`camera_to_grouping`](https://ruibings.github.io/MS1Match_Rpackage/reference/camera_to_grouping.md)
  or
  [`ramclustr_to_grouping`](https://ruibings.github.io/MS1Match_Rpackage/reference/ramclustr_to_grouping.md).

- tool_name:

  Label used in printed output and in
  [`compare_grouping`](https://ruibings.github.io/MS1Match_Rpackage/reference/compare_grouping.md).

- rt_window:

  Optional `c(min, max)` in seconds. Applied only to features outside
  the reference set (Ref-unmatched), so their counts are comparable
  between tools.

- feature_table:

  The feature table MS1Match was run on. The function stops unless
  [`check_feature_match`](https://ruibings.github.io/MS1Match_Rpackage/reference/check_feature_match.md)
  confirms that the grouping and MS1Match used the same table. If one is
  a subset of the other, only features present in both are evaluated.

- mz_tol, rt_tol:

  Passed to
  [`check_feature_match`](https://ruibings.github.io/MS1Match_Rpackage/reference/check_feature_match.md).

## Value

Object of class `grouping_evaluation`: a list with `tool`, `benchmark`
(tier table; `pct_grouped` is the share of reference features the tool
grouped at all, i.e. its coverage), `status_counts`, `full_table` (one
row per reference feature x compound, plus one per feature outside the
reference set), `anchors` and `feature_match`.

## Details

The grouping table must list *every* feature the tool processed,
including those it did not group (`group = NA`). A feature missing from
the grouping table is treated as outside the shared feature table and
excluded, not counted as Ref-unassigned.
