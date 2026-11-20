# Build a flat, all-compounds feature table with CAMERA grouping context

Long table: one row per feature within a compound instance, across EVERY
compound instance present in `eval_full` (not one queried compound at a
time, unlike
[`inspect_compound`](https://ruibings.github.io/MS1Match_Rpackage/reference/inspect_compound.md)).
Reuses the same per-feature anchor-co-grouping logic as
`inspect_compound`'s internal `long_table`, and adds:

- `CAMERA_isotopes`, `CAMERA_adduct`:

  CAMERA's own native annotation columns, kept SEPARATE (not combined
  into one string) so downstream parsing of CAMERA's adduct hypotheses
  (e.g. neutral-mass cross-checks) is not obstructed by a merged column.

- `CAMERA_anno_correct`:

  TRUE/FALSE: does CAMERA's own adduct call on THIS feature match the
  compound's true neutral mass (within `anno_accuracy`'s ppm tolerance),
  established from that compound's own Tier1a/Tier1a_exp evidence.
  Applies across ALL tiers (a Tier1b/1c/1d/Tier2 feature's CAMERA adduct
  call is checked against the same compound-level mass), not only
  Tier1a/Tier1a_exp. `NA` when it cannot be evaluated: CAMERA assigned
  no adduct annotation at all, or the compound has no resolvable neutral
  mass.

- `CAMERA_Anchor_pcgroup_purity`:

  Number of distinct reference compounds whose Tier1a anchor falls in
  this compound's anchor pcgroup. How many OTHER compounds is CAMERA
  merging this one with. Shown only on the row(s) co-grouped with the
  anchor.

- `CAMERA_Anchor_pcgroup_Refconfirmed_pct`:

  This compound's own anchor-co-grouped feature count (X) / total anchor
  pcgroup feature count (Y) x 100. Shown only on the row(s) co-grouped
  with the anchor.

## Usage

``` r
build_compound_instance_table(
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
  The `all_pcgroups` element is used for purity / feature-count context.

- CAMERA_annotated_pool_FT:

  `data.frame`. CAMERA-annotated feature table. Used ONLY as a fallback
  source for pcgroup/adduct/isotopes if `eval_full` is missing those
  native columns – normally `eval_full`'s own `pcgroup`,
  `CAMERA_adduct`, and `isotopes` columns are used directly (the same
  columns
  [`build_anno_accuracy()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md)
  reads), so this table and `CAMERA_anno_correct` can never silently
  disagree for the same row.

- anno_accuracy:

  `list`. Output of
  [`build_anno_accuracy`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md).
  Its `feature_level` element is joined in to populate
  `CAMERA_anno_correct`.

## Value

A `data.frame`, one row per feature within a compound instance, with
columns: `compound_instance`, `MS1Match_lib_group`,
`MS1Match_reference_tier`, `feature_name`, `MS1Match_anno`,
`MS1Match_iso`, `CAMERA_group_status`, `CAMERA_isotopes`,
`CAMERA_adduct`, `CAMERA_anno_correct`, `CAMERA_pcgroup`,
`CAMERA_Anchor_pcgroup`, `CAMERA_Anchor_pcgroup_purity`,
`CAMERA_Anchor_pcgroup_Refconfirmed_pct`.

## See also

[`inspect_compound`](https://ruibings.github.io/MS1Match_Rpackage/reference/inspect_compound.md),
[`build_pcgroup_purity`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_pcgroup_purity.md),
[`build_anno_accuracy`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md)
