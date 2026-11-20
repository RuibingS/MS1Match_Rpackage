# Classify features by reference-relative evaluation status

Evaluation status is defined *relative to the compound-centric reference
feature set*, not relative to an assumed ground truth. For each
compound, the Tier 1a `[M+H]+` feature is the anchor and defines the
compound's expected group.

- **Ref-confirmed**:

  Feature is assigned to the same group as its compound's anchor.

- **Ref-conflicting**:

  Feature is assigned to a different group than its compound's anchor.

- **Ref-unmatched, co-grouped**:

  Feature has no reference entry but shares a group with a compound's
  anchor (possible over-grouping).

- **Ref-unmatched, unrelated**:

  Feature has no reference entry, and its group contains no compound
  anchor.

- **Ref-unassigned**:

  Feature has a reference entry, but the tool assigned it no group. Zero
  for CAMERA, which groups every feature.

Rows that fit none of these statuses (e.g. a reference feature whose
compound has no resolvable anchor) get `NA` and trigger a warning.

## Usage

``` r
classify_ref_status(eval_full, anchors = NULL)
```

## Arguments

- eval_full:

  `data.frame`. Must contain `Lib_group`, `pcgroup`, `CAMERA_grouped`,
  `CAMERA_split`; ideally `anchor_pcgroup`.

- anchors:

  Optional `data.frame` (output of
  [`build_anchors`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anchors.md)).
  Used to define the set of anchor groups for splitting Ref-unmatched.
  If `NULL`, `eval_full$anchor_pcgroup` is used.

## Value

`eval_full` with an ordered factor column `ref_status`.
