# Assess Tier 1a adduct annotation accuracy

For each Tier 1a (and optionally Tier 1a_exp) feature, checks whether
CAMERA's adduct neutral mass assignment matches the MS1Match's annotated
neutral mass within a given ppm tolerance.

## Usage

``` r
build_anno_accuracy(eval_full, anchors, ppm_tol = 10, mda_tol = 5)
```

## Arguments

- eval_full:

  `data.frame`. Output of
  [`build_eval_full`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_eval_full.md).
  Must contain: `reference_tier`, `CAMERA_adduct`, and at least one of
  `grp_anno_match` or `exp_anno_match`.

- anchors:

  `data.frame`. Output of
  [`build_anchors`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anchors.md).

- ppm_tol:

  `numeric(1)`. ppm tolerance for neutral mass matching. Default 10.

- mda_tol:

  `numeric(1)`. Absolute mass tolerance floor in millidaltons. A match
  is accepted if either the ppm difference is within `ppm_tol` OR the
  absolute difference is within `mda_tol` mDa. Default 5 mDa, which
  covers the worst-case truncation error from CAMERA reporting neutral
  masses to 3 decimal places.

## Value

A named `list`:

- `per_compound`:

  Per `Lib_group` annotation accuracy counts (Tier1a/Tier1a_exp only —
  the established methodology), with columns `n_total`, `n_unannotated`,
  `n_annotated`, `n_correct`, `n_wrong`, `coverage_pct`, `accuracy_pct`.

- `summary`:

  Two-row summary table (one per tier) plus a TOTAL row, with the same
  columns. Unaffected by the all-tier extension to `feature_level`
  below.

- `feature_level`:

  Row-level table covering EVERY tier (not only Tier1a/Tier1a_exp), one
  row per feature x compound instance, with `feature_name`, `Lib_group`,
  `Compound_name`, `tier_simple`, `has_CAMERA_annotation`, and
  `annotation_match`. Every feature belonging to a compound shares that
  compound's one true neutral mass (established from its own
  Tier1a/Tier1a_exp evidence), so CAMERA's adduct call on ANY tier's
  feature can be checked against it. `annotation_match` is `NA` when it
  cannot be evaluated (no CAMERA annotation present, or the compound has
  no resolvable neutral mass at all), and TRUE/FALSE only when CAMERA's
  call was actually checked — joinable into
  [`build_compound_instance_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_instance_table.md).

## Details

1.  **Tier-aware adduct column**: Tier1a features are evaluated against
    `grp_anno_match`; Tier1a_exp features against `exp_anno_match`.

2.  **Multi-annotation aware**: CAMERA sometimes assigns multiple adduct
    annotations to one feature in a single string, e.g.
    `"[M+K+HCOOH]+ 63.128 [M+Na+HCOOH]+ 79.102 [M+H+HCOOH]+ 101.084"`.The
    updated code splits the string on adduct boundaries and checks
    whether *any* of CAMERA's candidate neutral mass values matches the
    MS1Match neutral mass within `ppm_tol`.
