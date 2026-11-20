# Build the feature × compound evaluation table

Joins the MS1Match's tier-annotated feature-compound summary onto the
CAMERA annotated feature table, so every row carries both the MS1Match
annotation tier and the CAMERA pcgroup / adduct assignments.

## Usage

``` r
build_eval_full(feature_compound_summary, CAMERA_ft)
```

## Arguments

- feature_compound_summary:

  `data.frame`. The `export_feature_compound_summary` element of
  [`join_annotation_results()`](https://ruibings.github.io/MS1Match_Rpackage/reference/join_annotation_results.md).
  Must contain: `feature_name`, `Compound_name`, `reference_tier`,
  `annotation_phase`, `Lib_group`, `adduct_match`, `rt`.

- CAMERA_ft:

  `data.frame`. CAMERA-annotated feature table (e.g.
  `CAMERA_annotated_pool_FT_sub`). Must contain: `feature_name`,
  `pcgroup`, `adduct` (CAMERA's adduct/isotope string), `mz`.

## Value

A `data.frame` with one row per feature × compound combination,
containing all MS1Match annotation columns plus CAMERA's `pcgroup`,
`adduct` (renamed `CAMERA_adduct`), and `mz`.
