# Inspect a compound's CAMERA grouping performance

For a given compound name, presents a compound-centric evaluation of
CAMERA's pcgroup assignments using pre-built evaluation objects from the
MS1Match functions. Handles compounds with multiple chromatographic
instances (suffix `_1`, `_2`, etc.) automatically.

Unlike the MS1Match functions which operate across all compounds, this
function filters the pre-built objects to the queried compound and
presents results in a human-readable format. The only computation
performed internally is building the feature-level `long_table`, which
is unique to the per-compound view.

## Usage

``` r
inspect_compound(
  compound_query,
  eval_full,
  anchors,
  compound_eval_wide,
  anno_accuracy,
  purity,
  CAMERA_annotated_pool_FT,
  ppm_tol = 10,
  print_detail = TRUE
)
```

## Arguments

- compound_query:

  `character(1)`. Compound name to search for. Matches exact name or
  name with numeric suffix (e.g. `"L-Norvaline"` matches
  `"L-Norvaline"`, `"L-Norvaline_1"`, `"L-Norvaline_2"`). Special regex
  characters are automatically escaped.

- eval_full:

  `data.frame`. Master evaluation table from
  [`build_eval_full`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_eval_full.md).
  Required columns: `feature_name`, `Compound_name`, `Lib_group`,
  `reference_tier`, `mz`, `rt`, `pcgroup`, `camera_adduct`,
  `adduct_match`.

- anchors:

  `data.frame`. Anchor table from
  [`build_anchors`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anchors.md).
  Required columns: `Lib_group`, `anchor_feature`, `anchor_pcgroup`.

- compound_eval_wide:

  `data.frame`. Per-compound wide evaluation table from
  [`build_compound_eval`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_eval.md).
  Used for tier summary counts.

- anno_accuracy:

  `list`. Output of
  [`build_anno_accuracy`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md).
  The `per_compound` element is used for annotation accuracy counts.

- purity:

  `list`. Output of
  [`build_pcgroup_purity`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_pcgroup_purity.md).
  The `all_pcgroups` element is used for pcgroup context statistics.

- CAMERA_annotated_pool_FT:

  `data.frame`. CAMERA-annotated feature table. Must contain:
  `feature_name`, `pcgroup`. Optionally: `adduct` (or `camera_adduct`),
  `isotopes`.

- ppm_tol:

  `numeric(1)`. ppm tolerance for annotation accuracy check. Default
  `10`.

- print_detail:

  `logical(1)`. If `TRUE` (default), prints a formatted summary to the
  console / Rmd output.

## Value

Invisibly returns a named `list` with three elements:

- `long_table`:

  Feature-level detail table for the queried compound. Columns:
  `Compound_name`, `Lib_group`, `feature_name`, `mz`, `rt`,
  `tier_simple`, `our_annotation`, `camera_pcgroup`, `camera_location`,
  `camera_adduct`, `camera_isotopes`.

- `tier_summary`:

  Filtered rows from `compound_eval_wide` with annotation accuracy
  columns joined from `anno_accuracy$per_compound`.

- `pcgroup_purity`:

  Filtered rows from `purity$all_pcgroups` for the anchor pcgroup(s) of
  this compound.

## Tier definitions

|            |                                                    |
|------------|----------------------------------------------------|
| **Tier**   | **Description**                                    |
| Tier1a     | Confirmed adduct match via direct library matching |
| Tier1a_exp | Adduct match via Phase 2 full-FT expansion         |
| Tier1b     | Isotope peaks (Phase 1)                            |
| Tier1b_exp | Isotope peaks found in Phase 2 expansion           |
| Tier1c     | In-source MS2 fragments                            |
| Tier1d     | Neutral loss pairs                                 |
| Tier2      | Unexplained co-eluting features                    |

## See also

[`build_eval_full`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_eval_full.md),
[`build_anchors`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anchors.md),
[`build_compound_eval`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_compound_eval.md),
[`build_anno_accuracy`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anno_accuracy.md),
[`build_pcgroup_purity`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_pcgroup_purity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
eval_full          <- build_eval_full(anno_result$export_feature_compound_summary,
                                      CAMERA_annotated_pool_FT)
anchors            <- build_anchors(eval_full)
purity             <- build_pcgroup_purity(CAMERA_annotated_pool_FT,
                                           eval_full, anchors)
compound_eval_wide <- build_compound_eval(eval_full, anchors, purity)
anno_accuracy      <- build_anno_accuracy(eval_full, anchors, ppm_tol = 10)

result <- inspect_compound(
  compound_query           = "L-Norvaline",
  eval_full                = eval_full,
  anchors                  = anchors,
  compound_eval_wide       = compound_eval_wide,
  anno_accuracy            = anno_accuracy,
  purity                   = purity,
  CAMERA_annotated_pool_FT = CAMERA_annotated_pool_FT
)

write.csv(result$long_table, "L-Norvaline_detail.csv", row.names = FALSE)
} # }
```
