# Join lists of annotation results to library match output

Combine direct library match results (MS1/MS2) with downstream
annotation results (adduct/isotope/MS2 fragment/NL) onto the feature
table.

Tier naming convention:

- **1a** — adduct match (Phase 1, within library group)

- **1b** — isotope (Phase 1, within library group)

- **1c** — MS2 fragment match (Phase 1: feature already in a Phase 1
  library group)

- **1d** — neutral loss (Phase 1, pairwise within group)

- **1a_exp / 1b_exp / 1c_exp / 1d_exp** — Phase 2 expansions (feature
  NOT assigned in Phase 1; derived from a Phase 1 confirmed M)

- **2** — unexplained (in group, no rule fires)

Tiers are mutually exclusive — one tier per feature, highest priority
wins: **1a \> 1b \> 1c \> 1d \> 2 \> 1a_exp \> 1b_exp \> 1c_exp \>
1d_exp**. `Tier2` (Unexplained) is a Phase 1 assignment and therefore
outranks all Phase 2 expansions.

The `annotation_phase` column takes exactly two values: `"Phase1"` (any
Tier 1a/1b/1c/1d or Tier 2 assignment) and `"Phase2"` (any `_exp` tier).
Sub-labels such as `Phase1_MS2frag` or `Phase2_IsoExpansion` are not
used; `reference_tier` already encodes the specific evidence type.

## Usage

``` r
join_annotation_results(
  FT,
  direct_lib_list,
  lib_anno_list = NULL,
  target_match_list = NULL,
  isf_result = NULL,
  include_phase2 = TRUE,
  verbose = TRUE
)
```

## Arguments

- FT:

  Required data.frame. Feature table.

- direct_lib_list:

  Required list.

- lib_anno_list:

  Optional list/data.frame. Output of annotate().

- target_match_list:

  Optional list/data.frame. Target list annotations.

- isf_result:

  Optional data.frame. Output of `annotateISF()` /
  [`run_isf_annotation()`](https://ruibings.github.io/MS1Match_Rpackage/reference/run_isf_annotation.md).
  Adds Tier 1c / 1c_exp columns and resolves `reference_tier` by
  priority.

- include_phase2:

  Logical. If `FALSE`, Phase 2 expansion rows (`Tier1a_exp`,
  `Tier1b_exp`, `Tier1d_exp`) from `lib_anno_list` are excluded from the
  output. Default `TRUE`.

- verbose:

  Logical. If `FALSE`, suppresses all progress messages. Default `TRUE`.
