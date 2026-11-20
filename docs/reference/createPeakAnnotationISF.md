# Create a PeakAnnotationISF object

User-facing constructor for
[`PeakAnnotationISF-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationISF-class.md).

Tier naming:

- **Tier 1c** — MS2 fragment match, intra-group (Phase 1)

- **Tier 1c_ext** — MS2 fragment match, RT-gated features (Phase 2)

Priority: 1a \> 1b \> 1c \> 1d \> 2. MS2 tier (1c) only assigned when no
higher-priority tier already fired.

Parent features are identified by `parent_col`/`parent_val`:

- `parent_val = NULL` (default): any non-empty, non-NA string in
  `parent_col` is a parent (e.g. `PI_match` column containing compound
  name strings).

- Otherwise: exact match `parent_col == parent_val`.

## Usage

``` r
createPeakAnnotationISF(
  phase1_result,
  ms2_list,
  FT,
  parent_col = "PI_match",
  parent_val = NULL,
  ppm_isf = 10,
  mz_threshold_isf = 0.005,
  rt_isf_threshold = 5
)
```

## Arguments

- phase1_result:

  `data.frame`. Output of
  [`annotateLibMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/annotateLibMatch.md).
  Required columns: `feature_name`, `Lib_group`, `mz`, and `parent_col`.
  Optional: `rt`.

- ms2_list:

  `list`. From
  [`read_ms2_consensus()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_ms2_consensus.md).

- FT:

  `data.frame`. Full MS1 feature tabel. Required: `feature_name`,
  `mzmed`/`mz`. Optional: `rtmed`/`rt`.

- parent_col:

  Character. Column identifying parent features. Default `"PI_match"`.

- parent_val:

  Value marking a parent. `NULL` = any non-empty string. Default `NULL`.

- ppm_isf:

  Numeric. ppm tolerance. Default `10`.

- mz_threshold_isf:

  Numeric. Absolute m/z floor (Da). Default `0.005`.

- rt_isf_threshold:

  Numeric. RT window (seconds). Default `5`.

## Value

A
[`PeakAnnotationISF-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationISF-class.md)
object.

## See also

`annotateISF()`,
[`PeakAnnotationISF-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationISF-class.md),
[`read_ms2_consensus()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_ms2_consensus.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ms2_list <- read_ms2_consensus(ms2_consensus, feature_table = FT)

isf_obj <- createPeakAnnotationISF(
  phase1_result    = pool_annotate_res_ms1,
  ms2_list         = ms2_list,
  FT          = FT,
  parent_col       = "PI_match",
  ppm_isf          = 10,
  rt_isf_threshold = 5.0
)
isf_result <- annotateISF(isf_obj)
} # }
```
