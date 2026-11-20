# Constructor for PeakAnnotationISF

Constructor for PeakAnnotationISF

## Usage

``` r
PeakAnnotationISF(
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

A `PeakAnnotationISF` object.

## See also

[`createPeakAnnotationISF()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createPeakAnnotationISF.md)
