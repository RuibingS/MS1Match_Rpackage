# Run MS2 fragment annotation (Phase 1 + Phase 2)

R wrapper around `isf_annotate_cpp()`. Identifies features whose `mzmed`
matches a fragment ion of a confirmed precursor.

## Usage

``` r
run_isf_annotation(
  phase1_result,
  ms2_list,
  FT,
  ppm_isf = 10,
  mz_threshold_isf = 0.005,
  rt_isf_threshold = 30
)
```

## Arguments

- phase1_result:

  `data.frame`. Output of
  [`annotateLibMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/annotateLibMatch.md).
  Must contain columns `feature_name`, `Lib_group`, `reference_tier`,
  and `Compound_M`.

- ms2_list:

  `list`. From
  [`read_ms2_consensus()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_ms2_consensus.md).

- FT:

  `data.frame`. Full MS1 feature pool.

- ppm_isf:

  Numeric. ppm tolerance for fragment matching. Default `10`.

- mz_threshold_isf:

  Numeric. Absolute m/z floor (Da). Default `0.005`.

- rt_isf_threshold:

  Numeric. RT window (s). Default `30`.

## Value

`data.frame` — see `annotateISF()`.
