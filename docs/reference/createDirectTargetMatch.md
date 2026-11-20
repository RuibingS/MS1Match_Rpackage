# Create a DirectTargetMatch object

Constructor for the `DirectTargetMatch` S4 class. Validates inputs and
provides sensible defaults before calling `new()`.

## Usage

``` r
createDirectTargetMatch(
  FT,
  Comp_data,
  PIon = "[M+H]+",
  diff_mz_thr = 0.005,
  ppm = 10,
  rt_thr_compound = 30
)
```

## Arguments

- FT:

  A `data.frame` of detected features. Must contain a `feature_name`
  column and an m/z column (`mzmed`, `mz`, or `m/z`). An RT column
  (`rtmed` or `rt`) is optional but enables RT gating when the target
  list also has RT values.

- Comp_data:

  A `data.frame` target list. Accepted formats: HMDB (`MONO_MASS`,
  `[M+H]+`, `[M-H]-`), plate-based (`Neutral mass`, `M+H`, `M-H`,
  `CNAME`), PA14 (`exactMass`, `[M+H]+`, `[M-H]-`, `RT`), or any format
  supported by
  [`normalize_target_list()`](https://ruibings.github.io/MS1Match_Rpackage/reference/normalize_target_list.md).

- PIon:

  Primary ion: `"[M+H]+"` (positive mode) or `"[M-H]-"` (negative mode).
  Default `"[M+H]+"`.

- diff_mz_thr:

  Numeric. Absolute m/z tolerance in Da. Default `0.005` (5 mDa).

- ppm:

  Numeric. ppm m/z tolerance. Default `10`.

- rt_thr_compound:

  Numeric. RT window in seconds for feature vs compound RT gating.
  Default `30`. Set to `Inf` to disable entirely.

## Value

A `DirectTargetMatch` S4 object ready to pass to
[`PImatch_target()`](https://ruibings.github.io/MS1Match_Rpackage/reference/PImatch_target.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Works with any supported target list format
obj <- createDirectTargetMatch(
  FT              = xcms_output_ft,
  Comp_data       = plate_target_list,  # or PA14_target_list,  hmdb_target_list
  PIon            = "[M+H]+",
  diff_mz_thr     = 0.005,
  ppm             = 10,
  rt_thr_compound = 30
)
obj <- PImatch_target(obj)
df  <- flatten(obj)
} # }
```
