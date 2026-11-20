# Flatten DirectLibMatch results to a data.frame

Convert the `results` slot (a list of per-feature/per-hit lists) from a
[DirectLibMatch](https://ruibings.github.io/MS1Match_Rpackage/reference/DirectLibMatch-class.md)
object into a single `data.frame`. Optionally collapse vector fields
`Compound_peak_mz_full` and `Compound_peak_int_full` into
comma-separated strings, either overwriting the originals or writing to
new `*_str` columns.

## Usage

``` r
results_df(object, collapse_full = TRUE, remove_full = TRUE, suffix = "_str")
```

## Arguments

- object:

  A
  [DirectLibMatch](https://ruibings.github.io/MS1Match_Rpackage/reference/DirectLibMatch-class.md)
  object with `results` populated (e.g., after calling
  [`libAnnotate()`](https://ruibings.github.io/MS1Match_Rpackage/reference/libAnnotate.md)).

- collapse_full:

  Logical; collapse `*_full` vectors to strings? Default is `TRUE`.

- remove_full:

  Logical; if `TRUE`, remove the original `*_full` columns and write
  collapsed strings to new columns with `suffix`. If `FALSE`, overwrite
  the original `*_full` with the collapsed strings. Default is `FALSE`.

- suffix:

  Character; suffix for new columns when `remove_full = TRUE`. Default
  is `"_str"`.

## Value

A `data.frame` with one row per element of `object@results`.

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- DirectLibMatch(
  FT = result_PA14_FT_filtered,
  compounds_list = PA14_CBIO_library_MS1,
  diff_mz_thr = 0.005, diff_rt_thr = 60, ppm = 10,
  ms_level = "MS1", polarity = "pos", results = list()
)
obj <- libAnnotate(obj)
df  <- results_df(obj)
df2 <- results_df(obj, remove_full = TRUE)
as.data.frame(obj)
} # }
```
