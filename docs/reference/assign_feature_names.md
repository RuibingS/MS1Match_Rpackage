# Give a grouping result the feature ids of the MS1Match feature table

For a grouping result without feature ids. Each row gets the
`feature_name` of the feature in `feature_table` with the same m/z and
RT. Because both come from the same feature table, the values are
identical; the tolerances only absorb rounding. Rows that match no
feature, or more than one, get `NA`; two rows matching the same feature
both get `NA`.

## Usage

``` r
assign_feature_names(
  x,
  feature_table,
  mz_col = "mz",
  rt_col = "rt",
  rt_unit = c("sec", "min"),
  mz_tol = 0.001,
  rt_tol = 1
)
```

## Arguments

- x:

  Grouping result, one row per feature.

- feature_table:

  The feature table MS1Match was run on (`feature_name`, `mzmed`/`mz`,
  `rtmed`/`rt` in seconds).

- mz_col, rt_col:

  Names of the m/z and RT columns in `x`.

- rt_unit:

  Unit of `rt_col`: `"sec"` (default) or `"min"`.

- mz_tol, rt_tol:

  Rounding tolerances (Da, seconds). Defaults 0.001, 1.

## Value

`x` with a `feature_name` column and `rt_col` converted to seconds. A
message reports how many rows were matched.
