# Standardize a grouping tool's output for evaluation

Reduces a feature-level grouping result to the columns used by
[`evaluate_grouping`](https://ruibings.github.io/MS1Match_Rpackage/reference/evaluate_grouping.md):
`feature_name`, `group`, `mz` and `rt`. The grouping tool must have been
run on the same feature table as MS1Match; this function only renames
columns and cannot map features between different feature tables.

## Usage

``` r
prepare_grouping_table(
  x,
  group_col,
  feature_col = "feature_name",
  mz_col = NULL,
  rt_col = NULL,
  rt_unit = c("sec", "min"),
  unassigned = c("0", "")
)
```

## Arguments

- x:

  `data.frame` with one row per feature.

- group_col:

  Name of the column holding the group / cluster id.

- feature_col:

  Name of the column holding the feature id (same ids as in the feature
  table MS1Match was run on). Plain numbers (1, 2, ...) are converted to
  `FT0001`, `FT0002`, ..., as in
  [`read_feature_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_feature_table.md).
  Default `"feature_name"`.

- mz_col:

  Name of the m/z column. Needed to verify that the tool used the same
  feature table as MS1Match.

- rt_col:

  Optional name of the retention-time column.

- rt_unit:

  Unit of `rt_col`: `"sec"` (default) or `"min"`. Minutes are converted
  to seconds.

- unassigned:

  Values of `group_col` that mean "not grouped" (in addition to `NA`).
  Default `c("0", "")`.

## Value

`data.frame` with columns `feature_name`, `group` (character, `NA` =
unassigned), `mz` and `rt`.
