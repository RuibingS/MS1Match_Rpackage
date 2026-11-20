# Grouping table from an MZmine feature-table export

Converts a grouping stored in an MZmine 3 full feature table (e.g. ion
identity networking) to the format used by
[`evaluate_grouping`](https://ruibings.github.io/MS1Match_Rpackage/reference/evaluate_grouping.md).
Applies the same conversions as
[`read_mzmine_feature_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_mzmine_feature_table.md):
MZmine `id` becomes `feature_name = sprintf("FT%04d", id)` and RT is
converted from minutes to seconds.

## Usage

``` r
mzmine_to_grouping(
  x,
  group_col,
  keep_features = NULL,
  mz_col = "mz",
  rt_col = "rt",
  unassigned = ""
)
```

## Arguments

- x:

  `data.frame` read from the MZmine CSV
  (`read.csv(file, check.names = FALSE)`). Must contain `id` and
  `group_col`.

- group_col:

  Name of the column holding the group id (check your export, e.g. the
  ion identity network id column).

- keep_features:

  Optional character vector of `feature_name` values to keep, typically
  `read_mzmine_feature_table(file)$feature_name`, so that features
  removed by the import's detection filter are dropped.

- mz_col:

  Name of the m/z column. Default `"mz"`.

- rt_col:

  Name of the RT column in minutes. Default `"rt"`; set to `NULL` if
  absent.

- unassigned:

  Values of `group_col` meaning "not grouped". Default `""`: MZmine
  numbers its groups from 0, so 0 is a real group.

## Value

See
[`prepare_grouping_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/prepare_grouping_table.md).
