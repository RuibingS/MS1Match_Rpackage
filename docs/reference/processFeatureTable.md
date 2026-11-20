# Post-process feature table

1.  Moves rownames to a leading `feature_name` column.

2.  Sorts rows by `rtmed`.

3.  Appends `intensity_median` (median of per-sample intensities).

## Usage

``` r
processFeatureTable(ft, x_exp, rt_range = NULL, sort_by = "rtmed")
```

## Arguments

- ft:

  data.frame from
  [`buildFeatureTable()`](https://ruibings.github.io/MS1Match_Rpackage/reference/buildFeatureTable.md).

- x_exp:

  XcmsExperiment used to produce `ft`.

- rt_range:

  Optional numeric vector of length 2 for RT filtering (e.g.
  `c(60, 1200)`).

- sort_by:

  Optional column name to sort by (e.g. `"rtmed"`, `"mzmed"`,
  `"feature_name"`). Defaults to `"rtmed"`.

## Value

Processed data.frame.
