# Compare several grouping evaluations side by side

Compare several grouping evaluations side by side

## Usage

``` r
compare_grouping(
  ...,
  metrics = c("pct_no_anchor", "pct_unassigned", "pct_confirmed", "pct_conflicting")
)
```

## Arguments

- ...:

  Objects returned by
  [`evaluate_grouping`](https://ruibings.github.io/MS1Match_Rpackage/reference/evaluate_grouping.md).

- metrics:

  Benchmark columns to show per tool.

## Value

`data.frame`: one row per tier, the reference size, and one column per
tool x metric.
