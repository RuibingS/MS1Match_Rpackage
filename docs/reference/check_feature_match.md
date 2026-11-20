# Check that MS1Match and a grouping tool used the same feature table

A yes/no check. The answer is `TRUE` only if

1.  every MS1Match reference feature exists in `feature_table` with the
    same m/z, and

2.  every feature id shared by the grouping table and `feature_table`
    has the same m/z (and, where available, the same RT), and at least
    one id is shared.

Feature ids alone are not enough: separate xcms runs all name their
features `FT0001`, `FT0002`, ...

One table may be a subset of the other, e.g. MS1Match run on an
RT-truncated feature table and the grouping tool on the full one. Ids
present in only one of them do not fail the check; they are counted and
excluded by
[`evaluate_grouping`](https://ruibings.github.io/MS1Match_Rpackage/reference/evaluate_grouping.md).

## Usage

``` r
check_feature_match(
  anno_result,
  grouping,
  feature_table,
  mz_tol = 0.001,
  rt_tol = 1
)
```

## Arguments

- anno_result:

  Output of
  [`join_annotation_results()`](https://ruibings.github.io/MS1Match_Rpackage/reference/join_annotation_results.md).

- grouping:

  Output of
  [`prepare_grouping_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/prepare_grouping_table.md)
  or a converter; must contain `mz`.

- feature_table:

  The feature table MS1Match was run on (e.g. `out$ft` for xcms, or
  [`read_mzmine_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_mzmine_feature_table.md)
  for MZmine). Must contain `feature_name` and `mzmed` or `mz`.

- mz_tol:

  Allowed m/z difference in Da, only to absorb rounding (e.g.
  RAMClustR's `mzdec`). Default 0.001.

- rt_tol:

  Allowed RT difference in seconds. Default 1.

## Value

Object of class `feature_table_check`: a list with `same_feature_table`
(`TRUE`/`FALSE`), `summary` (counts) and `mismatches` (one row per
offending feature and reason).
