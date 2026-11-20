# Grouping table from a CAMERA peak list

Grouping table from a CAMERA peak list

## Usage

``` r
camera_to_grouping(peaklist, feature_names = NULL)
```

## Arguments

- peaklist:

  Output of
  [`CAMERA::getPeaklist()`](https://rdrr.io/pkg/CAMERA/man/getPeaklist-methods.html)
  (or a saved copy). Must contain `pcgroup`.

- feature_names:

  Optional character vector of xcms feature ids, one per row of
  `peaklist`. If `NULL` and `peaklist` has no `feature_name` column, ids
  are generated in xcms style (`FT0001`, ...), assuming the rows are in
  xcms feature order.

## Value

See
[`prepare_grouping_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/prepare_grouping_table.md).
