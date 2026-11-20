# Convert long annotations to wide or collapsed format

Convert long annotations to wide or collapsed format

## Usage

``` r
long_to_wide.fun(df, col_name1, col_name2, col_name3, collapse = TRUE)
```

## Arguments

- df:

  A data frame with annotations.

- col_name1:

  The feature name column (e.g., Feature_name).

- col_name2:

  The compound name column (e.g., Comp_name_PAMDB).

- col_name3:

  The ion name + exact mass column (e.g., PI_name_PAMDB).

- collapse:

  Logical. If TRUE (default), collapse all annotations into one column
  separated by "; ".

## Value

A data frame with wide or collapsed annotations.
