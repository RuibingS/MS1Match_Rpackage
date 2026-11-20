# Read a pre-processed feature table for MS1Match

Imports a feature table from any software (CSV file or `data.frame`)
into the format MS1Match uses: `feature_name`, `mzmed`, `rtmed`
(seconds), `intensity_mean`, `intensity_median`, followed by one column
per sample. A group / cluster column is not required. Any other columns
(e.g. such a group index) are kept unchanged.

## Usage

``` r
read_feature_table(
  x,
  mz_col = "mz",
  rt_col = "rt",
  rt_unit = c("sec", "min"),
  sample_cols,
  id_col = NULL
)
```

## Arguments

- x:

  Path to a CSV file, or a `data.frame`.

- mz_col, rt_col:

  Names of the m/z and RT columns.

- rt_unit:

  Unit of `rt_col`: `"sec"` (default) or `"min"`. Converted to seconds.

- sample_cols:

  Per-sample intensity columns: a character vector of column names, or a
  single regular expression matched against the column names (e.g.
  `"\\.mzML$"`).

- id_col:

  Optional feature-id column, used as `feature_name`. Plain numbers (1,
  2, ...) are converted to `FT0001`, `FT0002`, ... other ids are kept as
  they are. If `NULL`, ids `FT0001`, ... are assigned after sorting by
  RT then m/z.

## Value

`data.frame`, one row per retained feature.
