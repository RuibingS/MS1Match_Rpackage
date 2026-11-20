# Draw one mirror-plot panel (query up, reference down)

Draw one mirror-plot panel (query up, reference down)

## Usage

``` r
plot_mirror(
  df,
  xlim_range = NULL,
  title = "",
  precursor_mz = NULL,
  mz_digits = 1,
  mz_label_mode = "all"
)
```

## Arguments

- df:

  Data frame from
  [`build_mirror_df()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_mirror_df.md)

- xlim_range:

  Length-2 numeric or `NULL`

- title:

  Panel title

- precursor_mz:

  Numeric scalar or `NULL`

- mz_digits:

  Integer

- mz_label_mode:

  One of `"all"`, `"primary"`, `"adduct"`

## Value

A `ggplot` object
