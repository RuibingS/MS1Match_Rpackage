# Build per-pcgroup purity summary for the over-grouping scatter plot

Extracts neutral masses from CAMERA's adduct annotation column for every
feature, then summarises per pcgroup: how many distinct neutral masses
were assigned by CAMERA within that pcgroup, and what is the median RT.
The result is the `plot2_data` object consumed by the
`plot-over-grouping` chunk in `CAMERA_evaluation_analysis.Rmd`.

CAMERA writes adduct annotations as e.g. `"[M+H]+ 194.057909"` or
`"[M+Na]+ 194.057909; [M+K]+ 194.057909"`. This function extracts the
trailing numeric mass (the last decimal number on each semicolon-
delimited entry), rounds to 4 dp to absorb floating point noise, and
counts distinct values per pcgroup.

## Usage

``` r
build_plot2_data(CAMERA_ft, mass_digits = 4L)
```

## Arguments

- CAMERA_ft:

  `data.frame`. CAMERA-annotated feature table. Must contain: `pcgroup`,
  `rt`, and a CAMERA adduct column. The adduct column is detected
  automatically: the first column whose name matches `adduct`
  (case-insensitive) is used, falling back to a column named
  `CAMERA_adduct`. Rows with missing or empty adduct strings are
  excluded.

- mass_digits:

  `integer(1)`. Number of decimal places used when rounding neutral
  masses before counting distinct values. Default `4`.

## Value

A `data.frame` with one row per pcgroup (those with at least one CAMERA
adduct annotation), containing:

- `pcgroup`:

  CAMERA pcgroup identifier.

- `median_rt`:

  Median RT (seconds) of all features in the pcgroup.

- `n_distinct_neutralM`:

  Number of distinct neutral masses found within the pcgroup.

- `grouping_class`:

  Ordered factor: `"1 neutral mass"`, `"2 neutral masses"`,
  `"3-5 neutral masses"`, `">5 neutral masses"`.
