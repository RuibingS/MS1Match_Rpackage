# Read and Filter an MZmine Full Feature Table

## Usage

``` r
read_mzmine_feature_table(
  file,
  min_detect_fraction = 0.5,
  area_col = c("area", "height"),
  clean_names = TRUE,
  keep_meta = TRUE,
  zero_as_na = FALSE,
  verbose = FALSE
)
```

## Arguments

- file:

  `character(1)`. Path to the MZmine full feature table CSV (e.g.
  `"test_1_full_feature_table.csv"`).

- min_detect_fraction:

  `numeric(1)`. Minimum fraction of samples (0–1) in which a feature
  must have a non-zero area to be retained. Defaults to `0.5` (50 \\

  area_col`character(1)`. Which per-sample intensity metric to extract.
  One of `"area"` (default) or `"height"`.

  clean_names`logical(1)`. If `TRUE` (default), strip the `"datafile:"`
  prefix and `".mzXML:area"` (or `":height"`) suffix from column names,
  leaving only the bare sample names.

  keep_meta`logical(1)`. If `TRUE` (default), prepend the consensus
  columns `id`, `mz`, and `rt` to the output. Set to `FALSE` to return a
  pure sample-by-feature matrix.

  zero_as_na`logical(1)`. If `TRUE`, replace all zero area values with
  `NA` before returning. Useful for log-transformation workflows.
  Defaults to `FALSE`.

  verbose`logical(1)`. If `TRUE` (default), print a short summary of
  features before and after filtering, plus the number of samples
  detected.

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
one row per retained feature. If `keep_meta = TRUE` the first three
columns are `id`, `mz`, and `rt`, followed by one column per sample. If
`keep_meta = FALSE` only the sample columns are returned (with row names
set to the feature `id` values, if the result is coerced to a
`data.frame`). Imports an MZmine 3 `*_feature_table.csv` file, extracts
per-sample peak-area columns, optionally renames them to bare sample
names, and filters out features that were not detected in a minimum
fraction of samples. MZmine 3 names every per-sample column with the
pattern `datafile:<sample_filename>.<ext>:<metric>`. This function uses
a regex on the column suffix (`:area` or `:height`) to identify sample
columns, so it is robust to mixed file extensions (`.mzML`, `.mzXML`,
`.raw`, etc.).Features with `feature_state == "ESTIMATED"` (gap-filled)
are treated identically to `"DETECTED"` features for filtering purposes;
only cells containing `NA` or `0` count as "not detected".
