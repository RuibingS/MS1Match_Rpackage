# Extract precursor m/z from a direct_match entry

Tries several common field name spellings in order.

## Usage

``` r
extract_precursor_mz(entry)
```

## Arguments

- entry:

  A named list (one element of direct_match_res_MS2)

## Value

Numeric precursor m/z, or `NA_real_` if not found
