# Normalize a target list data frame to a standard schema

Accepts target lists from different sources and maps heterogeneous
column names to a unified schema with columns: `Name`, `Formula`,
`exactMass`, `[M+H]+`, `[M-H]-`, and optionally `RT`.

Supported input column name variants:

- Name:

  CNAME, Name

- Formula:

  FORMULA, Formula, Formula (as seen in ...), any column starting with
  "Formula"

- exactMass:

  exactMass, MONO_MASS, Neutral mass, Neutral.mass

- M+H+:

  M+H+, M+H, MH, MplusH

- M-H-:

  M-H-, M-H, MH_neg, MminusH

- RT:

  RT, RetentionTime, rt

If `[M+H]+` or `[M-H]-` are absent but `exactMass` is present, the ion
m/z values are computed automatically.

## Usage

``` r
normalize_target_list(target_list)
```

## Arguments

- target_list:

  A `data.frame` or `tibble` target list

## Value

A `data.frame` with standardized column names
