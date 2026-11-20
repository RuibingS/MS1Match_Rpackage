# Build the label string for a single query peak

Build the label string for a single query peak

## Usage

``` r
build_query_label(mz, ion_type, adduct_match, mz_digits, mz_label_mode)
```

## Arguments

- mz:

  Numeric scalar

- ion_type:

  Character scalar

- adduct_match:

  Character scalar, possibly semicolon-separated

- mz_digits:

  Integer

- mz_label_mode:

  One of `"all"`, `"primary"`, `"adduct"`

## Value

Character scalar label, or `""`
