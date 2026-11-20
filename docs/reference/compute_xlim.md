# Compute a shared x-axis range for the MS1 + MS2 mirror plot pair

Compute a shared x-axis range for the MS1 + MS2 mirror plot pair

## Usage

``` r
compute_xlim(
  ref_mz_ms1_full,
  ref_mz_ms2_full,
  precursor_mz = NULL,
  padding_frac = 0.05,
  precursor_extra_frac = 0.04
)
```

## Arguments

- ref_mz_ms1_full:

  Numeric vector; all MS1 library peak m/z

- ref_mz_ms2_full:

  Numeric vector; all MS2 library peak m/z

- precursor_mz:

  Numeric scalar or `NULL`

- padding_frac:

  Fractional padding added to each side

- precursor_extra_frac:

  Extra right-side padding for precursor marker

## Value

Length-2 numeric vector `c(xlo, xhi)`
