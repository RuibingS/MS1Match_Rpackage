# Get MS2 peaks with feature_name assigned

Wrapper around `chromPeakSpectra()` that maps `chrom_peak_id` to
`feature_name` and adds it as a spectra variable.

## Usage

``` r
get_ms2_peaks(x_exp, expand_rt = 10, expand_mz = 0.02, msms_files = NULL, ...)
```

## Arguments

- x_exp:

  A processed `XcmsExperiment` object.

- expand_rt:

  Numeric. RT expansion in seconds. Default `10`.

- expand_mz:

  Numeric. m/z expansion in Da. Default `0.02`.

- ...:

  Additional arguments passed to `chromPeakSpectra()`.

## Value

A `Spectra` object with `feature_name` added.
