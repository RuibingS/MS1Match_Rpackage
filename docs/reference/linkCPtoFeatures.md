# Link chromPeaks to feature IDs

Link chromPeaks to feature IDs

## Usage

``` r
linkCPtoFeatures(x_exp)
```

## Arguments

- x_exp:

  XcmsExperiment.

## Value

data.frame: all chromPeak columns + `feature_name` + `sample_name`.
Unassigned peaks have `NA` in `feature_name`.
