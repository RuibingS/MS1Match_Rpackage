# Build simple feature table

@description Feature definitions cbind'd with the per-sample `into`
intensity matrix.

## Usage

``` r
buildFeatureTable(x_exp)
```

## Arguments

- x_exp:

  XcmsExperiment after fillChromPeaks.

## Value

data.frame, one row per feature, rownames = feature IDs.
