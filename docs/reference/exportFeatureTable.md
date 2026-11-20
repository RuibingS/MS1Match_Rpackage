# Export feature table to CSV

Writes the processed feature table to a CSV file. The `peakidx` list
column is collapsed to a single semicolon-separated string per row so
the file remains rectangular and readable in Excel / any CSV reader.

## Usage

``` r
exportFeatureTable(ft, file, sep = ",", ...)
```

## Arguments

- ft:

  data.frame from
  [`processFeatureTable()`](https://ruibings.github.io/MS1Match_Rpackage/reference/processFeatureTable.md)
  or `runxcmsFT()$ft`.

- file:

  Output file path, e.g. `"results/feature_table.csv"`. The parent
  directory is created automatically if it does not exist.

- sep:

  Field separator. Default `","`.

- ...:

  Additional arguments passed to
  [`write.csv()`](https://rdrr.io/r/utils/write.table.html).

## Value

`ft` invisibly (so the call can be piped).

## Details

Any other residual list columns (unexpected) are similarly collapsed as
a safety measure.

## Examples

``` r
if (FALSE) { # \dontrun{
out <- runxcmsFT("data/mzXML/", sample_name = c("S1","S2","S3"))
exportFeatureTable(out$ft, "results/feature_table.csv")
} # }
```
