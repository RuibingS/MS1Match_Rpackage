# Pipeline function feature table pipeline

Pipeline function feature table pipeline

## Usage

``` r
runxcmsFT(
  files_dir,
  cwp = NULL,
  pdp = NULL,
  obi = ObiwarpParam(),
  sample_name = NULL,
  sample_group = NULL,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- files_dir:

  Passed to
  [`runXCMSProcessing()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runXCMSProcessing.md).

- cwp, pdp, obi, sample_name, sample_group, BPPARAM:

  Passed to
  [`runXCMSProcessing()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runXCMSProcessing.md).

## Value

Named list:

- result:

  Output of
  [`runXCMSProcessing()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runXCMSProcessing.md)
  (`x_exp`, `pd`).

- cp:

  Raw chromPeaks data.frame.

- cp_linked:

  chromPeaks with `feature_name` column.

- ft:

  Processed feature table (sorted, `intensity_median`).

## Examples

``` r
if (FALSE) { # \dontrun{
out      <- runxcmsFT("data/mzXML/",
                        sample_name  = c("WT_1","WT_2","WT_3","KO_1","KO_2","KO_3"),
                        sample_group = c(rep("WT",3), rep("KO",3)))
out_ft  <- out$ft
} # }
```
