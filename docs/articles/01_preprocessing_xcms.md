# 1. Preprocessing raw files with xcms

------------------------------------------------------------------------

## Overview

This tutorial walks through the full XCMS-only pipeline from raw
mzXML/mzML files to a processed feature table. This step is optional:
start here only if you have raw files. Its feature table is the input
for vignette 2a.

The pipeline uses five step functions, or
[`runxcmsFT()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runxcmsFT.md)
to run them all:

| Step | Function | Output |
|:--:|:---|:---|
| 0 | runXCMSProcessing() | XcmsExperiment + sample metadata |
| 1 | getChromPeaks() | Raw chromPeaks data.frame |
| 2 | linkCPtoFeatures() | chromPeaks + feature_id column |
| 3 | buildFeatureTable() | Feature defs + intensity matrix |
| 4 | processFeatureTable() | Sorted, feature_name, intensity_median, intensity_mean |
| 5 | runxcmsFT() | Steps 0-4 in one call |

------------------------------------------------------------------------

## 1 Load packages

``` r

library(MS1Match)
library(xcms)

library(MsExperiment)
library(BiocParallel)
library(S4Vectors)
```

------------------------------------------------------------------------

## 2 Set XCMS parameters

Parameters are defined as standard xcms `Param` objects and passed into
the pipeline. All three are optional – if omitted, the pipeline falls
back to built-in defaults tuned for Bruker maXis HD / UHPLC / positive
ESI.

### 2.1 Peak picking - CentWave

``` r

cwp <- CentWaveParam(
  peakwidth          = c(5, 20),   
  ppm                = 5,          
  snthresh           = 10,         
  mzdiff             = 0.001,      
  prefilter          = c(3, 1000), 
  noise              = 500,        
  verboseBetaColumns = TRUE        # retains beta_cor / beta_snr in chromPeaks
)
```

> **Note:** For detailed guidance on parameter selection, please refer
> to the [xcms
> documentation](https://bioconductor.org/packages/release/bioc/html/xcms.html)
> and the
> [`?CentWaveParam`](https://rdrr.io/pkg/xcms/man/findChromPeaks-centWave.html)
> help page.

### 2.2 Feature grouping - PeakDensity

``` r

pdp <- PeakDensityParam(
  sampleGroups = c("3ul"),
  bw           = 5,     
  minFraction  = 0.5,   
  binSize      = 0.01  
)
```

### 2.3 RT alignment - Obiwarp

``` r

obi <- ObiwarpParam(binSize = 0.5)
```

------------------------------------------------------------------------

## 3 Run the full pipeline

### 3.1 Pipeline function

[`runxcmsFT()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runxcmsFT.md)
runs all five steps and returns a named list.

``` r

xcms_out <- runxcmsFT(
  files_dir    = system.file("extdata", "mix84_data", package = "MS1Match"), 
  cwp          = cwp,
  pdp          = pdp,
  obi          = obi,
  sample_name  = c("3ul"),
  sample_group = c("3ul"),
  BPPARAM      = SnowParam(workers = 4)   
)
```

Retrieve xcms_outputs:

``` r

xcms_output_result <- xcms_out$result     # list(x_exp, pd)
xcms_output_cp     <- xcms_out$cp         # all detected chromPeaks
xcms_output_cpft   <- xcms_out$cp_linked  # chromPeaks with feature_id
xcms_output_ft     <- xcms_out$ft         # final processed feature table
```

### 3.2 Step-by-step approach

Users can inspect or modify intermediate objects.

``` r

# Step 0 -- XCMS processing
result <- runXCMSProcessing(
  files_dir    = "yourdata/mzXML/",
  cwp          = cwp,
  pdp          = pdp,
  obi          = obi,
  sample_name  = c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"),
  sample_group = c(rep("WT", 3), rep("KO", 3)),
  BPPARAM      = SnowParam(workers = 4)
)

# Step 1 -- raw chromPeaks
cp <- getChromPeaks(result$x_exp)

# Step 2 -- chromPeaks linked to feature IDs
cp_linked <- linkCPtoFeatures(result$x_exp)

# Step 3 -- feature table (feature defs + intensity matrix)
ft_raw <- buildFeatureTable(result$x_exp)

# Step 4 -- post-process: sort by RT, add feature_name + intensity_median, intensity_mean
xcms_output_ft <- processFeatureTable(ft_raw, result$x_exp)
```

------------------------------------------------------------------------

## 4 Inspect xcms_outputs

### 4.1 chromPeaks (Step 1)

One row per detected peak across all samples.

``` r

dim(cp)
head(cp)
```

Key columns: `mz`, `mzmin`, `mzmax`, `rt`, `rtmin`, `rtmax`, `into`
(integrated area), `maxo` (apex intensity), `sample`, `sample_name`.

### 4.2 chromPeaks linked to features (Step 2)

Same as above with `feature_id` appended. `NA` means the peak was not
grouped into any feature.

``` r

# How many peaks were assigned to a feature?
table(is.na(cp_linked$feature_id))

# Peaks for a specific feature
cp_linked[cp_linked$feature_id == "FT0001" & !is.na(cp_linked$feature_id), ]
```

### 4.3 Feature table (Step 4)

One row per feature. Columns: feature metadata, one intensity column per
sample, and `intensity_median` and `intensity_mean`.

``` r

dim(xcms_output_ft)
head(xcms_output_ft[, 1:10])


summary(xcms_output_ft$intensity_median)
```

------------------------------------------------------------------------

## 5 Export

``` r

exportFeatureTable(xcms_output_ft, "xcms_output_feature_table.csv")
```

------------------------------------------------------------------------

## 6 Parameter defaults reference

If no parameters are passed,
[`runXCMSProcessing()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runXCMSProcessing.md)
/
[`runxcmsFT()`](https://ruibings.github.io/MS1Match_Rpackage/reference/runxcmsFT.md)
use:

| Parameter     | Default       | Notes                                     |
|:--------------|:--------------|:------------------------------------------|
| peakwidth     | c(5, 20)      | UHPLC; use c(10, 30) for standard HPLC    |
| ppm           | 5             | maXis HD; use 10 for generic QTOF         |
| snthresh      | 10            | Lower to 5 if standards are missed        |
| mzdiff        | 0.001         | High-res; use 0.01 for low-res            |
| noise         | 500           | Set from diagnoseParams()                 |
| bw            | 5             | Increase to 10 for larger RT drift        |
| minFraction   | 0.5           | Use 0.25 for exploratory runs             |
| binSize (pdp) | 0.01          | High-res; use 0.1 for low-res             |
| BPPARAM       | SerialParam() | Use SnowParam(workers = 4) to parallelise |

------------------------------------------------------------------------

## Session info

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26100)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.7.3        svglite_2.2.2      cli_3.6.6          knitr_1.52        
    ##  [5] rlang_1.3.0        xfun_0.62          stringi_1.8.9      otel_0.2.0        
    ##  [9] textshaping_1.0.5  jsonlite_2.0.0     glue_1.8.1         htmltools_0.5.9   
    ## [13] ragg_1.5.2         sass_0.4.10        scales_1.4.0       rmarkdown_2.32    
    ## [17] evaluate_1.0.5     jquerylib_0.1.4    kableExtra_1.4.1   fastmap_1.2.0     
    ## [21] yaml_2.3.12        lifecycle_1.0.5    stringr_1.6.0      compiler_4.6.1    
    ## [25] RColorBrewer_1.1-3 fs_2.1.0           htmlwidgets_1.6.4  rstudioapi_0.19.0 
    ## [29] farver_2.1.2       systemfonts_1.3.2  digest_0.6.39      viridisLite_0.4.3 
    ## [33] R6_2.6.1           magrittr_2.0.5     bslib_0.12.0       tools_4.6.1       
    ## [37] xml2_1.5.2         pkgdown_2.2.1      cachem_1.1.0       desc_1.4.3
