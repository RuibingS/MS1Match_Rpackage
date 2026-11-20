# Run XCMS preprocessing

Run XCMS preprocessing

## Usage

``` r
runXCMSProcessing(
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

  Path to directory containing mzXML or mzML files.

- cwp:

  [`CentWaveParam`](https://rdrr.io/pkg/xcms/man/findChromPeaks-centWave.html).
  NULL = defaults for Bruker maXis HD / UHPLC / positive ESI.

- pdp:

  [`PeakDensityParam`](https://rdrr.io/pkg/xcms/man/groupChromPeaks.html).
  NULL = defaults.

- obi:

  [`ObiwarpParam`](https://rdrr.io/pkg/xcms/man/adjustRtime.html) for RT
  alignment. Pass NULL to skip (single file or pre-aligned data).

- sample_name:

  Character vector, one name per file. NULL = basenames.

- sample_group:

  Character vector of group labels, parallel to `sample_name`. Example:
  `c(rep("WT",3), rep("KO",3))`.

- BPPARAM:

  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html).
  Default `SerialParam()`. Parallelise with `SnowParam(workers = 4)`.

## Value

Named list: `x_exp` (XcmsExperiment) and `pd` (sample metadata
DataFrame).
