# Annotate features using library matching and cosine similarity

This function matches features to compounds using direct m/z and RT
comparison, followed optionally by cosine similarity scoring. Cosine
scoring is computed in parallel using BiocParallel.

Two RT thresholds are used:

- `rt_thr_compound`:

  Loose threshold for the precursor gate: how far a feature's RT may
  deviate from the library compound's RT (e.g. 30 s).

- `rt_thr_coelution`:

  Strict threshold for the co-elution pool: only features within this
  window of the matched precursor feature are considered as fragment /
  isotope candidates (e.g. 5 s). Applies to both MS1 and MS2 matching.

## Usage

``` r
libAnnotate(
  object,
  compute_cosine = TRUE,
  ppm = 10,
  BPPARAM = BiocParallel::bpparam()
)

# S4 method for class 'DirectLibMatch'
libAnnotate(
  object,
  compute_cosine = TRUE,
  ppm = 10,
  BPPARAM = BiocParallel::bpparam()
)
```

## Arguments

- object:

  A `DirectLibMatch` object containing the feature table and compound
  list

- compute_cosine:

  Logical. Should cosine similarity scores be computed? Default is
  `TRUE`

- ppm:

  Numeric. Parts-per-million (ppm) tolerance used in cosine similarity
  scoring. Default is `10`

- BPPARAM:

  A `BiocParallelParam` object controlling parallelization. Default is
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).

## Value

Updated `DirectLibMatch` object

An updated `DirectLibMatch` object with added `@results` slot containing
cosine similarity and matched features
