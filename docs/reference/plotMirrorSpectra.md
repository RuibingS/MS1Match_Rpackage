# Plot mirror plots for a list of spectral matches

Plot mirror plots for a list of spectral matches

## Usage

``` r
plotMirrorSpectra(object, matches, output_dir = NULL, ppm = 10)
```

## Arguments

- object:

  A `DirectLibMatch` object containing the compound library

- matches:

  A `list` of match results, as returned by
  [`libAnnotate()`](https://ruibings.github.io/MS1Match_Rpackage/reference/libAnnotate.md)

- output_dir:

  Optional directory to save PNG files, otherwise plots are shown
  interactively

- ppm:

  Numeric value for ppm tolerance used in mirror plot matching
