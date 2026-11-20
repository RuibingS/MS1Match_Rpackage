# Plot a mirror plot between a query feature and its reference compound

Plot a mirror plot between a query feature and its reference compound

## Usage

``` r
plotMirrorSpectraSingle(match, object, output_dir = NULL, ppm = 10)
```

## Arguments

- match:

  A single entry (list) from a DirectLibMatch result.

- object:

  The DirectLibMatch object containing @compounds_list.

- output_dir:

  Optional directory to save the mirror plot as PNG.

- ppm:

  m/z tolerance (default 10 ppm).

## Value

Invisibly returns the ggplot object.
