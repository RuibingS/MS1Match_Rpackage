# Create a DirectLibMatch object

Create a DirectLibMatch object

## Usage

``` r
createDirectLibMatch(
  FT,
  compounds_list,
  rt_thr_compound,
  rt_thr_coelution = 5,
  ppm,
  ms_level = "MS1",
  polarity = "pos"
)
```

## Arguments

- FT:

  A data.frame of features

- compounds_list:

  A list of compounds

- rt_thr_compound:

  Numeric loose RT threshold for precursor gate vs library RT (e.g. 30
  s)

- rt_thr_coelution:

  Numeric strict RT threshold for co-elution pool anchored to the
  matched precursor feature (e.g. 5 s). Defaults to `5`.

- ppm:

  Numeric ppm tolerance

- ms_level:

  Either `"MS1"` or `"MS2"`

- polarity:

  Either `"pos"` or `"neg"`

## Value

A DirectLibMatch S4 object
