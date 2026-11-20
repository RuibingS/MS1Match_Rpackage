# Class for Direct Library Matching

Represents a configuration for matching detected features against a
compound library, using m/z and RT thresholds, ppm tolerance, MS level,
and ion polarity.

## Slots

- `FT`:

  A `data.frame` of detected features (must include m/z and RT columns)

- `compounds_list`:

  A `list` of known compounds with metadata

- `rt_thr_compound`:

  Numeric loose RT threshold for precursor gate vs library RT (e.g. 30
  s)

- `rt_thr_coelution`:

  Numeric strict RT threshold for co-elution features to the matched
  precursor feature (e.g. 5 s).

- `ppm`:

  Numeric parts-per-million tolerance for m/z

- `ms_level`:

  Character, either `"MS1"` or `"MS2"`

- `polarity`:

  Character, either `"pos"` or `"neg"`

- `results`:

  List of matched results
