# Constructor for PeakAnnotationFromLibraryMatch

Thin wrapper around
[`createPeakAnnotationFromLibraryMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createPeakAnnotationFromLibraryMatch.md),
the canonical user-facing constructor. See that function for full
parameter documentation.

## Usage

``` r
PeakAnnotationFromLibraryMatch(
  library_matches,
  FT = data.frame(),
  adducts_df,
  neutral_loss_table,
  polarity = "positive",
  mz_threshold = 0.005,
  ppm = 30,
  nl_ppm = 5,
  rt_threshold = 3,
  rt_exp_threshold = 5,
  ppm_exp = 5,
  iso_threshold = 1.003355,
  max_charge = 3L,
  max_isotope = 4L
)
```

## Arguments

- library_matches:

  A `list` of spectral library match results (e.g. output of
  `direct_match_MS1` or `direct_match_MS2`). Each element must contain:

  - `Compound_name`

  - `Compound_M`

  - `Ion_type`

  - `Lib_group`

  - `Matched_features_to_peaks`

  - `Matched_features_mz`

  Optional (used when available):

  - `Matched_features_rt`

  - `Matched_features_int`

  - `Precursor_mz`

  - `inner_cosine`, `outer_cosine`, `Match_peak_score`

- FT:

  `data.frame`. Full feature pool table (all features, not only
  library-matched ones). Required columns: `feature_name`, one of
  `mzmed`/`mz`. Optional: `rtmed`/`rt`, intensity column
  (`maxo`/`into`/`intensity`). Pass
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html) to run Phase
  1 only (no M-expansion).

- adducts_df:

  `data.frame`. Adduct definition table with columns `Ion.name`,
  `Ion.polarity`, `Ion.mass`.

- neutral_loss_table:

  `data.frame`. Neutral loss table with columns `Accurate.Mass`,
  `Formula` (or `Neutral.Loss`), `Pos`, `Neg`.

- polarity:

  `"positive"` or `"negative"` (normalised internally to `"pos"` /
  `"neg"`).

- mz_threshold:

  Numeric. Absolute m/z tolerance in Da; acts as floor for all
  OR-condition ppm gates. Default `0.005`.

- ppm:

  Numeric. ppm gate for Phase 1 adduct matching. Effective tolerance is
  `max(mz_threshold, mz * ppm * 1e-6)`. Controls confidence tier
  assignment: `<=5` Confirmed, `<=10` Probable, `<=20` Tentative, `<=30`
  Questionable, `>30` rejected. Default `30`.

- nl_ppm:

  Numeric. ppm gate for neutral-loss matching in both phases. Default
  `5`.

- rt_threshold:

  Numeric. RT tolerance (seconds) for Phase 1 within-group NL matching.
  Default `3`.

- rt_exp_threshold:

  Numeric. RT window (seconds) around the compound reference RT for
  Phase 2 pool search. Default `5`.

- ppm_exp:

  Numeric. ppm gate for Phase 2 adduct search. Default `5`.

- iso_threshold:

  Numeric. Isotope mass delta (12C/13C spacing). Default `1.003355`.

- max_charge:

  Integer. Maximum charge state for isotope detection. Default `3L`.

- max_isotope:

  Integer. Maximum isotope peaks per charge state. Default `4L`.

## Value

A `PeakAnnotationFromLibraryMatch` object.

## See also

[`createPeakAnnotationFromLibraryMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createPeakAnnotationFromLibraryMatch.md)
