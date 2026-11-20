# Create a PeakAnnotationFromLibraryMatch object

User-facing constructor for
[`PeakAnnotationFromLibraryMatch-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationFromLibraryMatch-class.md).
Bundles spectral library match results, the full feature pool, and all
annotation parameters into a single object ready for `annotate()`.

Annotation runs in two phases:

- **Phase 1** — within each library-match group: adduct match (Tier 1a),
  isotope lookup (Tier 1b), pairwise neutral loss (Tier 1c),
  precursor/MS2 match (Tier 1d), unexplained (Tier 2). Isotope labels
  are computed *once* on the union of all library features and the full
  pool before any per-group work.

- **Phase 2** — M-expansion: for each confirmed neutral mass M from
  Phase 1, pool features not yet annotated are scanned within
  `rt_exp_threshold` seconds and `ppm_exp` ppm (Tier 1a_exp, 1b_exp,
  1c_exp). Skipped when `FT` is not supplied.

RT is used automatically when present; if RT values are `NA` the RT gate
does not fire — no separate with-RT / without-RT paths needed.

## Usage

``` r
createPeakAnnotationFromLibraryMatch(
  library_matches,
  FT = data.frame(),
  adducts_df,
  neutral_loss_table,
  polarity = "positive",
  mz_threshold = 0.005,
  ppm = 30,
  nl_ppm = 5,
  rt_threshold = 3,
  pi_rt_threshold = 30,
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

- pi_rt_threshold:

  Numeric. RT tolerance (seconds) for matching the precursor
  `[M+H]+`/`[M-H]-` feature against the library compound RT. Features
  outside this window are not considered as precursor matches. Default
  `30`.

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

A
[`PeakAnnotationFromLibraryMatch-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationFromLibraryMatch-class.md)
object. Pass to `annotate()` to run annotation.

## See also

`annotate()`,
[`PeakAnnotationFromLibraryMatch-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationFromLibraryMatch-class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- createPeakAnnotationFromLibraryMatch(
  library_matches    = my_library_matches,
  FT            = my_feature_table,
  adducts_df         = adduct_table,
  neutral_loss_table = nl_table,
  polarity           = "positive"
)
result <- annotate(obj)
} # }
```
