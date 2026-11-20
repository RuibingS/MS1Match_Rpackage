# Annotate features from spectral library match results

Two-phase annotation dispatched on a
[`PeakAnnotationFromLibraryMatch-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationFromLibraryMatch-class.md)
object.

**Phase 1** — `annotate_librarymatch_phase1_cpp()`:

- Runs `Check_Iso_Charge` *once* on the union of all library-matched
  features and the full pool. Produces a global isotope label table
  shared by both phases — isotope is never re-computed.

- For each library-match group: adduct match (Tier 1a), isotope lookup
  (Tier 1b), pairwise NL within group (Tier 1c), precursor/MS2 match
  (Tier 1d), unexplained (Tier 2).

- RT is used automatically when available; if RT values are `NA` the RT
  gate simply does not fire.

- Returns annotated rows plus `confirmed_M` (confirmed neutral masses
  with reference RT) and `iso_labels` (global isotope table) for Phase
  2.

**Phase 2** — `annotate_librarymatch_phase2_cpp()` (only when `FT` slot
is non-empty):

- For each confirmed M from Phase 1, scans pool features not already
  annotated within `rt_exp_threshold` s of the compound reference RT and
  `ppm_exp` ppm of any adduct mass (Tier 1a_exp).

- Looks up iso labels from the Phase 1 global table (Tier 1b_exp).

- Runs pairwise NL within the same-M expansion group (Tier 1c_exp).

## Usage

``` r
# S4 method for class 'PeakAnnotationFromLibraryMatch'
annotateLibMatch(object, cutoff_MH = FALSE, ...)
```

## Arguments

- object:

  A
  [`PeakAnnotationFromLibraryMatch-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationFromLibraryMatch-class.md)
  object.

- cutoff_MH:

  logical. If `TRUE`, adds a `cutoff_to_MH` column to the result that
  classifies each `Tier2_Unexplained` feature as `"above_MH"` or
  `"below_MH"` relative to the theoretical \\\[M+H\]^+\\ / \\\[M-H\]^-\\
  ion of the matched compound. All other tiers and rows with unknown
  `Compound_M` receive `NA`. Default `FALSE` — the column is not added
  unless explicitly requested.

- ...:

  Currently unused.

## Value

A `data.frame` with columns: `feature_name`, `mz`, `rt`, `intensity`,
`Compound_name`, `Lib_group`, `Ion_type`, `Compound_M`, `PI_match`,
`grp_anno_match`, `grp_anno_ppm`, `grp_anno_mda`, `grp_anno_confidence`,
`exp_anno_match`, `exp_anno_ppm`, `exp_anno_mda`, `exp_anno_confidence`,
`NL_match`, `iso_annotation`, `reference_tier`, `annotation_phase`,
`Match_peak_score`, `inner_cosine`, `outer_cosine`. When
`cutoff_MH = TRUE`, the additional column `cutoff_to_MH` is appended
(see `cutoff_MH` parameter above).

## See also

[`createPeakAnnotationFromLibraryMatch()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createPeakAnnotationFromLibraryMatch.md)
