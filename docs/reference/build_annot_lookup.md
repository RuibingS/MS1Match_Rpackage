# Build a per-feature annotation lookup table

Indexes an annotate_res data frame by `feature_name`. Only rows with
`annotation_phase == "Phase1"` are used (Phase2 `_exp` rows are dropped
here; call this on the full table and it will filter internally).
`adduct_match` is taken directly from `grp_anno_match` with no
group-level fallback: in this data model, `grp_anno_match` is already
populated whenever a feature genuinely has a Tier1a/Tier1b precursor or
isotope match, so a feature with an empty `grp_anno_match` should get no
label — not one inherited from the compound's overall ion type. (A prior
version promoted `group_ion_type` into `adduct_match` for any row
lacking its own annotation; that caused every unannotated fragment in a
group to share the same "M+H+" label text as the true precursor feature,
so the plot's per-label intensity-based dedup would silently keep
whichever row was tallest — not necessarily the real precursor — and
point the leader line at the wrong peak.) A `nl_label` field is also
stored: it is non-empty only when `grp_anno_match` itself already
encodes a neutral loss (e.g. `"[M+H-H2O]+ 274.131742"` -\>
`"[M+H-H2O]+"`). `NL_match` is not used for plot labels.

## Usage

``` r
build_annot_lookup(annotate_res, group_ion_type = "")
```

## Arguments

- annotate_res:

  Data frame of per-feature annotation rows (e.g.
  `StM14_annotate_res_ms1`). Expected columns include `feature_name`,
  `grp_anno_match`, `grp_anno_ppm`, `grp_anno_mda`,
  `grp_anno_confidence`, `iso_annotation`, `annotation_phase`.

- group_ion_type:

  Unused (kept for call-site compatibility).

## Value

Named list keyed by `feature_name`
