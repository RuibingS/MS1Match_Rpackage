# PeakAnnotationFromLibraryMatch S4 class

Container for two-phase feature-level annotation derived from spectral
library matching.

**Phase 1** annotates features within each library-match group (tiers
1a-1d, 2). **Phase 2** extends annotation to pool features absent from
the library by applying adduct rules to every confirmed M from Phase 1
(tiers 1a_exp, 1b_exp, 1c_exp).

Isotope annotation (`Check_Iso_Charge`) runs *once* on the union of all
library-matched features and the full pool before either phase. Both
phases only look up the resulting labels — isotope is never re-run.

## Slots

- `library_matches`:

  list. Output of MS1/MS2 spectral library matching.

- `FT`:

  data.frame. Full feature pool table (all features, not only matched
  ones). Required columns: `feature_name`, one of `mzmed`/`mz`.
  Optional: `rtmed`/`rt`, intensity column (`maxo`/`into`/`intensity`).
  Pass [`data.frame()`](https://rdrr.io/r/base/data.frame.html) to run
  Phase 1 only.

- `adducts_df`:

  data.frame. Adduct definition table with columns `Ion.name`,
  `Ion.polarity`, `Ion.mass`.

- `neutral_loss_table`:

  data.frame. Neutral loss table with columns `Accurate.Mass`, `Formula`
  (or `Neutral.Loss`), `Pos`, `Neg`.

- `polarity`:

  character. `"pos"` or `"neg"`.

- `mz_threshold`:

  numeric. Absolute m/z tolerance in Da (floor for all OR-condition ppm
  gates). Default 0.005.

- `ppm`:

  numeric. ppm gate for Phase 1 adduct matching. Default 30. Controls
  confidence tier: \<=5 Confirmed, \<=10 Probable, \<=20 Tentative,
  \<=30 Questionable, \>30 rejected.

- `nl_ppm`:

  numeric. ppm gate for neutral-loss matching in both phases. Default 5.

- `rt_threshold`:

  numeric. RT tolerance in seconds for Phase 1 within-group NL matching.
  Default 3.0.

- `rt_exp_threshold`:

  numeric. RT tolerance in seconds for Phase 2 pool search (window
  around compound reference RT). Default 5.0.

- `ppm_exp`:

  numeric. ppm gate for Phase 2 adduct search. Default 5.0.

- `iso_threshold`:

  numeric. Isotope mass delta (12C/13C). Default 1.003355.

- `max_charge`:

  integer. Maximum charge state for isotope detection. Default 3L.

- `max_isotope`:

  integer. Maximum isotope peaks per charge state. Default 4L.
