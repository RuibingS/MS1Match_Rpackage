# Class for Direct Target Matching

Holds all inputs and results for matching a feature table against a
compound target list using primary ion m/z (`"[M+H]+"` or `"[M-H]-"`).

The target list (`Comp_data`) can be in any supported format — HMDB,
plate-based, PA14, etc. Column normalisation is handled automatically
inside
[`PImatch_target()`](https://ruibings.github.io/MS1Match_Rpackage/reference/PImatch_target.md).

## Slots

- `FT`:

  A `data.frame` of detected features (needs `feature_name` and
  `mzmed`/`mz` columns; optionally `rtmed`/`rt` for RT gating).

- `Comp_data`:

  A `data.frame` target list in any supported format (HMDB, plate-based,
  PA14, …).

- `PIon`:

  Character. Primary ion selector: `"[M+H]+"` for positive mode,
  `"[M-H]-"` for negative mode.

- `diff_mz_thr`:

  Numeric. Absolute m/z tolerance in Da (floor applied when ppm
  tolerance is smaller, e.g. at low m/z).

- `ppm`:

  Numeric. Parts-per-million m/z tolerance.

- `rt_thr_compound`:

  Numeric. RT window in seconds for feature vs compound RT gating. Only
  applied when the target list contains an RT column and the feature
  table has RT values. Set to `Inf` to disable.

- `results`:

  List of raw match results (populated by
  [`PImatch_target()`](https://ruibings.github.io/MS1Match_Rpackage/reference/PImatch_target.md)).
