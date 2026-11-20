# Import a standardized NIST MSP spectral library

Full-featured importer for clean NIST-format `.msp` files. Supports
filtering by spectrum type, ion polarity, ionization method, collision
energy and instrument name. Applies smart deduplication — all
implemented in C++ for performance.

Use this function for:

- MoNA `.msp` downloads

- Files previously converted with
  [`convert_to_nist_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/convert_to_nist_msp.md)

- Any NIST-standard `.msp` file

For raw non-standard files (Bruker `.library`, CSU-PMF `.mspLib`), use
[`read_library`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_library.md)
→
[`convert_to_nist_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/convert_to_nist_msp.md)
first.

## Usage

``` r
import_msp(
  file_path,
  spectrum_type = NULL,
  ionPolarity = NULL,
  ioniMethod = NULL,
  colEnergy = NULL,
  instName = NULL,
  filter = NULL,
  dedup_by = c("name_adduct_mode_energy", "name_adduct_mode", "name_only", "none"),
  require_formula = FALSE,
  reserve_n = 60000L,
  max_records = -1L
)
```

## Arguments

- file_path:

  Path to the NIST `.msp` file.

- spectrum_type:

  Optional filter: `"MS1"` or `"MS2"`. MoNA entries labelled
  `Spectrum_type: MS` match `"MS1"`.

- ionPolarity:

  Optional filter: `"POSITIVE"` or `"NEGATIVE"` (also accepts `"pos"`,
  `"neg"`, `"P"`, `"N"`).

- ioniMethod:

  Optional filter for ionization method (string match).

- colEnergy:

  Optional filter for collision energy (string match).

- instName:

  Optional filter for instrument name (string match).

- filter:

  Optional filter for metadata fields. Accepts either:

  - A **single** named list with elements `field` and `pattern` — drops
    entries whose field does not contain the pattern:

        filter = list(field   = "Instrument_type",
                         pattern = "LC-ESI-QTOF")

  - A **list of** such named lists for multiple conditions (AND logic —
    all conditions must be satisfied):

        filter = list(
              list(field = "Instrument_type", pattern = "LC-ESI-QTOF"),
              list(field = "Instrument",      pattern = "Agilent 6520 QTOF")
            )

  Match is case-sensitive and fixed.

- dedup_by:

  Deduplication grouping key:

  `"name_adduct_mode_energy"`

  :   (default) One entry per Name + Adduct + IonMode + CollisionEnergy.
      Preserves fragmentation diversity across energies — recommended
      for MoNA MS2.

  `"name_adduct_mode"`

  :   One entry per Name + Adduct + IonMode.

  `"name_only"`

  :   One entry per compound name — recommended for CSU-PMF and other
      MS1 libraries where collision energy variants are not meaningful.
      Keeps the entry with the most peaks.

  `"none"`

  :   No deduplication — keep all entries. Recommended for small
      in-house libraries.

- require_formula:

  Logical. If `TRUE`, entries with a missing or invalid formula
  (`"unk"`, empty) are discarded after parsing. Default `FALSE`. Set to
  `TRUE` when exact mass is required for downstream matching or
  annotation — for example, CSU-PMF MS1 libraries where entries without
  a formula cannot be used.

- reserve_n:

  Memory hint: expected number of compounds (default 60000). Increase
  for very large files to reduce internal reallocations.

- max_records:

  Integer. Stop after this many records (-1 = no limit).

## Value

A named list of compounds. Each element contains standardized metadata
fields and a `Peaks` data frame (`mz`, `rel_intensity`).

## See also

[`read_library`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_library.md),
[`convert_to_nist_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/convert_to_nist_msp.md),
[`calculate_exact_mass`](https://ruibings.github.io/MS1Match_Rpackage/reference/calculate_exact_mass.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# MoNA MS1 QTOF library — single instrument type filter
mona_lib <- import_msp(
  system.file("extdata", "Massbank_MS1",
              "massbank_result.nist_msp_MS1_TOF.msp", package = "MS1Match"),
  spectrum_type = "MS1",
  ionPolarity   = "POSITIVE",
  dedup_by      = "name_only",
  filter        = list(field = "Instrument_type", pattern = "LC-ESI-QTOF")
)

# Filter by both instrument type AND specific instrument model
# Keeps only Agilent 6520 QTOF entries
mona_agilent <- import_msp(
  "ms1_tof_positive_filtered.msp",
  spectrum_type = "MS1",
  ionPolarity   = "pos",
  filter = list(
    list(field = "Instrument_type", pattern = "LC-ESI-QTOF"),
    list(field = "Instrument",      pattern = "Agilent 6520 QTOF")
  )
)

# Filter by instrument type only (keeps all LC-ESI-QTOF regardless of model)
mona_qtof <- import_msp(
  "ms1_tof_positive_filtered.msp",
  spectrum_type = "MS1",
  ionPolarity   = "pos",
  filter        = list(field = "Instrument_type", pattern = "LC-ESI-QTOF")
)

# No filter — import everything in the pre-filtered MSP
mona_all <- import_msp(
  "ms1_tof_positive_filtered.msp",
  spectrum_type = "MS1",
  ionPolarity   = "pos"
)

# In-house MS1 library (after conversion from Bruker .library)
inhouse_lib <- import_msp(
  system.file("extdata", "CBIO", "CBIO_Standards.msp", package = "MS1Match"),
  spectrum_type = "MS1",
  dedup_by      = "name_only"
)

# CSU-PMF MS1 library — one entry per compound, formula required
csu_lib <- import_msp(
  system.file("extdata", "CSUPMF_C8_lib", "CSUPMF_library.msp",
              package = "MS1Match"),
  spectrum_type   = "MS1",
  dedup_by        = "name_only",
  require_formula = TRUE
)

# Access individual compounds
inhouse_lib[["Caffeine"]]$ExactMass
inhouse_lib[["Caffeine"]]$Peaks
} # }
```
