# Export a compound library list to NIST MSP format

Writes a named list of compounds (as returned by
[`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md))
to a NIST-standard `.msp` file. Each compound is written as a
self-contained block with metadata fields followed by its peak list.

Use this function for:

- Saving a filtered or deduplicated library back to disk

- Exporting processed libraries for use in other tools (e.g. NIST MS
  Search, MZmine, MS-DIAL)

- Archiving subsets produced by
  [`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md)

## Usage

``` r
export_msp(
  library,
  file_path = "library_export.msp",
  fields = c("Name", "PrecursorMZ", "Adduct", "Formula", "ExactMass", "InChIKey",
    "InChI", "SMILES", "IonMode", "CollisionEnergy", "InstrumentType", "Instrument",
    "IonizationType", "Spectrum_type", "Comments"),
  mz_digits = 4L,
  intensity_digits = 0L
)
```

## Arguments

- library:

  A named list of compounds as returned by
  [`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md).
  Each element must contain a `Peaks` data frame with columns `mz` and
  `rel_intensity`, plus any metadata fields to be written.

- file_path:

  Path to the output `.msp` file. Parent directory must exist. Existing
  files are overwritten.

- fields:

  Character vector of metadata field names to write for each compound,
  in the order they will appear in the output. Fields absent from a
  compound are silently skipped. Defaults to all standard NIST MSP
  fields:
  `c("Name", "PrecursorMZ", "Adduct", "Formula", "ExactMass", "InChIKey", "InChI", "SMILES", "IonMode", "CollisionEnergy", "InstrumentType", "Instrument", "IonizationType", "Spectrum_type", "Comments")`.

- mz_digits:

  Integer. Number of decimal places for `mz` values in the peak list.
  Default `4`.

- intensity_digits:

  Integer. Number of decimal places for intensity values in the peak
  list. Default `0` (integer-like output, consistent with NIST
  convention).

## Value

Invisibly returns `file_path`. Called primarily for its side effect of
writing the `.msp` file.

## See also

[`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md),
[`convert_to_nist_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/convert_to_nist_msp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Round-trip: import, filter, export
CBIO_library <- import_msp(
  system.file("extdata", "CBIO", "CBIO_Standards.msp", package = "MS1Match"),
  ionPolarity   = "pos",
  spectrum_type = "MS1",
  filter        = list(field = "Comments", pattern = "MSMLS")
)

export_msp(
  library   = CBIO_library,
  file_path = "CBIO_MSMLS_filtered.msp"
)

# Custom field order and precision
export_msp(
  library          = CBIO_library,
  file_path        = "CBIO_MSMLS_filtered.msp",
  fields           = c("Name", "Formula", "ExactMass", "Adduct",
                       "IonMode", "Comments"),
  mz_digits        = 6,
  intensity_digits = 2
)
} # }
```
