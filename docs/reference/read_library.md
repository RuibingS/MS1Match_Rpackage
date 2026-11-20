# Read any MSP-like spectral library file

Generic reader for spectral library files regardless of field naming
conventions. Accepts Bruker `.library`, CSU-PMF `.mspLib`, MoNA `.msp`,
NIST `.msp`, or any MSP-like plain text format.

This function is permissive — it reads whatever fields are present
without normalizing or filtering. Use it as **Step 1** of the workflow
before
[`convert_to_nist_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/convert_to_nist_msp.md).

For importing already-clean NIST MSP files into your analysis pipeline,
use
[`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md)
directly instead.

## Usage

``` r
read_library(
  file_path,
  dedup_by = "none",
  reserve_n = 10000L,
  max_records = -1L
)
```

## Arguments

- file_path:

  Path to the library file (`.msp`, `.library`, `.mspLib`, `.txt`).

- dedup_by:

  Deduplication strategy. Default `"none"` — keep all entries so nothing
  is lost before conversion. Options: `"name_adduct_mode_energy"`,
  `"name_adduct_mode"`, `"name_only"`, `"none"`.

- reserve_n:

  Memory hint for expected number of compounds (default 10000).

- max_records:

  Integer. Stop after this many records (-1 = no limit).

## Value

A named list of compounds, each containing raw metadata fields and a
`Peaks` data frame (`mz`, `rel_intensity`).

## See also

[`convert_to_nist_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/convert_to_nist_msp.md),
[`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Step 1a: Read Bruker .library (non-standard fields: IonPolarity, InstName, etc.)
bruker_raw <- read_library("inhouse_library.library")

# Step 1b: Read CSU-PMF .mspLib (SYNON-encoded metadata)
csu_raw <- read_library("CSUPMF_C8_lib.mspLib")

# Step 2: Convert both to clean NIST MSP (run once, then use the .msp files)
convert_to_nist_msp(bruker_raw, "inhouse_library.msp",
                    source_lab = "Your Lab, Your University")
convert_to_nist_msp(csu_raw, "CSUPMF_C8.msp",
                    source_lab = "Broeckling Lab, Colorado State University")
} # }
```
