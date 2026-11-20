# Convert a parsed library to standard NIST MSP format

One-time conversion of any in-house library (Bruker `.library` or
CSU-PMF `.mspLib`) to a clean NIST MSP file suitable for Zenodo
deposition and downstream use with
[`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md).

What this function does:

- Normalizes field names to NIST standard (`IonPolarity: pos` →
  `Ion_mode: POSITIVE`, `InstName` → `Instrument`, etc.)

- Decodes CSU-PMF `SYNON` codes to proper fields (`$:11Pos` →
  `Ion_mode: POSITIVE`, `$:05` → `RetentionTime`, etc.)

- Infers `Spectrum_type: MS1` when the `SYNON $:00` value is `in-source`
  and no `MSMS` field is present

- Calculates `ExactMass` from `Formula` via
  [`calculate_exact_mass`](https://ruibings.github.io/MS1Match_Rpackage/reference/calculate_exact_mass.md)
  if absent (common in Bruker files)

- Writes peaks one pair per line (NIST standard)

- Drops Bruker-specific internal fields (`ReserveDWord`, `AnalName`,
  `BasePeak`, `PeakWidth`, etc.)

- Skips optional identifiers (InChIKey, InChI, SMILES) if absent —
  entries do not need to match MoNA completeness

## Usage

``` r
convert_to_nist_msp(
  library,
  output_path,
  source_lab = NULL,
  source_instrument = NULL
)
```

## Arguments

- library:

  A parsed library list — output of
  [`read_library`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_library.md).

- output_path:

  Path for the output `.msp` file.

- source_lab:

  Optional contributor/lab name appended to `Comments:` (e.g.
  `"Your Lab, Your University"`).

- source_instrument:

  Optional fallback instrument string used when no instrument field is
  present in an entry.

## Value

Invisibly returns `output_path`.

## Details

**Three-step workflow:**

1.  [`read_library`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_library.md)
    — tolerates any MSP-like format

2.  `convert_to_nist_msp` — one-time normalization; run once and save
    the resulting `.msp` for all future use

3.  [`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md)
    — strict NIST MSP importer with filtering, entropy QC, and
    deduplication for downstream analysis

## See also

[`read_library`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_library.md),
[`import_msp`](https://ruibings.github.io/MS1Match_Rpackage/reference/import_msp.md),
[`calculate_exact_mass`](https://ruibings.github.io/MS1Match_Rpackage/reference/calculate_exact_mass.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ---- Bruker .library → NIST MSP ----
bruker_raw <- read_library("inhouse_library.library")
convert_to_nist_msp(bruker_raw, "inhouse_library.msp",
                    source_lab = "Your Lab, Your University")

# ---- CSU-PMF .mspLib → NIST MSP ----
csu_raw <- read_library("CSUPMF_C8_lib.mspLib")
convert_to_nist_msp(csu_raw, "CSUPMF_C8.msp",
                    source_lab = "Broeckling Lab, Colorado State University")

# ---- From now on, always load the clean .msp files via import_msp() ----
inhouse_lib <- import_msp("inhouse_library.msp",
                           spectrum_type = "MS1",
                           dedup_by      = "none")
csu_lib     <- import_msp("CSUPMF_C8.msp",
                           spectrum_type = "MS1",
                           dedup_by      = "none")
} # }
```
