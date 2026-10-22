#' PeakAnnotationFromLibraryMatch S4 class
#'
#' @description
#' Container for two-phase feature-level annotation derived from spectral
#' library matching.
#'
#' \strong{Phase 1} annotates features within each library-match group
#' (tiers 1a-1d, 2).  \strong{Phase 2} extends annotation to pool features
#' absent from the library by applying adduct rules to every confirmed M from
#' Phase 1 (tiers 1a_exp, 1b_exp, 1c_exp).
#'
#' Isotope annotation (\code{Check_Iso_Charge}) runs \emph{once} on the union
#' of all library-matched features and the full pool before either phase.
#' Both phases only look up the resulting labels — isotope is never re-run.
#'
#' @slot library_matches list. Output of MS1/MS2 spectral library matching.
#' @slot FT data.frame. Full feature pool table (all features, not only
#'   matched ones).  Required columns: \code{feature_name}, one of
#'   \code{mzmed}/\code{mz}.  Optional: \code{rtmed}/\code{rt}, intensity
#'   column (\code{maxo}/\code{into}/\code{intensity}).
#'   Pass \code{data.frame()} to run Phase 1 only.
#' @slot adducts_df data.frame. Adduct definition table with columns
#'   \code{Ion.name}, \code{Ion.polarity}, \code{Ion.mass}.
#' @slot neutral_loss_table data.frame. Neutral loss table with columns
#'   \code{Accurate.Mass}, \code{Formula} (or \code{Neutral.Loss}),
#'   \code{Pos}, \code{Neg}.
#' @slot polarity character. \code{"pos"} or \code{"neg"}.
#' @slot mz_threshold numeric. Absolute m/z tolerance in Da (floor for all
#'   OR-condition ppm gates). Default 0.005.
#' @slot ppm numeric. ppm gate for Phase 1 adduct matching. Default 30.
#'   Controls confidence tier: <=5 Confirmed, <=10 Probable, <=20 Tentative,
#'   <=30 Questionable, >30 rejected.
#' @slot nl_ppm numeric. ppm gate for neutral-loss matching in both phases.
#'   Default 5.
#' @slot rt_threshold numeric. RT tolerance in seconds for Phase 1 within-group
#'   NL matching. Default 3.0.
#' @slot rt_exp_threshold numeric. RT tolerance in seconds for Phase 2 pool
#'   search (window around compound reference RT). Default 5.0.
#' @slot ppm_exp numeric. ppm gate for Phase 2 adduct search. Default 5.0.
#' @slot iso_threshold numeric. Isotope mass delta (12C/13C). Default 1.003355.
#' @slot max_charge integer. Maximum charge state for isotope detection.
#'   Default 3L.
#' @slot max_isotope integer. Maximum isotope peaks per charge state.
#'   Default 4L.
#'
#' @name PeakAnnotationFromLibraryMatch-class
#' @rdname PeakAnnotationFromLibraryMatch-class
#' @exportClass PeakAnnotationFromLibraryMatch
setClass(
  "PeakAnnotationFromLibraryMatch",
  slots = c(
    library_matches    = "list",
    FT            = "data.frame",
    adducts_df         = "data.frame",
    neutral_loss_table = "data.frame",
    polarity           = "character",
    mz_threshold       = "numeric",
    ppm                = "numeric",
    nl_ppm             = "numeric",
    rt_threshold       = "numeric",
    pi_rt_threshold    = "numeric",      # RT gate for [M+H]+/[M-H]- precursor match
    rt_exp_threshold   = "numeric",
    ppm_exp            = "numeric",
    iso_threshold      = "numeric",
    max_charge         = "integer",
    max_isotope        = "integer"
  )
)

#' Constructor for PeakAnnotationFromLibraryMatch
#'
#' @description
#' Thin wrapper around \code{\link{createPeakAnnotationFromLibraryMatch}()},
#' the canonical user-facing constructor.  See that function for full parameter
#' documentation.
#'
#' @inheritParams createPeakAnnotationFromLibraryMatch
#' @return A \code{PeakAnnotationFromLibraryMatch} object.
#' @seealso \code{\link{createPeakAnnotationFromLibraryMatch}()}
#' @export
PeakAnnotationFromLibraryMatch <- function(
    library_matches,
    FT            = data.frame(),
    adducts_df,
    neutral_loss_table,
    polarity           = "positive",
    mz_threshold       = 0.005,
    ppm                = 30,
    nl_ppm             = 5.0,
    rt_threshold       = 3.0,
    rt_exp_threshold   = 5.0,
    ppm_exp            = 5.0,
    iso_threshold      = 1.003355,
    max_charge         = 3L,
    max_isotope        = 4L
) {
  createPeakAnnotationFromLibraryMatch(
    library_matches    = library_matches,
    FT            = FT,
    adducts_df         = adducts_df,
    neutral_loss_table = neutral_loss_table,
    polarity           = polarity,
    mz_threshold       = mz_threshold,
    ppm                = ppm,
    nl_ppm             = nl_ppm,
    rt_threshold       = rt_threshold,
    rt_exp_threshold   = rt_exp_threshold,
    ppm_exp            = ppm_exp,
    iso_threshold      = iso_threshold,
    max_charge         = max_charge,
    max_isotope        = max_isotope
  )
}
