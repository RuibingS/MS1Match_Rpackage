#' Create a PeakAnnotationFromLibraryMatch object
#'
#' @description
#' User-facing constructor for
#' \code{\link{PeakAnnotationFromLibraryMatch-class}}.
#' Bundles spectral library match results, the full feature pool, and all
#' annotation parameters into a single object ready for \code{\link{annotate}()}.
#'
#' Annotation runs in two phases:
#' \itemize{
#'   \item \strong{Phase 1} — within each library-match group: adduct match
#'     (Tier 1a), isotope lookup (Tier 1b), pairwise neutral loss (Tier 1c),
#'     precursor/MS2 match (Tier 1d), unexplained (Tier 2).  Isotope labels
#'     are computed \emph{once} on the union of all library features and the
#'     full pool before any per-group work.
#'   \item \strong{Phase 2} — M-expansion: for each confirmed neutral mass M
#'     from Phase 1, pool features not yet annotated are scanned within
#'     \code{rt_exp_threshold} seconds and \code{ppm_exp} ppm (Tier 1a_exp,
#'     1b_exp, 1c_exp).  Skipped when \code{FT} is not supplied.
#' }
#'
#' RT is used automatically when present; if RT values are \code{NA} the RT
#' gate does not fire — no separate with-RT / without-RT paths needed.
#'
#' @param library_matches A \code{list} of spectral library match results
#'   (e.g. output of \code{direct_match_MS1} or \code{direct_match_MS2}).
#'   Each element must contain:
#'   \itemize{
#'     \item \code{Compound_name}
#'     \item \code{Compound_M}
#'     \item \code{Ion_type}
#'     \item \code{Lib_group}
#'     \item \code{Matched_features_to_peaks}
#'     \item \code{Matched_features_mz}
#'   }
#'   Optional (used when available):
#'   \itemize{
#'     \item \code{Matched_features_rt}
#'     \item \code{Matched_features_int}
#'     \item \code{Precursor_mz}
#'     \item \code{inner_cosine}, \code{outer_cosine}, \code{Match_peak_score}
#'   }
#' @param FT \code{data.frame}. Full feature pool table (all features,
#'   not only library-matched ones).  Required columns: \code{feature_name},
#'   one of \code{mzmed}/\code{mz}.  Optional: \code{rtmed}/\code{rt},
#'   intensity column (\code{maxo}/\code{into}/\code{intensity}).
#'   Pass \code{data.frame()} to run Phase 1 only (no M-expansion).
#' @param adducts_df \code{data.frame}. Adduct definition table with columns
#'   \code{Ion.name}, \code{Ion.polarity}, \code{Ion.mass}.
#' @param neutral_loss_table \code{data.frame}. Neutral loss table with columns
#'   \code{Accurate.Mass}, \code{Formula} (or \code{Neutral.Loss}),
#'   \code{Pos}, \code{Neg}.
#' @param polarity \code{"positive"} or \code{"negative"} (normalised
#'   internally to \code{"pos"} / \code{"neg"}).
#' @param mz_threshold Numeric. Absolute m/z tolerance in Da; acts as floor
#'   for all OR-condition ppm gates. Default \code{0.005}.
#' @param ppm Numeric. ppm gate for Phase 1 adduct matching.  Effective
#'   tolerance is \code{max(mz_threshold, mz * ppm * 1e-6)}.
#'   Controls confidence tier assignment:
#'   \code{<=5} Confirmed, \code{<=10} Probable, \code{<=20} Tentative,
#'   \code{<=30} Questionable, \code{>30} rejected. Default \code{30}.
#' @param nl_ppm Numeric. ppm gate for neutral-loss matching in both phases.
#'   Default \code{5}.
#' @param rt_threshold Numeric. RT tolerance (seconds) for Phase 1
#'   within-group NL matching. Default \code{3}.
#' @param pi_rt_threshold Numeric. RT tolerance (seconds) for matching the
#'   precursor \code{[M+H]+}/\code{[M-H]-} feature against the library
#'   compound RT. Features outside this window are not considered as
#'   precursor matches. Default \code{30}.
#' @param rt_exp_threshold Numeric. RT window (seconds) around the compound
#'   reference RT for Phase 2 pool search. Default \code{5}.
#' @param ppm_exp Numeric. ppm gate for Phase 2 adduct search. Default \code{5}.
#' @param iso_threshold Numeric. Isotope mass delta (12C/13C spacing).
#'   Default \code{1.003355}.
#' @param max_charge Integer. Maximum charge state for isotope detection.
#'   Default \code{3L}.
#' @param max_isotope Integer. Maximum isotope peaks per charge state.
#'   Default \code{4L}.
#'
#' @return A \code{\link{PeakAnnotationFromLibraryMatch-class}} object.
#'   Pass to \code{\link{annotate}()} to run annotation.
#'
#' @examples
#' \dontrun{
#' obj <- createPeakAnnotationFromLibraryMatch(
#'   library_matches    = my_library_matches,
#'   FT            = my_feature_table,
#'   adducts_df         = adduct_table,
#'   neutral_loss_table = nl_table,
#'   polarity           = "positive"
#' )
#' result <- annotate(obj)
#' }
#'
#' @seealso \code{\link{annotate}()},
#'   \code{\link{PeakAnnotationFromLibraryMatch-class}}
#' @export
createPeakAnnotationFromLibraryMatch <- function(
    library_matches,
    FT            = data.frame(),
    adducts_df,
    neutral_loss_table,
    polarity           = "positive",
    mz_threshold       = 0.005,
    ppm                = 30,
    nl_ppm             = 5.0,
    rt_threshold       = 3.0,
    pi_rt_threshold    = 30.0,
    rt_exp_threshold   = 5.0,
    ppm_exp            = 5.0,
    iso_threshold      = 1.003355,
    max_charge         = 3L,
    max_isotope        = 4L
) {
  if (!is.list(library_matches))
    stop("`library_matches` must be a list of spectral library match results.")
  if (!is.data.frame(FT))
    stop("`FT` must be a data.frame (use data.frame() to skip Phase 2).")
  if (!is.data.frame(adducts_df))
    stop("`adducts_df` must be a data.frame.")
  if (!is.data.frame(neutral_loss_table))
    stop("`neutral_loss_table` must be a data.frame.")

  new("PeakAnnotationFromLibraryMatch",
      library_matches    = library_matches,
      FT            = FT,
      adducts_df         = adducts_df,
      neutral_loss_table = neutral_loss_table,
      polarity           = polarity,
      mz_threshold       = mz_threshold,
      ppm                = ppm,
      nl_ppm             = nl_ppm,
      rt_threshold       = rt_threshold,
      pi_rt_threshold    = pi_rt_threshold,
      rt_exp_threshold   = rt_exp_threshold,
      ppm_exp            = ppm_exp,
      iso_threshold      = iso_threshold,
      max_charge         = as.integer(max_charge),
      max_isotope        = as.integer(max_isotope)
  )
}
