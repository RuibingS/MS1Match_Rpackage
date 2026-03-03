#' Create a PeakAnnotationFromLibraryMatch object
#'
#' @description
#' Construct a [PeakAnnotationFromLibraryMatch-class] object for annotating features
#' **after a spectral library match has already been performed** (MS1 or MS2).
#'
#' Annotation is performed **within each matched group** using:
#' \itemize{
#'   \item The neutral mass `Compound_M` to compute theoretical adduct m/z values.
#'   \item Pairwise m/z differences among `Matched_features_mz` for neutral-loss detection.
#'   \item Within-group isotope detection.
#' }
#'
#' This avoids rescanning the full library against the full feature table.
#'
#' @param library_matches A `list` of spectral library match results
#'   (e.g. `direct_match_MS1`, `direct_match_MS2`).
#'   Each element should contain:
#'   \itemize{
#'     \item `Compound_name`
#'     \item `Compound_M`
#'     \item `Ion_type`
#'     \item `Lib_group`
#'     \item `Matched_features_to_peaks`
#'     \item `Matched_features_mz`
#'   }
#'   Optional:
#'   \itemize{
#'     \item `Matched_features_rt`
#'     \item `Matched_features_int`
#'   }
#' @param adducts_df Adduct definition table (`Ion.name`, `Ion.polarity`, `Ion.mass`, optional `Charge`)
#' @param neutral_loss_table Neutral loss table (`Accurate.Mass`, `Formula`, `Pos`, `Neg`)
#' @param polarity `"positive"` or `"negative"` (normalized internally)
#' @param mz_threshold Absolute m/z tolerance (Da). Default `0.005`
#' @param ppm ppm tolerance for adduct matching (effective tolerance is `max(mz_threshold, ppm*mz*1e-6)`)
#' @param rt_threshold RT tolerance (seconds) for within-group neutral-loss and isotope detection
#' @param iso_threshold Isotope mass delta. Default `1.003355`
#' @param max_charge Maximum charge state. Default `3L`
#' @param max_isotope Maximum isotope number. Default `4L`
#'
#' @return An S4 object of class [PeakAnnotationFromLibraryMatch-class].
#'
#' @seealso [annotate()]
#'
#' @export
createPeakAnnotationFromLibraryMatch <- function(library_matches,
                                                 adducts_df,
                                                 neutral_loss_table,
                                                 polarity = "positive",
                                                 mz_threshold = 0.005,
                                                 ppm = 20,
                                                 rt_threshold = 3,
                                                 iso_threshold = 1.003355,
                                                 max_charge = 3L,
                                                 max_isotope = 4L) {

  if (!is.list(library_matches)) {
    stop("`library_matches` must be a list of spectral library match results.")
  }
  if (!is.data.frame(adducts_df)) stop("`adducts_df` must be a data.frame.")
  if (!is.data.frame(neutral_loss_table)) stop("`neutral_loss_table` must be a data.frame.")

  new("PeakAnnotationFromLibraryMatch",
      library_matches = library_matches,
      adducts_df = adducts_df,
      neutral_loss_table = neutral_loss_table,
      polarity = polarity,
      mz_threshold = mz_threshold,
      ppm = ppm,
      rt_threshold = rt_threshold,
      iso_threshold = iso_threshold,
      max_charge = as.integer(max_charge),
      max_isotope = as.integer(max_isotope))
}
