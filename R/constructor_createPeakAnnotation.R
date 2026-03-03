#' Create a PeakAnnotation object
#'
#' @description
#' Constructor for the [`PeakAnnotation-class`] S4 class.
#' This function bundles the feature table, compound list, adduct definitions,
#' and neutral loss definitions into a single object for downstream annotation.
#'
#' @param compounds A named `list` of parsed compounds.
#'   Each element should represent one compound and may contain metadata
#'   (e.g., name, formula, exact mass) and peak data (m/z and intensity).
#' @param FT A `data.frame` of feature data.
#'   Must include at least a `mzmed` column; optionally `rtmed` or `rt` for
#'   retention time-constrained annotation.
#' @param adducts_df A `data.frame` of adduct definitions.
#' @param neutral_loss_table A `data.frame` of neutral loss definitions.
#' @param polarity Character string indicating ionization mode: `"positive"` or `"negative"`.
#' @param intensity Numeric vector of feature intensities (mean or representative values).
#' @param mz_threshold Numeric m/z tolerance for matching (in Da).
#' @param PI_rt_threshold Numeric retention time tolerance (in seconds) for protonated ion detection.
#' @param rt_threshold Numeric retention time tolerance (in seconds) for general matching.
#' @param iso_threshold Numeric mass difference for isotopes (default ~1.003355).
#' @param max_charge Integer. Maximum charge state to consider.
#' @param max_isotope Integer. Maximum number of isotope peaks to consider.
#'
#' @return An S4 object of class [`PeakAnnotation-class`].
#' @seealso [annotate()] for running the annotation workflow.
#'
#' @examples
#' \dontrun{
#' pa <- createPeakAnnotation(
#'   compounds = my_compound_list,
#'   FT = my_feature_table,
#'   adducts_df = adduct_table,
#'   neutral_loss_table = neutral_loss_table,
#'   polarity = "positive",
#'   intensity = my_feature_table$intensity_mean
#' )
#' annotated <- annotate(pa)
#' }
#'
#' @export
createPeakAnnotation <- function(compounds, FT, adducts_df, neutral_loss_table,
                                 polarity, intensity,
                                 mz_threshold = 0.005, PI_rt_threshold = 5,
                                 rt_threshold = 3, iso_threshold = 1.003355,
                                 max_charge = 3L, max_isotope = 4L) {
  if (!is.list(compounds)) {
    stop("`compounds` must be a list.")
  }

  new("PeakAnnotation",
      FT = FT,
      compounds_list = compounds,
      adducts_df = adducts_df,
      neutral_loss_table = neutral_loss_table,
      polarity = polarity,
      intensity = intensity,
      mz_threshold = mz_threshold,
      PI_rt_threshold = PI_rt_threshold,
      rt_threshold = rt_threshold,
      iso_threshold = iso_threshold,
      max_charge = as.integer(max_charge),
      max_isotope = as.integer(max_isotope))
}

