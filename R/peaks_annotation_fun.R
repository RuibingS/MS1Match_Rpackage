#' Annotate mass spec peaks
#' @description
#' This function annotates feature table peaks using known compounds, adducts,
#' and neutral loss information.
#'
#' @param compounds_list A list of processed compounds.
#' @param FT A data frame of the feature table, sorted by ascending retention time.
#' @param adducts_df A data frame of adduct table
#' @param neutral_loss_table A data frame of neutral loss table
#' @param polarity Either "pos" or "neg".
#' @param intensity A numeric vector of mean intensity. It is used to calculate the C13 isotope intensity ratio
#' @param mz_threshold m/z matching threshold (in Da).
#' @param PI_rt_threshold Precursor ion retention time tolerance.
#' @param rt_threshold General RT matching tolerance.
#' @param iso_threshold Isotope m/z tolerance.
#' @param max_charge Maximum charge considered for adducts.
#' @param max_isotope Maximum number of isotopes to consider.
#'
#' @details
#' The feature table `FT` must be sorted by `rtmed` or `rt` in ascending order.
#' If not, the function will auto-fix with a warning.
#'
#' @return A list with annotated results.
#' @export
annotate_peaks_fun <- function(compounds_list, FT, adducts_df, neutral_loss_table,
                               polarity, intensity,
                               mz_threshold, PI_rt_threshold,
                               rt_threshold, iso_threshold,
                               max_charge, max_isotope) {
  # placeholder content
  list()
}
