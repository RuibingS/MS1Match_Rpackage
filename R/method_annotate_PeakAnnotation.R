#' Create a PeakAnnotation object
#'
#' @description
#' Construct a [`PeakAnnotation-class`] S4 object bundling the feature table,
#' compound list, adduct definitions, and neutral-loss definitions for downstream
#' peak annotation. RT-constrained vs. unconstrained annotation is decided later
#' by [annotate()] based on the presence and type of an RT column in `FT`.
#'
#' @param compounds A named `list` of parsed compounds. Each element should
#'   represent one compound and may include fields such as `Name`, `ExactMass`,
#'   `RetTime`, and a `Peaks` data.frame with columns `mz` and `rel_intensity`.
#' @param FT A `data.frame` of feature data. Must include at least one m/z column
#'   among `mzmed`, `mz`, or `m/z`. An RT column (`rtmed`, `rt`, or `retention_time`)
#'   is optional but enables RT-constrained annotation.
#' @param adducts_df A `data.frame` of adduct definitions (e.g. columns
#'   `Ion.name`, `Ion.polarity`, `Ion.mass`).
#' @param neutral_loss_table A `data.frame` of neutral losses (e.g. columns
#'   `Accurate.Mass`, `Formula`, `Pos`, `Neg`).
#' @param polarity Character; `"positive"` or `"negative"`. (Will be normalized to
#'   `"pos"` / `"neg"` internally.)
#' @param mz_threshold Numeric; absolute m/z tolerance in Da. Default `0.005`.
#' @param PI_rt_threshold Numeric; RT tolerance (seconds) for protonated/deprotonated
#'   ion (`"[M+H]+"` or `"[M-H]-"`).
#' @param rt_threshold Numeric; RT tolerance (seconds) for general matching.
#'   Default `3`.
#' @param iso_threshold Numeric; isotopic mass delta. Default `1.003355`.
#' @param max_charge Integer; maximum charge state to consider. Default `3L`.
#' @param max_isotope Integer; maximum number of isotope peaks to consider. Default `4L`.
#'
#' @return An S4 object of class [`PeakAnnotation-class`].
#'
#' @seealso [annotate()] to run the annotation workflow; also
#'   the C++ backend `annotate_with_rt_cpp()` used by RT-constrained matching.
#'
#' @examples
#' \dontrun{
#' pa <- createPeakAnnotation(
#'   compounds = PA14_library,
#'   FT = result_PA14_FT,
#'   adducts_df = adduct_table,
#'   neutral_loss_table = neutral_loss_table,
#'   polarity = "positive",
#'   mz_threshold = 0.005,
#'   PI_rt_threshold = 60,
#'   rt_threshold = 3,
#'   iso_threshold = 1.003355,
#'   max_charge = 3L,
#'   max_isotope = 4L
#' )
#' res <- annotate(pa)
#' }
#'
#' @export
createPeakAnnotation <- function(compounds, FT, adducts_df, neutral_loss_table,
                                 polarity, intensity,
                                 mz_threshold = 0.005,
                                 PI_rt_threshold = 5,
                                 rt_threshold = 3,
                                 iso_threshold = 1.003355,
                                 max_charge = 3L,
                                 max_isotope = 4L) {
  if (!is.list(compounds)) stop("`compounds` must be a list.")
  if (!is.data.frame(FT)) stop("`FT` must be a data.frame.")


  new("PeakAnnotation",
      compounds_list = compounds,
      FT = FT,
      adducts_df = adducts_df,
      neutral_loss_table = neutral_loss_table,
      polarity = polarity,

      mz_threshold = mz_threshold,
      PI_rt_threshold = PI_rt_threshold,
      rt_threshold = rt_threshold,
      iso_threshold = iso_threshold,
      max_charge = as.integer(max_charge),
      max_isotope = as.integer(max_isotope))
}


#' @export
setMethod("annotate", "PeakAnnotation", function(object) {
  has_rt <- any(c("rtmed", "rt", "retention_time") %in% colnames(object@FT)) &&
    is.numeric(object@FT[[intersect(c("rtmed", "rt", "retention_time"), colnames(object@FT))[1]]])

  if (has_rt) {
    message("Using RT-constrained annotation.")
    annotate_has_rt(object)
  } else {
    message("No RT info. Running unconstrained annotation.")
    annotate_no_rt(object)
  }
})
