#' Create a DirectTargetMatch object
#'
#' @description
#' Constructor for the \code{DirectTargetMatch} S4 class. Validates inputs
#' and provides sensible defaults before calling \code{new()}.
#'
#' @param FT              A \code{data.frame} of detected features. Must contain
#'   a \code{feature_name} column and an m/z column (\code{mzmed}, \code{mz},
#'   or \code{m/z}). An RT column (\code{rtmed} or \code{rt}) is optional but
#'   enables RT gating when the target list also has RT values.
#' @param Comp_data       A \code{data.frame} target list. Accepted formats:
#'   HMDB (\code{MONO_MASS}, \code{[M+H]+}, \code{[M-H]-}), plate-based
#'   (\code{Neutral mass}, \code{M+H}, \code{M-H}, \code{CNAME}), PA14
#'   (\code{exactMass}, \code{[M+H]+}, \code{[M-H]-}, \code{RT}), or any
#'   format supported by \code{normalize_target_list()}.
#' @param PIon            Primary ion: \code{"[M+H]+"} (positive mode) or
#'   \code{"[M-H]-"} (negative mode). Default \code{"[M+H]+"}.
#' @param diff_mz_thr     Numeric. Absolute m/z tolerance in Da. Default
#'   \code{0.005} (5 mDa).
#' @param ppm             Numeric. ppm m/z tolerance. Default \code{10}.
#' @param rt_thr_compound Numeric. RT window in seconds for feature vs compound
#'   RT gating. Default \code{30}. Set to \code{Inf} to disable entirely.
#'
#' @return A \code{DirectTargetMatch} S4 object ready to pass to
#'   \code{PImatch_target()}.
#'
#' @examples
#' \dontrun{
#' # Works with any supported target list format
#' obj <- createDirectTargetMatch(
#'   FT              = xcms_output_ft,
#'   Comp_data       = plate_target_list,  # or PA14_target_list,  hmdb_target_list
#'   PIon            = "[M+H]+",
#'   diff_mz_thr     = 0.005,
#'   ppm             = 10,
#'   rt_thr_compound = 30
#' )
#' obj <- PImatch_target(obj)
#' df  <- flatten(obj)
#' }
#'
#' @export
createDirectTargetMatch <- function(FT,
                                    Comp_data,
                                    PIon            = "[M+H]+",
                                    diff_mz_thr     = 0.005,
                                    ppm             = 10,
                                    rt_thr_compound = 30) {

  # ---- input validation ----
  if (!is.data.frame(FT))
    stop("`FT` must be a data.frame.")
  if (!is.data.frame(Comp_data))
    stop("`Comp_data` must be a data.frame.")
  if (!PIon %in% c("[M+H]+", "[M-H]-"))
    stop("`PIon` must be \"[M+H]+\" or \"[M-H]-\".")
  if (!is.numeric(diff_mz_thr) || diff_mz_thr <= 0)
    stop("`diff_mz_thr` must be a positive number.")
  if (!is.numeric(ppm) || ppm <= 0)
    stop("`ppm` must be a positive number.")
  if (!is.numeric(rt_thr_compound) || rt_thr_compound <= 0)
    stop("`rt_thr_compound` must be a positive number (use Inf to disable RT gating).")
  if (!"feature_name" %in% colnames(FT))
    stop("`FT` must contain a 'feature_name' column.")

  new("DirectTargetMatch",
      FT              = as.data.frame(FT),
      Comp_data       = as.data.frame(Comp_data),
      PIon            = PIon,
      diff_mz_thr     = diff_mz_thr,
      ppm             = as.integer(ppm),
      rt_thr_compound = rt_thr_compound,
      results         = list()
  )
}
