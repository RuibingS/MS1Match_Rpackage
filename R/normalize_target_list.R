#' Normalize a target list data frame to a standard schema
#'
#' @description
#' Accepts target lists from different sources
#' and maps heterogeneous column names to a unified schema with columns:
#' \verb{Name}, \verb{Formula}, \verb{exactMass}, \verb{[M+H]+}, \verb{[M-H]-},
#' and optionally \verb{RT}.
#'
#' Supported input column name variants:
#' \describe{
#'   \item{Name}{CNAME, Name}
#'   \item{Formula}{FORMULA, Formula, Formula (as seen in ...), any column starting with "Formula"}
#'   \item{exactMass}{exactMass, MONO_MASS, Neutral mass, Neutral.mass}
#'   \item{[M+H]+}{[M+H]+, M+H, MH, MplusH}
#'   \item{[M-H]-}{[M-H]-, M-H, MH_neg, MminusH}
#'   \item{RT}{RT, RetentionTime, rt}
#' }
#'
#' If \verb{[M+H]+} or \verb{[M-H]-} are absent but \verb{exactMass} is present,
#' the ion m/z values are computed automatically.
#'
#' @param target_list A \code{data.frame} or \code{tibble} target list
#' @return A \code{data.frame} with standardized column names
#' @export
normalize_target_list <- function(target_list) {

  df <- as.data.frame(target_list, stringsAsFactors = FALSE)
  cols <- colnames(df)

  # ---- helper: find first matching column (case-insensitive, partial ok) ----
  find_col <- function(patterns) {
    for (p in patterns) {
      # exact match first (case-insensitive)
      hit <- which(tolower(cols) == tolower(p))
      if (length(hit) > 0) return(cols[hit[1]])
      # startsWith match
      hit <- which(startsWith(tolower(cols), tolower(p)))
      if (length(hit) > 0) return(cols[hit[1]])
    }
    return(NULL)
  }

  out <- data.frame(row.names = seq_len(nrow(df)))

  # ---- Name ----
  name_col <- find_col(c("CNAME", "Name", "name"))
  if (is.null(name_col))
    stop("Cannot find compound name column. Expected one of: CNAME, Name.")
  out[["Name"]] <- as.character(df[[name_col]])

  # ---- Formula ----
  formula_col <- find_col(c("FORMULA", "Formula"))
  if (!is.null(formula_col))
    out[["Formula"]] <- as.character(df[[formula_col]])
  else
    out[["Formula"]] <- NA_character_

  # ---- Exact mass ----
  mass_col <- find_col(c("exactMass", "MONO_MASS", "Neutral mass", "Neutral.mass",
                         "MonoisotopicMass", "mono_mass"))
  if (is.null(mass_col))
    stop("Cannot find exact mass column. Expected one of: exactMass, MONO_MASS, Neutral mass.")
  out[["exactMass"]] <- as.numeric(df[[mass_col]])

  # ---- [M+H]+ ----
  mh_pos_col <- find_col(c("[M+H]+", "M+H", "MplusH", "MH_pos"))
  if (!is.null(mh_pos_col)) {
    out[["[M+H]+"]] <- as.numeric(df[[mh_pos_col]])
  } else {
    message("'[M+H]+' column not found — computing from exactMass + 1.007276")
    out[["[M+H]+"]] <- out[["exactMass"]] + 1.007276
  }

  # ---- [M-H]- ----
  mh_neg_col <- find_col(c("[M-H]-", "M-H", "MminusH", "MH_neg"))
  if (!is.null(mh_neg_col)) {
    out[["[M-H]-"]] <- as.numeric(df[[mh_neg_col]])
  } else {
    message("'[M-H]-' column not found — computing from exactMass - 1.007276")
    out[["[M-H]-"]] <- out[["exactMass"]] - 1.007276
  }

  # ---- RT (optional) ----
  rt_col <- find_col(c("RT", "rt", "RetentionTime", "retention_time"))
  if (!is.null(rt_col))
    out[["RT"]] <- as.numeric(df[[rt_col]])
  else
    out[["RT"]] <- NA_real_

  out
}
