# =============================================================================
# Calculate exact mass
# Standalone exported wrapper around the C++ monoisotopic mass calculator.
# Kept in its own file because it is independently useful — it can be called
# directly by users outside of any library reading or conversion workflow.
# =============================================================================


#' Calculate monoisotopic exact mass from molecular formulae
#'
#' @description
#' Fast C++ implementation of monoisotopic mass calculation. No external
#' packages required — replaces the previous \code{MetaboCoreUtils} dependency.
#'
#' @param formulas Character vector of molecular formulae
#'   (e.g. \code{c("C6H12O6", "C8H10N4O2")}).
#'
#' @return Numeric vector of monoisotopic masses, in the same order as
#'   \code{formulas}. Returns \code{NA} for unknown elements or
#'   unparseable formulae.
#'
#' @examples
#' calculate_exact_mass(c("C6H12O6", "C8H10N4O2", "C31H36N2O11"))
#' # [1] 180.0634 194.0798 592.2375
#'
#' calculate_exact_mass("C9H17NO6")
#' # [1] 235.1056
#'
#' @export
calculate_exact_mass <- function(formulas) {
  calculate_exact_mass_cpp(as.character(formulas))
}
