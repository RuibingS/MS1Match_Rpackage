#' Compare two MSP library lists for overlapping compounds
#'
#' @param lib1 A list of compounds (from import_msp)
#' @param lib2 A second list of compounds (from import_msp)
#'
#' @return A data.frame with matching compounds containing:
#'   - Formula
#'   - Name from lib1
#'   - Synonyms from lib1
#'   - Name from lib2
#'   - Synonyms from lib2
#'   - Match_type ("Formula+Name" or "Formula_only")
#'
#' @export
compare_libraries <- function(lib1, lib2) {

  get_field <- function(x, field) if (!is.null(x[[field]])) x[[field]] else NA_character_

  results <- lapply(seq_along(lib1), function(i) {
    f1 <- get_field(lib1[[i]], "Formula")
    n1 <- get_field(lib1[[i]], "Name")
    s1 <- unlist(lib1[[i]][grep("^Synon", names(lib1[[i]]))])
    if (length(s1) == 0) s1 <- NA_character_

    matches <- lapply(seq_along(lib2), function(j) {
      f2 <- get_field(lib2[[j]], "Formula")
      n2 <- get_field(lib2[[j]], "Name")
      s2 <- unlist(lib2[[j]][grep("^Synon", names(lib2[[j]]))])
      if (length(s2) == 0) s2 <- NA_character_

      if (!is.na(f1) && !is.na(f2) && f1 == f2) {
        all_names1 <- tolower(na.omit(c(n1, s1)))
        all_names2 <- tolower(na.omit(c(n2, s2)))

        matched <- any(vapply(all_names1, function(a) {
          any(vapply(all_names2, function(b) {
            grepl(a, b, fixed = TRUE) || grepl(b, a, fixed = TRUE)
          }, logical(1)))
        }, logical(1)))

        match_type <- if (matched) "Formula+Name" else "Formula_only"

        return(data.frame(
          Formula = f1,
          Name_in_house = n1,
          Synonyms_in_house = paste(na.omit(s1), collapse = "; "),
          Name_Massbank = n2,
          Synonyms_Massbank = paste(na.omit(s2), collapse = "; "),
          Match_type = match_type,
          stringsAsFactors = FALSE
        ))
      } else {
        return(NULL)
      }
    })

    do.call(rbind, matches[!vapply(matches, is.null, logical(1))])
  })

  do.call(rbind, results[!vapply(results, is.null, logical(1))])
}

#' Compare up to three spectral libraries by compound identity
#'
#' @description
#' Compare compounds across up to three libraries (e.g., Inhouse, MoNA, CSU–PMF)
#' based on formula and name/synonym overlap.
#'
#' @param lib1,lib2,lib3 Lists of compounds parsed with `import_msp()` or `import_msp_csu()`.
#' @param names Optional character vector of length 3 giving custom library names.
#'
#' @return A data.frame with matched formulas and their presence across libraries.
#'
#' @details
#' - Matching is first done by exact formula.
#' - Then, within matching formulas, name/synonym overlap is checked (case-insensitive).
#' - The output includes columns showing names and synonyms from each library.
#'
#' @examples
#' \dontrun{
#' comp_tbl <- compare_libraries_3(inhouse_lib, mona_lib, csu_lib)
#' head(comp_tbl)
#' }
#' @export
compare_libraries_3 <- function(lib1, lib2, lib3 = NULL,
                                names = c("Library1", "Library2", "Library3")) {

  get_field <- function(x, field) if (!is.null(x[[field]])) x[[field]] else NA_character_

  extract_info <- function(lib, label) {
    data.frame(
      Formula = sapply(lib, function(x) get_field(x, "Formula")),
      Name = sapply(lib, function(x) get_field(x, "Name")),
      Synonyms = sapply(lib, function(x) {
        syn <- unlist(x[grep("^SYNON", names(x), ignore.case = TRUE)])
        if (length(syn) == 0) NA_character_ else paste(syn, collapse = "; ")
      }),
      Source = label,
      stringsAsFactors = FALSE
    )
  }

  df1 <- extract_info(lib1, names[1])
  df2 <- extract_info(lib2, names[2])
  df3 <- if (!is.null(lib3)) extract_info(lib3, names[3]) else NULL

  libs <- list(df1, df2, df3)
  libs <- libs[!sapply(libs, is.null)]

  # --- merge by formula
  merged <- Reduce(function(x, y) merge(x, y, by = "Formula", all = TRUE, suffixes = c("", ".y")), libs)

  # --- name/synonym overlap detection
  merged$Match_type <- apply(merged, 1, function(row) {
    names_all <- tolower(na.omit(c(row[grep("Name", names(row))], row[grep("Synonyms", names(row))])))
    uniq_names <- unique(unlist(strsplit(names_all, ";\\s*")))
    if (length(uniq_names) > 1) {
      "Formula+Name"
    } else {
      "Formula_only"
    }
  })

  # --- presence flags
  for (nm in names[seq_along(libs)]) {
    merged[[nm]] <- ifelse(!is.na(merged[[paste0("Name.", nm)]]) |
                             !is.na(merged[[paste0("Synonyms.", nm)]]), TRUE, FALSE)
  }

  merged
}
