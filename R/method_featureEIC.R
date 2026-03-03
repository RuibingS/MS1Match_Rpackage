#' @rdname feature_EIC_process
#' @return The updated EICData object with parsed_EIC slot which has chrom_peaks_df and one EIC_df with maximum "into"
#' @export
setMethod("feature_EIC_process", "EICData",
          function(object, nrows, ncols) {

            input_list <-  as.vector(object@eic_data@.Data, mode = "list")

            # Call the Rcpp function
            parsed <- feature_EIC_process_fun_cpp(input_list, nrows, ncols)
            # Name
            names(parsed) <- sprintf("FT%04d", seq_along(parsed))
            # Store result back into the object
            object@parsed_EIC <- parsed
            return(object)
          })
