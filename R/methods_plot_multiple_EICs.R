#' @description plot_multiple_EICs
setMethod("plot_multiple_EICs",
          signature(object = "EICData"),
          function(object, feature_indices) {
            library(ggplot2)
            library(dplyr)

            # Extract the parsed_EIC list
            all_EIC_list <- object@parsed_EIC

            # Combine data for selected features
            combined_df <- lapply(feature_indices, function(i) {
              df <- all_EIC_list[[i]]$EIC_df
              df$sample <- rownames(df)               # rownames -> sample column
              df$feature <- paste0("Feature ", i)      # label feature
              df
            }) %>%
              bind_rows() %>%
              filter(!is.na(intensity) & intensity > 0)

            # Plot
            p <- ggplot(combined_df, aes(
              x = rt,
              y = log10(intensity),
              color = feature,
              group = feature
            )) +
              geom_line() +
              geom_point(size = 1) +
              labs(
                title = paste("Overlayed EICs for Features",
                              paste(feature_indices, collapse = ", ")),
                x = "Retention Time (rt)",
                y = "log10(Intensity)",
                color = "Feature"
              ) +
              theme_minimal()
            p
          })
