# example data frame: L-Phenylalanine
# ms_plot_df <- data.frame(
#   mz = c(77.0385, 91.0539, 93.0696, 102.0463, 103.0542, 104.0575,
#          107.0490, 118.0650, 119.0729, 120.0106, 120.0810, 121.0650,
#          121.0841, 122.0874, 128.0194, 131.0491, 132.0525, 148.0756,
#          149.0598, 150.0633, 166.0866, 167.0897, 168.0918, 188.0679),
#   intensity = c(1, 3, 3, 2, 57, 5, 7, 2, 2, 1, 620, 1, 52, 2,
#                 1, 12, 1, 2, 18, 2, 999, 102, 8, 4)
# )
# ms_plot_df$rel_intensity <-
#   ms_plot_df$intensity / max(ms_plot_df$intensity) * 100
# plot(
#   ms_plot_df$mz,
#   ms_plot_df$rel_intensity,
#   type = "h",
#   lwd = 2,
#   xlab = "m/z",
#   ylab = "Relative Intensity (%)",
#   main = "L-Phenylalanine MS1 Spectrum"
# )
# png(
#   "L-Phenylalanine_ms_reference_spectrum.png",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 600
# )
# plot(
#   ms_plot_df$mz,
#   ms_plot_df$rel_intensity,
#   type = "h",
#   lwd = 2,
#   xlab = "m/z",
#   ylab = "Relative Intensity (%)",
#   main = "L-Phenylalanine MS1 Spectrum"
# )
# dev.off()
