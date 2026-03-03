#' XCMS preprocessing
#' @import MSnbase
#' @import xcms
#' @import Rcpp
#' @importFrom MSnbase readMSData
#' @param files_dir A file directory of mzXML or mzXL files.
#' @param peakwidth CentWave parameter, refer to XCMS.
#' @param ppm `CentWave parameter`, refer to XCMS.
#' @param prefilter `CentWave parameter`, refer to XCMS.
#' @param noise `CentWave parameter`, refer to XCMS.
#' @param bw `PeakDensityParam` parameter, refer to XCMS.
#' @param sample_name The name of samples. Please using "_A","_B" and "_C" or "_1","_2","_3" for replicates.
#' @param sample_group The sample group names.
#' @param polarity The ion mode. "positive" or "negative".
#' @param minFraction `PeakDensityParam`paramters, refer to XCMS.
#' @return It returns a list with a feature table and EIC object.
#' @examples
#' files_dir <- "C:/Users/Downloads/PA14_mzXML_files/MS1_files_mzXML"
#' peakwidth <- c(10, 25)
#' ppm <- 10
#' snthresh <- 10
#' mzdiff <- 0.01
#' prefilter <- c(3, 1000)
#' noise <- 100
#' bw <- 5
#' minFraction <- 0.5
#' sample_name <- c(
#'   paste0("AR351", "_A"), paste0("AR351", "_B"), paste0("AR351", "_C"),
#'   paste0("Cipro", "_A"), paste0("Cipro", "_B"), paste0("Cipro", "_C"),
#'   paste0("Coumermycin", "_A"), paste0("Coumermycin", "_B"), paste0("Coumermycin", "_C"),
#'   paste0("GT3", "_A"), paste0("GT3", "_B"), paste0("GT3", "_C"),
#'   paste0("Levo", "_A"), paste0("Levo", "_B"), paste0("Levo", "_C"),
#'   paste0("Lome", "_A"), paste0("Lome", "_B"), paste0("Lome", "_C"),
#'   paste0("Novo", "_A"), paste0("Novo", "_B"), paste0("Novo", "_C"),
#'   paste0("PA14wt", "_A"), paste0("PA14wt", "_B"), paste0("PA14wt", "_C")
#' )
#' sample_group <- c(
#'   rep("AR351", 3), rep("Cipro", 3), rep("Coumermycin", 3),
#'   rep("GT3", 3), rep("Levo", 3), rep("Lome", 3),
#'   rep("Novo", 3), rep("PA14wt", 3)
#' )
#' polarity <- "positive"


#' Run full XCMS preprocessing with beta + CPC metrics
#'
#' @param files_dir Directory with mzML/mzXML files.
#' @param cwp CentWaveParam (optional, default reasonable).
#' @param pdp PeakDensityParam (optional, default reasonable).
#' @param cpc_param CPC parameters (optional, default reasonable).
#' @param sample_name Sample names.
#' @param sample_group Sample groups.
#' @return List with objects and merged feature table.
#' @export
runXCMS <- function(
    files_dir,
    cwp = NULL,
    pdp = NULL,
    cpc_param = NULL,
    sample_name = NULL,
    sample_group = NULL
) {
  MS_files <- list.files(files_dir, pattern = "\\.(mzXML|mzML)$", full.names = TRUE)
  if (!length(MS_files)) stop("No mzXML/mzML files in ", files_dir)

  if (is.null(sample_group)) sample_group <- seq_along(MS_files)
  if (is.null(sample_name))  sample_name  <- basename(MS_files)

  pd <- data.frame(sample_name = sample_name,
                   sample_group = sample_group,
                   stringsAsFactors = FALSE)

  # defaults
  if (is.null(cwp)) {
    cwp <- CentWaveParam(peakwidth = c(10, 25), ppm = 10,
                         snthresh = 10, mzdiff = 0.01,
                         prefilter = c(3, 1000), noise = 200,
                         verboseBetaColumns = TRUE)
  }
  if (is.null(pdp)) {
    pdp <- PeakDensityParam(sampleGroups = sample_group,
                            bw = 5, minFraction = 0.5)
  }
  if (is.null(cpc_param)) {
    cpc_param <- cpc::cpcProcParam(min_sn = 10, min_pts = 10,
                                   min_intensity = 2000)
  }

  # --- XcmsExperiment pipeline (for beta) ---
  raw_exp <- MsExperiment::readMsExperiment(spectraFiles = MS_files, sampleData = pd)
  x_exp <- findChromPeaks(raw_exp, param = cwp, msLevel = 1L)
  if (length(MS_files) > 1) x_exp <- adjustRtime(x_exp, ObiwarpParam())
  x_exp <- groupChromPeaks(x_exp, param = pdp)
  x_exp <- fillChromPeaks(x_exp, ChromPeakAreaParam())

  fea_def <- featureDefinitions(x_exp)
  feature_chroms <- featureChromatograms(x_exp, features = seq_len(nrow(fea_def)))

  # --- XCMSnExp pipeline (for CPC) ---
  raw_nexp <- MSnbase::readMSData(files = MS_files,
                                  pdata = new("AnnotatedDataFrame", pd),
                                  mode = "onDisk")
  x_nexp <- findChromPeaks(raw_nexp, param = cwp, msLevel = 1L)
  cpc_obj <- cpc::filter_xcms_peaklist(xd = x_nexp,
                                       return_type = "cpc",
                                       param = cpc_param)

  # --- Merge metrics ---
  fd_beta <- link_feature_CP_wide(x_exp, append_intensity = FALSE)
  fd_cpc  <- link_feature_CPC_wide(x_nexp, cpc_obj, append_intensity = FALSE)
  fint    <- featureValues(x_exp, value = "into")



  # --- Merge metrics ---
  message("Merging feature-level metrics")

  fd_beta <- link_feature_CP_wide(x_exp, append_intensity = FALSE)
  message("fd_beta rows: ", nrow(fd_beta))

  fd_cpc  <- link_feature_CPC_wide(x_nexp, cpc_obj, append_intensity = FALSE)
  message("fd_cpc rows: ", nrow(fd_cpc))

  fint    <- featureValues(x_exp, value = "into")
  message("fint rows: ", nrow(fint))


  fd_final <- cbind(
    fd_beta,
    fd_cpc,
    fint
  )

  message("fd_final rows: ", nrow(fd_final))

  list(
    x_exp          = x_exp,
    x_nexp         = x_nexp,
    cpc            = cpc_obj,
    features       = fd_final,
    feature_chroms = feature_chroms
  )
}

# runXCMS <- function(files_dir,
#                     cwp = NULL,
#                     pdp = NULL,
#                     cpc_param = NULL,
#                     sample_name = NULL,
#                     sample_group = NULL) {
#   ## -----------------------------
#   ## 0. File handling
#   ## -----------------------------
#   MS_files <- list.files(files_dir,
#                          pattern = "\\.(mzXML|mzML)$",
#                          full.names = TRUE)
#   if (!length(MS_files)) stop("No mzXML/mzML files found in ", files_dir)
#
#   if (is.null(sample_group)) sample_group <- seq_along(MS_files)
#   if (is.null(sample_name))  sample_name  <- basename(MS_files)
#
#   pd <- data.frame(sample_name = sample_name,
#                    sample_group = sample_group,
#                    stringsAsFactors = FALSE)
#
#
#   ## -----------------------------
#   ## 1. Build parameter objects
#   ## -----------------------------
#   if (is.null(cwp)) {
#     message("Creating CentWaveParam with defaults")
#     cwp <- CentWaveParam(peakwidth = c(10, 25), ppm = 10,
#                          snthresh = 10, mzdiff = 0.01,
#                          prefilter = c(3, 1000), noise = 200,
#                          verboseBetaColumns = TRUE)
#   }
#   if (is.null(pdp)) {
#     message("Creating PeakDensityParam with sample groups")
#     pdp <- PeakDensityParam(sampleGroups = sample_group,
#                             bw = 5, minFraction = 0.5)
#   }
#   if (is.null(cpc_param)) {
#     message("Creating CPC parameters with defaults")
#     cpc_param <- cpc::cpcProcParam(min_sn = 10,
#                                    min_pts = 10,
#                                    min_intensity = 2000)
#   }
#
#   ## -----------------------------
#   ## 2. XcmsExperiment workflow
#   ## -----------------------------
#   message("Starting XcmsExperiment pipeline")
#
#   raw_exp <- readMsExperiment(spectraFiles = MS_files, sampleData = pd)
#
#   message(" - Peak detection")
#   x_exp <- findChromPeaks(raw_exp, param = cwp, msLevel = 1L)
#
#   if (length(MS_files) > 1) {
#     message(" - Retention time alignment (obiwarp)")
#     x_exp <- adjustRtime(x_exp, ObiwarpParam())
#   }
#
#   message(" - Grouping features")
#   x_exp <- groupChromPeaks(x_exp, param = pdp)
#
#   message(" - Gap filling")
#   x_exp <- fillChromPeaks(x_exp, ChromPeakAreaParam())
#
#   fea_def <- featureDefinitions(x_exp)
#   feature_chroms <- featureChromatograms(x_exp,
#                                          features = seq_len(nrow(fea_def)))
#
#   ## -----------------------------
#   ## 3. XCMSnExp workflow for CPC
#   ## -----------------------------
#   message("Running CPC (requires XCMSnExp backend)")
#
#   raw_nexp <- MSnbase::readMSData(files = MS_files,
#                                   pdata = new("AnnotatedDataFrame", pd),
#                                   mode = "onDisk")
#   x_nexp <- findChromPeaks(raw_nexp, param = cwp, msLevel = 1L)
#
#   #if (length(MS_files) > 1) {
#    # x_nexp <- adjustRtime(x_nexp, ObiwarpParam())
#   #}
#
#   #x_nexp <- groupChromPeaks(x_nexp, param = pdp)
#   #x_nexp <- fillChromPeaks(x_nexp, ChromPeakAreaParam())
#
#   message(" - Characterizing peaks with CPC")
#   cpc_obj <- cpc::filter_xcms_peaklist(
#     xd = x_nexp,
#     return_type = "cpc",
#     param = cpc_param
#   )
#   #cpc_tab <- cpc::cpt(cpc_obj)
#
#
#   ## -----------------------------
#   ## 4. Merge metrics
#   ## -----------------------------
#   message("Merging feature-level metrics")
#
#   fd_beta <- link_feature_CP_wide(
#     x_exp,
#     metrics = c("beta_cor", "beta_snr"),
#     append_intensity = FALSE
#   )
#
#   fd_cpc <- link_feature_CPC_wide(
#     x_nexp,
#     cpc_obj,
#     append_intensity = FALSE
#   )
#
#   fint <- featureValues(x_exp, value = "into")
#
#   fd_final <- cbind(
#     fd_beta,
#     fd_cpc[, setdiff(colnames(fd_cpc), colnames(fd_beta))],
#     fint
#   )
#
#   message("Pipeline finished")
#
#   ## -----------------------------
#   return(list(
#     x_exp          = x_exp,
#     x_nexp         = x_nexp,
#     cpc            = cpc_obj,
#     features       = fd_final,
#     feature_chroms = feature_chroms
#   ))
# }


