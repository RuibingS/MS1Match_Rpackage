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
#' @export
#'
#' @import xcms
#' @import MSnbase

runXCMS <- function(files_dir,
                    peakwidth,
                    ppm,
                    snthresh,
                    mzdiff,
                    prefilter,
                    noise,
                    bw,
                    sample_name,
                    sample_group,
                    minFraction) {



  MS_files <- list.files(files_dir,
                             recursive = FALSE,
                             pattern = "\\.mzXML$|\\.mzML$",
                             full.names = TRUE)


  if (is.null(sample_group)) {
    sample_group <- rep(1:length(MS_files))
  }
  if (is.null(sample_name)) {
    sample_name <- rep(1:length(MS_files))
  }

  pd <- data.frame(sample_name = sample_name,
                   sample_group = sample_group,
                   stringsAsFactors = FALSE)

  # Load raw MS data
  raw_data <- MSnbase::readMSData(files = MS_files,
                         pdata = new("AnnotatedDataFrame", pd),
                         mode = "onDisk")

  # Peak detection with CentWave
  cwp <- CentWaveParam(peakwidth = peakwidth, ppm = ppm,
                       snthresh = snthresh, mzdiff = mzdiff,
                       prefilter = prefilter, noise = noise)

  xdata <- findChromPeaks(raw_data, param = cwp)

  # Filter MS1 level
  xdata_MS1 <- filterMsLevel(xdata, msLevel = 1L)


  # Retention time correction alignment
  if (length(MS_files) > 1) {
    xdata_MS1 <- adjustRtime(xdata_MS1, param = ObiwarpParam())
  }


  # Peak correspondence across samples
  pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                          minFraction = minFraction, bw = bw)

  xdata_MS1 <- groupChromPeaks(xdata_MS1, param = pdp)

  # Fill missing peaks
  xdata_MS1 <- fillChromPeaks(xdata_MS1, param = ChromPeakAreaParam())


  # feature definition
  fea_def <- featureDefinitions(xdata_MS1)

  # feature intensity value
  fea_int <- featureValues(xdata_MS1, value = "into")

  # feature chromatograms
  feature_chroms <- featureChromatograms(xdata_MS1, features = 1:nrow(fea_def))

  return(list(xdata=xdata_MS1,feature_def=fea_def,feature_int=fea_int,feature_chroms=feature_chroms))
}
