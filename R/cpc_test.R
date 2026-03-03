#' # test with the new xcms function
#'
#' #
# PA14_xcms_res <- runXCMS(files_dir = "C:/Users/rsh19/Downloads/PA14_mzXML_files/MS1_files_mzXML",
#                        sample_name = c(
#                          paste0("AR351", "_A"), paste0("AR351", "_B"), paste0("AR351", "_C"),
#                          paste0("Cipro", "_A"), paste0("Cipro", "_B"), paste0("Cipro", "_C"),
#                           paste0("Coumermycin", "_A"), paste0("Coumermycin", "_B"), paste0("Coumermycin", "_C"),
#                          paste0("GT3", "_A"), paste0("GT3", "_B"), paste0("GT3", "_C"),
#                            paste0("Levo", "_A"), paste0("Levo", "_B"), paste0("Levo", "_C"),
#                           paste0("Lome", "_A"), paste0("Lome", "_B"), paste0("Lome", "_C"),
#                            paste0("Novo", "_A"), paste0("Novo", "_B"), paste0("Novo", "_C"),
#                            paste0("PA14wt", "_A"), paste0("PA14wt", "_B"), paste0("PA14wt", "_C")
#                          ),
#                          sample_group = c(
#                            rep("AR351", 3), rep("Cipro", 3), rep("Coumermycin", 3),
#                            rep("GT3", 3), rep("Levo", 3), rep("Lome", 3),
#                            rep("Novo", 3), rep("PA14wt", 3)
#                          )
#
#                          )




#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' #'
#' #' library(MSnbase)
#' #' library(xcms)
#' #' library(MsExperiment)
#' #' library(BiocParallel)
#' #'
#' #' #' files_dir <- "C:/Users/Downloads/PA14_mzXML_files/MS1_files_mzXML"
#' #' #' peakwidth <- c(10, 25)
#' #' #' ppm <- 10
#' #' #' snthresh <- 10
#' #' #' mzdiff <- 0.01
#' #' #' prefilter <- c(3, 1000)
#' #' #' noise <- 100
#' #' #' bw <- 5
#' #' #' minFraction <- 0.5
#' #' #' sample_name <- c(
#' #' #'   paste0("AR351", "_A"), paste0("AR351", "_B"), paste0("AR351", "_C"),
#' #' #'   paste0("Cipro", "_A"), paste0("Cipro", "_B"), paste0("Cipro", "_C"),
#' #' #'   paste0("Coumermycin", "_A"), paste0("Coumermycin", "_B"), paste0("Coumermycin", "_C"),
#' #' #'   paste0("GT3", "_A"), paste0("GT3", "_B"), paste0("GT3", "_C"),
#' #' #'   paste0("Levo", "_A"), paste0("Levo", "_B"), paste0("Levo", "_C"),
#' #' #'   paste0("Lome", "_A"), paste0("Lome", "_B"), paste0("Lome", "_C"),
#' #' #'   paste0("Novo", "_A"), paste0("Novo", "_B"), paste0("Novo", "_C"),
#' #' #'   paste0("PA14wt", "_A"), paste0("PA14wt", "_B"), paste0("PA14wt", "_C")
#' #' #' )
#' #' #' sample_group <- c(
#' #' #'   rep("AR351", 3), rep("Cipro", 3), rep("Coumermycin", 3),
#' #' #'   rep("GT3", 3), rep("Levo", 3), rep("Lome", 3),
#' #' #'   rep("Novo", 3), rep("PA14wt", 3)
#' #' #' )
#' #' #' polarity <- "positive"
#' #'
#' #' # get example raw data file paths from package
#' fp_raw <- list.files( "C:/Users/rsh19/Downloads/PA14_mzXML_files/MS1_files_mzXML",
#'                       recursive = FALSE,
#'                       pattern = "\\.mzXML$|\\.mzML$",
#'                       full.names = TRUE)
#' #'
#' #' # setup multi-core processing for XCMS
#' #' #register(bpstart(SnowParam()))
#' #'
#' #' # setup metadata
#' pd <- data.frame(sample_name =c(
#'      paste0("AR351", "_A"), paste0("AR351", "_B"), paste0("AR351", "_C"),
#'      paste0("Cipro", "_A"), paste0("Cipro", "_B"), paste0("Cipro", "_C"),
#'      paste0("Coumermycin", "_A"), paste0("Coumermycin", "_B"), paste0("Coumermycin", "_C"),
#'      paste0("GT3", "_A"), paste0("GT3", "_B"), paste0("GT3", "_C"),
#'      paste0("Levo", "_A"), paste0("Levo", "_B"), paste0("Levo", "_C"),
#'      paste0("Lome", "_A"), paste0("Lome", "_B"), paste0("Lome", "_C"),
#'      paste0("Novo", "_A"), paste0("Novo", "_B"), paste0("Novo", "_C"),
#'      paste0("PA14wt", "_A"), paste0("PA14wt", "_B"), paste0("PA14wt", "_C")
#'    ),
#'                   sample_group = c(
#'                       rep("AR351", 3), rep("Cipro", 3), rep("Coumermycin", 3),
#'                        rep("GT3", 3), rep("Levo", 3), rep("Lome", 3),
#'                       rep("Novo", 3), rep("PA14wt", 3)
#'                     ),
#'                   stringsAsFactors = F)
#' #'
#' #' # Create raw data object- "XCMSnExp"
#' msraw <- MSnbase::readMSData(files = fp_raw,
#'                               pdata = new("AnnotatedDataFrame", pd),
#'                               mode = "onDisk", msLevel. = 1, centroided. = T)
#'
#' #' # "XcmsExperiment" object
#' msraw_experiment <- readMsExperiment(spectraFiles = fp_raw, sampleData = pd)
#' msraw_experiment
#' #'
#' #'
#' #' # --- PEAK DETECTION ---
#' #'
#' #' # run peak detection
#' xdata_experiment <- xcms::findChromPeaks( msraw_experiment, #msraw,
#'
#'                                param = xcms::CentWaveParam(
#'                                peakwidth = c(10, 25),
#'                                ppm = 10,
#'                                snthresh = 10,
#'                                mzdiff = 0.01,
#'                                prefilter = c(3, 1000),
#'                                noise = 200,
#'                                verboseBetaColumns = TRUE),
#'                                msLevel = 1L)
#'
#' class(xdata_experiment)
#' #'
#' #'
#' #' eics <- chromPeakChromatograms(
#' #'   xdata, peaks = rownames(chromPeaks(xdata))[c(1,2)])
#' #'
#' #' peak_1 <- eics[1]
#' #' peak_2 <- eics[2]
#' #' par(mfrow = c(1, 2))
#' #' plot(peak_1)
#' #' plot(peak_2)
#' #'
#' #
#' beta_metrics <- chromPeakSummary(xdata_experiment, BetaDistributionParam())
#' head(beta_metrics)
#'
#' # after peak picking + (optionally) adjustRtime
#' pdp <- PeakDensityParam(sampleGroups = pd$sample_group,
#'                         minFraction = 0.5, bw = 5)
#'
#' xdata_experiment <- groupChromPeaks(xdata_experiment, param = pdp)
#'
#' # Optionally fill missing peaks
#' xdata_experiment <- fillChromPeaks(xdata_experiment, ChromPeakAreaParam())
#'
#' # Now featureDefinitions() is populated
#' beta_check <- link_feature_CP_id(xdata_experiment)
#'
#' # check if merge back to FT works
#'
#' beta_check_wide <- link_feature_CP_wide(x = xdata_experiment )
#'
#' head(beta_check_wide)
#'
#' #'
#' #'
#' #'
#' #'
#' #' ####################################################
#' #'
#' #'
#' #'
#' #' # Filter MS1 level
#' #' xdata_MS1 <- filterMsLevel(xdata_experiment, msLevel = 1L)
#' #'
#' #'
#' #'
#' #' # Retention time correction alignment
#' #'
#' #' xdata_MS1 <- adjustRtime(xdata_MS1, param = ObiwarpParam())
#' #'
#' #'
#' #' # Peak correspondence across samples
#' #' pdp <- PeakDensityParam(sampleGroups = pd$sample_group,
#' #'                         minFraction = 0.5, bw = 5)
#' #'
#' #' xdata_MS1 <- groupChromPeaks(xdata_MS1, param = pdp)
#' #'
#' #' # Fill missing peaks
#' #' xdata_MS1 <- fillChromPeaks(xdata_MS1, param = ChromPeakAreaParam())
#' #'
#' #'
#' #' # feature definition
#' #' fea_def <- featureDefinitions(xdata_MS1)
#' #'
#' #' # feature intensity value
#' #' fea_int <- featureValues(xdata_MS1, value = "into")
#' #'
#' #'
#' #'
#' #'
#' #'
#' #'
#' #'
#' #'
#' #'
#' #' #-----------------------------------------------#
#' #' # cpc
#' #'
#' #'
#' #'
#' #' # Create raw data object- "XCMSnExp"
#' msraw <- MSnbase::readMSData(files = fp_raw,
#'                              pdata = new("AnnotatedDataFrame", pd),
#'                              mode = "onDisk", msLevel. = 1, centroided. = T)
#'
#'
#'
#' x_nexp <-xcms::findChromPeaks( msraw,
#'
#'                                param = xcms::CentWaveParam(
#'                                  peakwidth = c(10, 25),
#'                                  ppm = 10,
#'                                  snthresh = 10,
#'                                  mzdiff = 0.01,
#'                                  prefilter = c(3, 1000),
#'                                  noise = 200,
#'                                  verboseBetaColumns = TRUE),
#'                                msLevel = 1L)
#'
#' # after peak picking + (optionally) adjustRtime
#' pdp <- PeakDensityParam(sampleGroups = pd$sample_group,
#'                         minFraction = 0.5, bw = 5)
#'
#' x_nexp <- groupChromPeaks(x_nexp, param = pdp)
#'
#' # Optionally fill missing peaks
#' x_nexp <- fillChromPeaks(x_nexp, ChromPeakAreaParam())
#'
#'
#' cpc <- cpc::filter_xcms_peaklist(xd = x_nexp , return_type = "cpc",
#'                                  param = cpc::cpcProcParam(min_sn = 10,
#'                                                            min_pts = 10,
#'                                                            min_intensity = 2000))
#'
#'
#'
#' cpcPeaktable <- cpc::cpt(cpc)
#'
#' str(cpcPeaktable)
#' head(cpcPeaktable)
#'
#'
#' dim(cpcPeaktable)
#'
#'
#' # check cpc wide
#'
#' # cpc_wide <- link_feature_CPC_wide(x_nexp, cpc)
#' #
#' # dim(cpc_wide)
#'
#' #' ############################################################################
#' #'
#' #' # get the original xcms object using getOriginalXCMS()
#' #' xdOrig <- cpc::getOriginalXCMS(cpc)
#' #'
#' #' xdOrig
#' #'
#' #' # get the original peak table
#' #' xcmsPeaktableOrig <- data.frame(xcms::chromPeaks(xdOrig))
#' #'
#' #' # get the filtered xcms object using getFilteredXCMS()
#' #' xdFilt <- cpc::getFilteredXCMS(cpc)
#' #'
#' #' xdFilt
#' #'
#' #' # get the filtered peak table
#' #' xcmsPeaktableFilt <- data.frame(xcms::chromPeaks(xdFilt))
#' #'
#' #' # get the filtered cpc peak table
#' #' cpcPeaktableFilt <- cpcPeaktable[row.names(xcmsPeaktableFilt), ]
#' #'
#' #'
#' #'
#' #' rbind(Original = table(xcmsPeaktableOrig[, "sample"]),
#' #'       Filtered = table(xcmsPeaktableFilt[, "sample"]),
#' #'       PercentKept = round(100*table(xcms::chromPeaks(xdFilt)[, "sample"]) /
#' #'                             table(xcms::chromPeaks(xdOrig)[, "sample"]), 1))
#' #'
#' #'
#' #'
#' outcomes <- cpc::getFilterOutcomes(cpc)
#'
#' head(outcomes)
#' #'
#' #'
#' #'
#' #'
#' #' # setup a list for the Venn diagram
#' #' outcomesList <- lapply(outcomes, function(x) which(!x))
#' #'
#' #' # output a Venn diagram for the used filters
#' #' ggvenn::ggvenn(data = outcomesList[c(1, 2, 3, 6)], show_percentage = FALSE,
#' #'                stroke_size = 0.5, set_name_size = 4, text_size = 2.8,
#' #'                fill_color = c("#868686FF", "#0073C2FF", "#CD534CFF",
#' #'                               "#EFC000FF", "#073b4cFF"))
#' #'
#' #'
#' #'
#' #' # set a seed for sample() to make it reproducible
#' #' seed <- 12345
#' #'
#' #' # as an example we will select 2 peaks from each filter outcome
#' #' # (1) Not detected as a peak,
#' #' # (2) Too low signal-to-noise,
#' #' # (3) Too low intensity (only 1 peak here)
#' #' # (4) Too narrow
#' #' set.seed(seed);
#' randomRemoved <- list(
#'   nd = sample(which(apply(outcomes, 1, function(x) { # (1)
#'     return(!x[1] &
#'              x[2] &
#'              x[3] &
#'              x[4] &
#'              x[5] &
#'              x[6])
#'   })), 2),
#'   lowsn = sample(which(apply(outcomes, 1, function(x) { # (2)
#'     return(x[1] &
#'              !x[2] &
#'              x[3] &
#'              x[4] &
#'              x[5] &
#'              x[6])
#'   })), 2),
#'   lowint = which(apply(outcomes, 1, function(x) { # (3)
#'     return(x[1] &
#'              !x[2] &
#'              x[3] &
#'              x[4] &
#'              x[5] &
#'              !x[6])
#'   })),
#'   toonarrow = sample(which(apply(outcomes, 1, function(x) { # (4)
#'     return(x[1] &
#'              x[2] &
#'              !x[3] &
#'              x[4] &
#'              x[5] &
#'              x[6])
#'   })), 2)
#' )
#' #'
#' #' randomRemoved
#' #' cpc::plotPeaks(cpc, peakIdx = randomRemoved[[2]][[1]])
#' #'
#' cpc::plotPeaks(cpc, peakIdx = randomRemoved[[2]][[1]])
#'
#'
#'
#'
#'
#' cpc::plotPeaks(cpc, peakIdx = 1)
#'
#'
#'
#'
#' cpc::plotPeaks(cpc, peakIdx = 189)
#'
#'
#'
