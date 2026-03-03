
#' Import PAMD Metabolite Data and Calculate Adducts
#'
#' @description Reads a CSV file, removes all-NA columns, computes exact mass from formula,
#' and calculates adduct m/z values using MetaboCoreUtils.
#'
#' @param target_list_path Path to the target list CSV file containing metabolite data.
#' @importFrom enviPat isopattern
#' @importFrom MetaboCoreUtils mass2mz
#' @importFrom utils data
#' @return A data.frame with original data, exact masses, and adduct m/z.
#' @export
PAMD_target_list_import <- function(target_list_path, IonPolarity){
  data(isotopes, package = "enviPat")
  metabo_data <- read.csv(target_list_path)

  #remove columns where all vales are NA
  metabo_data <- metabo_data[,colSums(is.na(metabo_data))<nrow(metabo_data)]

  # calulate exact mass from formula
  exactMass <- unname(MetaboCoreUtils::calculateMass(metabo_data$Formula))




  #calculate primary Adducts "[M+H]+"
 metabo_data_add_pos <- MetaboCoreUtils::formula2mz(metabo_data$Formula, adduct=c("[M+H]+"))


  #calculate primary Adducts [M-H]-

 metabo_data_add_neg <- MetaboCoreUtils::formula2mz(metabo_data$Formula, adduct=c("[M-H]-"))


  return(cbind(metabo_data, exactMass, "[M+H]+" = metabo_data_add_pos, "[M-H]-"= metabo_data_add_neg))
}




#' Import HUMD Metabolite Data and Calculate Adducts
#'
#' Reads a CSV file, removes all-NA columns, computes exact mass from formula,
#' and calculates adduct m/z values using MetaboCoreUtils.
#'
#' @param target_list_path Path to the target list CSV file containing metabolite data.
#' @param IonPolarity "pos" or "neg" to calculate both PI adducts.
#' @importFrom enviPat isopattern
#' @importFrom MetaboCoreUtils mass2mz
#' @importFrom utils data
#' @return A data.frame with original data, exact masses, and adduct m/z.
#' @export

#'
HMDB_target_list_import <- function(target_list_path, IonPolarity){

  metabo_data <- read.csv(target_list_path)

  #remove columns where all vales are NA
  metabo_data <- metabo_data[,colSums(is.na(metabo_data))<nrow(metabo_data)]

  # calulate exact mass from formula
  exactMass <- metabo_data$MONO_MASS


  #calculate primary Adducts "[M+H]+"

  metabo_data_add_pos <- MetaboCoreUtils::mass2mz(exactMass, adduct=c("[M+H]+"))


  #calculate primary Adducts [M-H]-

  metabo_data_add_neg <- MetaboCoreUtils::mass2mz(exactMass, adduct=c("[M-H]-"))



  return(cbind(metabo_data, exactMass, "[M+H]+" = metabo_data_add_pos, "[M-H]-" = metabo_data_add_neg))
}
