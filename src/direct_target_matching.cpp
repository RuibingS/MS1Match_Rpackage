#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <map>
#include <limits>
#include "utils_check.h"
#include "utils_sort.h"
using namespace Rcpp;


// HMDB target list csv

// [[Rcpp::export]]
Rcpp::List PImatch_target_rcpp(Rcpp::DataFrame FT,
                               Rcpp::DataFrame Comp_data,
                               Rcpp::CharacterVector PIon,
                               double diff_mz_thr,
                               std::string polarity,
                               int ppm) {
  int nrow_FT = FT.nrow();
  int nrow_Comp = Comp_data.nrow();

  NumericVector FT_mz = get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
  CharacterVector feature_name = get_required_column<CharacterVector>(FT, "feature_name");
  NumericVector FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention time");

  Rcpp::CharacterVector Comp_name = Comp_data["Name"];
  Rcpp::NumericVector Comp_M = get_column_with_fallback<Rcpp::NumericVector>(
    Comp_data, "exactMass", "MONO_MASS", "Compound exact mass"
  );



  Rcpp::NumericVector comp_mz;

  // Choose correct column based on polarity
  if (polarity == "pos") {
    comp_mz = Comp_data["[M+H]+"];
  } else if (polarity == "neg") {
    comp_mz = Comp_data["[M-H]-"];
  } else {
    stop("Invalid polarity: must be 'pos' or 'neg'");
  }


  String PIon_name = PIon[0];

  Rcpp::List output(nrow_FT);

  for (size_t  i = 0; i < nrow_FT; i++) {
    CharacterVector feaNames;
    CharacterVector CompNames;
    CharacterVector PInames;
    NumericVector fea_mz;
    NumericVector M_vec;

    for (size_t  j = 0; j < nrow_Comp; j++) {
      double mz_diff = std::fabs(FT_mz[i] - comp_mz[j]);
      double ppm_tol = comp_mz[j] * ppm * 1e-6;

      if (mz_diff <= ppm_tol || mz_diff <= diff_mz_thr) {
        feaNames.push_back(feature_name[i]);
        fea_mz.push_back(FT_mz[i]);
        CompNames.push_back(Comp_name[j]);
        PInames.push_back(PIon_name);
        M_vec.push_back(Comp_M[j]);
      }
    }

    output[i] = List::create(
      _["Feature_name"] = feaNames,
      _["Feature_mz"] = fea_mz,
      _["Comp_name"] = CompNames,
      _["PI_name"] = PInames,
      _["M"] = M_vec
    );
  }

  return output;
}
