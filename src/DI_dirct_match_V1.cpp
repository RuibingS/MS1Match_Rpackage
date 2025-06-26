#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <map>
#include <limits>
using namespace Rcpp;


// HMDB target list csv

// duplicated: // [[Rcpp::export]]
// Rcpp::List PImatch_hmdb_fun(Rcpp::DataFrame FT,
//                             Rcpp::DataFrame Comp_data,
//                             Rcpp::CharacterVector PIon,
//                             double diff_mz_thr,
//                             int ppm) {
//   int nrow_FT = FT.nrow();
//   int nrow_Comp = Comp_data.nrow();
//
//   Rcpp::CharacterVector feature_name = FT["feature_name"];
//   Rcpp::NumericVector FT_mz = FT["mz"];
//
//   Rcpp::CharacterVector Comp_name = Comp_data["NAME"];
//   Rcpp::NumericVector Comp_M = Comp_data["MONO_MASS"];
//
//   // Assume PIon is e.g., "[M+H]+" — only using the first entry for now
//   String PIon_name = PIon[0];
//
//   // Calculate PI m/z for each compound
//   Rcpp::NumericVector comp_mz = Comp_M + 1.007276;
//
//   Rcpp::List output(nrow_FT);
//
//   for (int i = 0; i < nrow_FT; i++) {
//     CharacterVector feaNames;
//     CharacterVector CompNames;
//     CharacterVector PInames;
//     NumericVector fea_mz;
//     NumericVector M_vec;
//
//     for (int j = 0; j < nrow_Comp; j++) {
//       double mz_diff = std::fabs(FT_mz[i] - comp_mz[j]);
//       double ppm_tol = comp_mz[j] * ppm * 1e-6;
//
//       if (mz_diff <= ppm_tol || mz_diff <= diff_mz_thr) {
//         feaNames.push_back(feature_name[i]);
//         fea_mz.push_back(FT_mz[i]);
//         CompNames.push_back(Comp_name[j]);
//         PInames.push_back(PIon_name);
//         M_vec.push_back(Comp_M[j]);
//       }
//     }
//
//     output[i] = List::create(
//       _["Feature_name"] = feaNames,
//       _["Feature_mz"] = fea_mz,
//       _["Comp_name"] = CompNames,
//       _["PI_name"] = PInames,
//       _["M"] = M_vec
//     );
//   }
//
//   return output;
// }


// [[Rcpp::export]]
Rcpp::List DI_dirct_match_fun(Rcpp::DataFrame FT,
                           Rcpp::List compounds_list,
                           double diff_mz_thr,
                           int ppm) {

  int nrow_FT = FT.nrow();
  int nrow_Comp = compounds_list.length();

  Rcpp::CharacterVector feature_name = FT["feature_name"];
  Rcpp::NumericVector FT_mz = FT["mz"];
  Rcpp::NumericVector FT_intensity = FT.containsElementNamed("int") ? FT["int"] : NumericVector(nrow_FT, NA_REAL);

  Rcpp::List output_list;
  std::map<std::string, std::pair<double, Rcpp::List>> closest_matches;

  for (int i = 0; i < nrow_Comp; i++) {
    List compound = compounds_list[i];

    if (!compound.containsElementNamed("Name") ||
        !compound.containsElementNamed("Formula") ||
        !compound.containsElementNamed("ExactMass") ||
        !compound.containsElementNamed("Peaks")) {
        Rcpp::Rcout << "Skipping compound at index " << i << " due to missing required fields.\n";
      continue;
    }

    std::string compound_name = compound["Name"];
    double ExactMass_val = compound["ExactMass"];
    double PI_H_val = ExactMass_val + 1.007276;
    double PI_Na_val = ExactMass_val + 22.989218;

    Rcpp::DataFrame Peaks = Rcpp::as<Rcpp::DataFrame>(compound["Peaks"]);
    NumericVector peak_mz = Peaks["mz"];
    NumericVector peak_rel_int = Peaks["rel_intensity"];
    int total_peaks = peak_mz.size();

    for (int j = 0; j < nrow_FT; j++) {
      double feature_mz = FT_mz[j];

      double diff_mz_H = std::fabs(feature_mz - PI_H_val);
      double diff_mz_Na = std::fabs(feature_mz - PI_Na_val);
      double ppm_thr_H = PI_H_val * ppm * 1e-6;
      double ppm_thr_Na = PI_Na_val * ppm * 1e-6;

      bool match_H = (diff_mz_H <= diff_mz_thr || diff_mz_H <= ppm_thr_H);
      bool match_Na = (diff_mz_Na <= diff_mz_thr || diff_mz_Na <= ppm_thr_Na);

      if (match_H || match_Na) {
        std::string ion_type = match_H ? "[M+H]+" : "[M+Na]+";

        int peak_matches = 0;
        std::vector<std::string> matched_peak_features;
        Rcpp::NumericVector matched_feature_mz;
        Rcpp::NumericVector matched_feature_int;
        Rcpp::NumericVector match_peak_mz;
        Rcpp::NumericVector match_peak_int;

        for (int k = 0; k < nrow_FT; k++) {
          double neighbor_mz = FT_mz[k];

          for (int p = 0; p < peak_mz.size(); p++) {
            if (std::fabs(neighbor_mz - peak_mz[p]) <= diff_mz_thr) {
              peak_matches++;
              match_peak_mz.push_back(peak_mz[p]);
              match_peak_int.push_back(peak_rel_int[p]);
              matched_peak_features.push_back(std::string(feature_name[k]));
              matched_feature_mz.push_back(FT_mz[k]);
              matched_feature_int.push_back(FT_intensity[k]);
              break;
            }
          }
        }

        double match_score = peak_matches > 0 && total_peaks > 0
        ? static_cast<double>(peak_matches) / total_peaks
        : NA_REAL;

        Rcpp::CharacterVector matched_features_vec;
        if (matched_peak_features.size() > 0) {
          matched_features_vec = wrap(matched_peak_features);
        } else {
          matched_features_vec = Rcpp::CharacterVector::create(NA_STRING);
        }

        Rcpp::NumericVector matched_feature_mz_vec =
          matched_feature_mz.size() > 0 ? matched_feature_mz : Rcpp::NumericVector::create(NA_REAL);

       Rcpp::NumericVector matched_feature_int_vec =
          matched_feature_int.size() > 0 ? matched_feature_int : Rcpp::NumericVector::create(NA_REAL);

        Rcpp::List feature_match = Rcpp::List::create(
          _["feature_name"] = std::string(feature_name[j]),
          _["feature_mz"] = feature_mz,
          _["feature_intensity"] = FT_intensity[j],
                                               _["Compound_name"] = compound_name,
                                               _["Compound_peak_mz"] = match_peak_mz.size() > 0 ? match_peak_mz : Rcpp::NumericVector::create(NA_REAL),
                                               _["Compound_peak_int"] = match_peak_int.size() > 0 ? match_peak_int : Rcpp::NumericVector::create(NA_REAL),
                                               _["Ion_type"] = ion_type,
                                               _["Matched_peaks"] = peak_matches,
                                               _["Match_peak_score"] = match_score,
                                               _["Total_library_peaks"] = total_peaks,
                                               _["Matched_features_to_peaks"] = matched_features_vec,
                                               _["Matched_features_mz"] = matched_feature_mz_vec,
                                               _["Matched_features_int"] = matched_feature_int_vec
        );


        double mz_diff = std::min(diff_mz_H, diff_mz_Na);
        std::string current_feature_name = std::string(feature_name[j]);

        if (closest_matches.find(current_feature_name) == closest_matches.end() ||
            mz_diff < closest_matches[current_feature_name].first) {
          closest_matches[current_feature_name] = {mz_diff, feature_match};
        }
      }
    }
  }

  for (const auto &entry : closest_matches) {
    output_list.push_back(entry.second.second);
  }

  return output_list;
}
