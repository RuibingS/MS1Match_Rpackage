// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include "utils_check.h"
#include "utils_sort.h"
#include "utils_match.h"

using namespace Rcpp;


// --------------------------------------------------------
// Helper: extract RT from compound
// --------------------------------------------------------
inline bool get_rt(const List& compound, double& rt_val) {
  if (!compound.containsElementNamed("RetTime")) return false;
  try {
    rt_val = std::stod(as<std::string>(compound["RetTime"]));
    return true;
  } catch (...) {
    return false;
  }
}


// --------------------------------------------------------
// Helper: convert string to lowercase
// --------------------------------------------------------
inline std::string lower_str(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  return s;
}


// --------------------------------------------------------
// Compute primary ion m/z from exact mass, ion mode, and adduct
// --------------------------------------------------------
inline double compute_ion_mz(double exact_mass,
                             const std::string& ion_mode_raw,
                             const std::string& adduct_raw) {

  std::string ion_mode = lower_str(ion_mode_raw);
  std::string adduct   = lower_str(adduct_raw);

  const double proton    = 1.007276;
  const double sodium    = 22.989218;
  const double potassium = 38.963158;

  // POSITIVE MODE
  if (ion_mode == "positive" || ion_mode == "pos") {
    if (adduct == "m+h" || adduct == "[m+h]+")
      return exact_mass + proton;
    if (adduct == "m+na" || adduct == "[m+na]+")
      return exact_mass + sodium;
    if (adduct == "m+k" || adduct == "[m+k]+")
      return exact_mass + potassium;

    // default POS adduct
    return exact_mass + proton;
  }

  // NEGATIVE MODE
  if (ion_mode == "negative" || ion_mode == "neg") {
    if (adduct == "m-h" || adduct == "[m-h]-")
      return exact_mass - proton;

    // default NEG adduct
    return exact_mass - proton;
  }

  // Unknown ion mode → assume [M+H]+
  return exact_mass + proton;
}



// --------------------------------------------------------
// Main
// --------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List direct_match_MS1Lib_cpp(Rcpp::DataFrame FT,
                                   Rcpp::List compounds_list,
                                   double diff_mz_thr,
                                   double diff_rt_thr,
                                   int ppm) {

  int nrow_FT = FT.nrow();

  // Extract feature table fields
  NumericVector FT_mz =
    get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");

  NumericVector FT_int =
    get_column_with_fallback<NumericVector>(FT, "intensity_mean",
                                            "mean_intensity", "intensity");

  CharacterVector feature_name =
    get_required_column<CharacterVector>(FT, "feature_name");

  // Optional RT
  NumericVector FT_rt;
  bool has_rt_FT = false;
  try {
    FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention_time");
    has_rt_FT = (FT_rt.size() == FT_mz.size());
  } catch (...) {
    has_rt_FT = false;
  }

  // Sort feature list by m/z
  if (!std::is_sorted(FT_mz.begin(), FT_mz.end())) {
    std::vector<int> idx = compute_sorted_index(FT_mz);
    FT_mz = reorder_vector(FT_mz, idx);
    FT_int = reorder_vector(FT_int, idx);
    feature_name = reorder_vector(feature_name, idx);
    if (has_rt_FT) FT_rt = reorder_vector(FT_rt, idx);
  }

  Rcpp::List output_list;
  int group_index = 1;


  for (int i = 0; i < compounds_list.size(); ++i) {

    List compound = compounds_list[i];

    if (!compound.containsElementNamed("Name") ||
        !compound.containsElementNamed("Peaks") ||
        !compound.containsElementNamed("ExactMass"))
        continue;

    std::string compound_name = as<std::string>(compound["Name"]);
    double exact_mass = as<double>(compound["ExactMass"]);

    std::string ion_mode = "POSITIVE";
    if (compound.containsElementNamed("Ion_mode"))
      ion_mode = as<std::string>(compound["Ion_mode"]);
    else if (compound.containsElementNamed("IonPolarity"))
      ion_mode = as<std::string>(compound["IonPolarity"]);

    std::string adduct = "M+H";
    if (compound.containsElementNamed("Adduct"))
      adduct = as<std::string>(compound["Adduct"]);

    double mz_pi = compute_ion_mz(exact_mass, ion_mode, adduct);

    double mz_tol = std::max(diff_mz_thr, mz_pi * ppm * 1e-6);

    double compound_rt = 0.0;
    bool has_rt_comp = get_rt(compound, compound_rt);

    // Extract peaks
    DataFrame Peaks = as<DataFrame>(compound["Peaks"]);
    NumericVector peak_mz  = Peaks["mz"];
    NumericVector peak_int = Peaks["rel_intensity"];
    int total_peaks = peak_mz.size();

    // Sort peak list
    if (!std::is_sorted(peak_mz.begin(), peak_mz.end())) {
      std::vector<int> idx = compute_sorted_index(peak_mz);
      peak_mz = reorder_vector(peak_mz, idx);
      peak_int = reorder_vector(peak_int, idx);
    }


    // --------------------------------------------------------
    //  matching the PI ion
    // --------------------------------------------------------
    for (int j = 0; j < nrow_FT; ++j) {

      double fmz = FT_mz[j];
      double frt = has_rt_FT ? FT_rt[j] : NA_REAL;

      // m/z check
      if (std::fabs(fmz - mz_pi) > mz_tol) continue;

      // RT check
      if (has_rt_comp && has_rt_FT &&
          std::fabs(frt - compound_rt) > diff_rt_thr)
        continue;


      // --------------------------------------------------------
      // within RT window matching
      // --------------------------------------------------------
      std::vector<int> candidates;
      for (int k = 0; k < nrow_FT; ++k) {

        double rt_k = has_rt_FT ? FT_rt[k] : NA_REAL;

        if (has_rt_comp && has_rt_FT &&
            std::fabs(rt_k - compound_rt) > diff_rt_thr)
          continue;

        candidates.push_back(k);
      }

      NumericVector sub_mz(candidates.size());
      NumericVector sub_rt(candidates.size());
      NumericVector sub_int(candidates.size());
      CharacterVector sub_name(candidates.size());

      for (size_t idx = 0; idx < candidates.size(); ++idx) {
        int k = candidates[idx];
        sub_mz[idx] = FT_mz[k];
        sub_rt[idx] = has_rt_FT ? FT_rt[k] : NA_REAL;
        sub_int[idx] = FT_int[k];
        sub_name[idx] = feature_name[k];
      }


      auto matches = find_one_to_one_matches(sub_mz, peak_mz, diff_mz_thr);
      int matched_peaks = matches.size();

      std::vector<std::string> matched_peak_features;
      NumericVector matched_feature_mz, matched_feature_rt, matched_feature_int;
      NumericVector match_peak_mz, match_peak_int;

      for (const auto& match : matches) {

        int feat_idx = std::get<0>(match);
        int peak_idx = std::get<1>(match);

        matched_peak_features.push_back(as<std::string>(sub_name[feat_idx]));
        matched_feature_mz.push_back(sub_mz[feat_idx]);
        matched_feature_rt.push_back(sub_rt[feat_idx]);
        matched_feature_int.push_back(sub_int[feat_idx]);

        match_peak_mz.push_back(peak_mz[peak_idx]);
        match_peak_int.push_back(peak_int[peak_idx]);
      }

      double score = matched_peaks > 0 ?
      static_cast<double>(matched_peaks) / total_peaks :
        NA_REAL;

      std::string lib_group_id = "Lib_group_" + std::to_string(group_index++);


      CharacterVector matched_peak_features_char;
      if (matched_peak_features.empty())
        matched_peak_features_char = CharacterVector::create(NA_STRING);
      else
        matched_peak_features_char = wrap(matched_peak_features);

      output_list.push_back(List::create(
          _["feature_name"] = as<std::string>(feature_name[j]),
          _["Lib_group"] = lib_group_id,
          _["feature_mz"] = fmz,
          _["feature_rt"] = frt,

          _["Compound_name"] = compound_name,
          _["Compound_M"] = exact_mass,
          _["Compound_rt"] = has_rt_comp ? compound_rt : NA_REAL,

          _["Compound_peak_mz"]  = match_peak_mz.size() ? match_peak_mz  : NumericVector::create(NA_REAL),
          _["Compound_peak_int"] = match_peak_int.size() ? match_peak_int : NumericVector::create(NA_REAL),

          _["Compound_peak_mz_full"] = peak_mz,
          _["Compound_peak_int_full"] = peak_int,

          _["Total_library_peaks"] = total_peaks,
          _["Ion_type"] = adduct,
          _["Matched_peaks"] = matched_peaks,
          _["Match_peak_score"] = score,

          _["Matched_features_to_peaks"] = matched_peak_features_char,
          _["Matched_features_mz"] = matched_feature_mz.size() ? matched_feature_mz : NumericVector::create(NA_REAL),
          _["Matched_features_rt"] = matched_feature_rt.size() ? matched_feature_rt : NumericVector::create(NA_REAL),
          _["Matched_features_int"] = matched_feature_int.size() ? matched_feature_int : NumericVector::create(NA_REAL)
      ));
    }
  }

  return output_list;
}
















// [[Rcpp::export]]
Rcpp::List direct_match_MS2Lib_cpp(Rcpp::DataFrame FT,
                                   Rcpp::List compounds_list,
                                   double diff_mz_thr,
                                   double diff_rt_thr,
                                   int ppm) {
  int nrow_FT = FT.nrow();
  int n_compounds = compounds_list.length();

  Rcpp::NumericVector FT_mz = get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
  Rcpp::CharacterVector feature_name = get_required_column<CharacterVector>(FT, "feature_name");
  Rcpp::NumericVector FT_int = get_column_with_fallback<NumericVector>(FT, "intensity_mean", "mean_intensity", "intensity");

  Rcpp::NumericVector FT_rt;
  bool has_rt = false;
  try {
    FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention_time");
    has_rt = FT_rt.size() == FT_mz.size();
  } catch (...) {
    has_rt = false;
  }

  if (!std::is_sorted(FT_mz.begin(), FT_mz.end())) {
    std::vector<int> idx = compute_sorted_index(FT_mz);
    FT_mz = reorder_vector(FT_mz, idx);
    if (has_rt) FT_rt = reorder_vector(FT_rt, idx);
    feature_name = reorder_vector(feature_name, idx);
    FT_int = reorder_vector(FT_int, idx);
  }

  Rcpp::List output_list;
  int group_index = 1;

  for (int i = 0; i < n_compounds; ++i) {
    Rcpp::List compound = compounds_list[i];
    if (!compound.containsElementNamed("Name") || !compound.containsElementNamed("Peaks")) continue;

    std::string compound_name = Rcpp::as<std::string>(compound["Name"]);
    double exact_mass = as<double>(compound["ExactMass"]);

    double pre_mz = NA_REAL;
    try {
      pre_mz = std::stod(Rcpp::as<std::string>(compound["PreIon"]));
    } catch (...) {
      continue;
    }

    double compound_rt = 0.0;
    bool has_rt_comp = get_rt(compound, compound_rt);
    double mz_tol = std::max(diff_mz_thr, pre_mz * ppm * 1e-6);

    Rcpp::DataFrame Peaks = Rcpp::as<Rcpp::DataFrame>(compound["Peaks"]);
    Rcpp::NumericVector peak_mz = Peaks["mz"];
    Rcpp::NumericVector peak_int = Peaks["rel_intensity"];
    int total_peaks = peak_mz.size();

    if (!std::is_sorted(peak_mz.begin(), peak_mz.end())) {
      std::vector<int> idx = compute_sorted_index(peak_mz);
      peak_mz = reorder_vector(peak_mz, idx);
      peak_int = reorder_vector(peak_int, idx);
    }

    for (int j = 0; j < nrow_FT; ++j) {
      double fmz = FT_mz[j];
      double frt = has_rt ? FT_rt[j] : NA_REAL;

      if (std::fabs(fmz - pre_mz) > mz_tol) continue;
      if (has_rt && has_rt_comp && std::fabs(frt - compound_rt) > diff_rt_thr) continue;

      std::vector<std::string> matched_peak_features;
      Rcpp::NumericVector matched_feature_mz, matched_feature_int, matched_feature_rt;
      Rcpp::NumericVector match_peak_mz, match_peak_int;
      int matched_peaks = 0;

      std::vector<bool> peak_used(total_peaks, false);
      std::vector<bool> ft_used(nrow_FT, false);

      for (int p = 0; p < total_peaks; ++p) {
        double lib_mz = peak_mz[p];
        double best_diff = diff_mz_thr + 1.0;
        int best_k = -1;

        for (int k = 0; k < nrow_FT; ++k) {
          if (ft_used[k]) continue;

          double mz_k = FT_mz[k];
          double rt_k = has_rt ? FT_rt[k] : NA_REAL;

          if (has_rt && has_rt_comp && std::fabs(rt_k - compound_rt) > diff_rt_thr) continue;

          double mz_diff = std::fabs(lib_mz - mz_k);
          if (mz_diff <= diff_mz_thr && mz_diff < best_diff) {
            best_diff = mz_diff;
            best_k = k;
          }
        }

        if (best_k >= 0) {
          peak_used[p] = true;
          ft_used[best_k] = true;

          matched_peaks++;
          match_peak_mz.push_back(lib_mz);
          match_peak_int.push_back(peak_int[p]);
          matched_peak_features.push_back(Rcpp::as<std::string>(feature_name[best_k]));
          matched_feature_mz.push_back(FT_mz[best_k]);
          matched_feature_rt.push_back(has_rt ? FT_rt[best_k] : NA_REAL);
          matched_feature_int.push_back(FT_int[best_k]);
        }
      }

      double score = matched_peaks > 0 ? (double)matched_peaks / total_peaks : NA_REAL;
      std::string lib_group_id = "Lib_group_" + std::to_string(group_index++);

      Rcpp::CharacterVector matched_peak_features_char;
      if (!matched_peak_features.empty()) {
        matched_peak_features_char = Rcpp::wrap(matched_peak_features);
      } else {
        matched_peak_features_char = Rcpp::CharacterVector::create(NA_STRING);
      }

      output_list.push_back(List::create(
          _["feature_name"] = Rcpp::as<std::string>(feature_name[j]),
          _["Lib_group"] = lib_group_id,
          _["feature_mz"] = fmz,
          _["feature_rt"] = frt,
          _["Compound_name"] = compound_name,
          _["Compound_M"] = exact_mass,
          _["Compound_rt"] = has_rt_comp ? compound_rt : NA_REAL,
          _["Precursor_mz"] = pre_mz,
          _["Compound_peak_mz"] = match_peak_mz.size() ? match_peak_mz : NumericVector::create(NA_REAL),
          _["Compound_peak_int"] = match_peak_int.size() ? match_peak_int : NumericVector::create(NA_REAL),
          _["Compound_peak_mz_full"] = peak_mz,
          _["Compound_peak_int_full"] = peak_int,
          _["Matched_peaks"] = matched_peaks,
          _["Match_peak_score"] = score,
          _["Total_library_peaks"] = total_peaks,
          _["Ion_type"] = "Precursor",
          _["Matched_features_to_peaks"] = matched_peak_features_char,
          _["Matched_features_mz"] = matched_feature_mz.size() ? matched_feature_mz : NumericVector::create(NA_REAL),
          _["Matched_features_rt"] = matched_feature_rt.size() ? matched_feature_rt : NumericVector::create(NA_REAL),
          _["Matched_features_int"] = matched_feature_int.size() ? matched_feature_int : NumericVector::create(NA_REAL)
      ));
    }
  }

  return output_list;
}



// Rcpp::List direct_match_MS1Lib_cpp(Rcpp::DataFrame FT,
//                                    Rcpp::List compounds_list,
//                                    double diff_mz_thr,
//                                    double diff_rt_thr,
//                                    int ppm) {
//   int nrow_FT = FT.nrow();
//   NumericVector FT_mz = get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
//   NumericVector FT_int = get_column_with_fallback<NumericVector>(FT, "intensity_mean", "mean_intensity", "intensity");
//   CharacterVector feature_name = get_required_column<CharacterVector>(FT, "feature_name");
//
//   NumericVector FT_rt;
//   bool has_rt_FT = false;
//   try {
//     FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention_time");
//     has_rt_FT = FT_rt.size() == FT_mz.size();
//   } catch (...) {
//     has_rt_FT = false;
//   }
//
//   if (!std::is_sorted(FT_mz.begin(), FT_mz.end())) {
//     std::vector<int> idx = compute_sorted_index(FT_mz);
//     FT_mz = reorder_vector(FT_mz, idx);
//     FT_int = reorder_vector(FT_int, idx);
//     feature_name = reorder_vector(feature_name, idx);
//     if (has_rt_FT) FT_rt = reorder_vector(FT_rt, idx);
//   }
//
//   Rcpp::List output_list;
//   int group_index = 1;
//
//   for (int i = 0; i < compounds_list.size(); ++i) {
//     List compound = compounds_list[i];
//
//     if (!compound.containsElementNamed("Name") ||
//         !compound.containsElementNamed("Peaks") ||
//         !compound.containsElementNamed("ExactMass")) continue;
//
//         std::string compound_name = as<std::string>(compound["Name"]);
//         double exact_mass = compound["ExactMass"];
//         double mz_pi = exact_mass + 1.007276;
//         double mz_tol = std::max(diff_mz_thr, mz_pi * ppm * 1e-6);
//         double compound_rt = 0.0;
//         bool has_rt_comp = get_rt(compound, compound_rt);
//
//         DataFrame Peaks = as<DataFrame>(compound["Peaks"]);
//         NumericVector peak_mz = Peaks["mz"];
//         NumericVector peak_int = Peaks["rel_intensity"];
//         int total_peaks = peak_mz.size();
//
//         if (!std::is_sorted(peak_mz.begin(), peak_mz.end())) {
//           std::vector<int> idx = compute_sorted_index(peak_mz);
//           peak_mz = reorder_vector(peak_mz, idx);
//           peak_int = reorder_vector(peak_int, idx);
//         }
//
//         for (int j = 0; j < nrow_FT; ++j) {
//           double fmz = FT_mz[j];
//           double frt = has_rt_FT ? FT_rt[j] : NA_REAL;
//           if (std::fabs(fmz - mz_pi) > mz_tol) continue;
//           if (has_rt_FT && has_rt_comp && std::fabs(frt - compound_rt) > diff_rt_thr) continue;
//
//           int matched_peaks = 0;
//           std::vector<std::string> matched_peak_features;
//           NumericVector matched_feature_mz, matched_feature_rt, matched_feature_int;
//           NumericVector match_peak_mz, match_peak_int;
//
//           for (int k = 0; k < nrow_FT; ++k) {
//             double mz_k = FT_mz[k];
//             double rt_k = has_rt_FT ? FT_rt[k] : NA_REAL;
//             if (has_rt_FT && has_rt_comp && std::fabs(rt_k - compound_rt) > diff_rt_thr) continue;
//
//             for (int p = 0; p < total_peaks; ++p) {
//               if (std::fabs(mz_k - peak_mz[p]) <= diff_mz_thr) {
//                 matched_peaks++;
//                 matched_peak_features.push_back(as<std::string>(feature_name[k]));
//                 matched_feature_mz.push_back(mz_k);
//                 matched_feature_rt.push_back(rt_k);
//                 matched_feature_int.push_back(FT_int[k]);
//                 match_peak_mz.push_back(peak_mz[p]);
//                 match_peak_int.push_back(peak_int[p]);
//                 break;
//               }
//             }
//           }
//
//           double score = matched_peaks > 0 ? (double)matched_peaks / total_peaks : NA_REAL;
//           std::string lib_group_id = "Lib_group_" + std::to_string(group_index++);
//
//           CharacterVector matched_peak_features_char;
//           if (!matched_peak_features.empty()) {
//             matched_peak_features_char = wrap(matched_peak_features);
//           } else {
//             matched_peak_features_char = CharacterVector::create(NA_STRING);
//           }
//
          // output_list.push_back(List::create(
          //     _["feature_name"] = as<std::string>(feature_name[j]),
          //     _["Lib_group"] = lib_group_id,
          //     _["feature_mz"] = fmz,
          //     _["feature_rt"] = frt,
          //     _["Compound_name"] = compound_name,
          //     _["Compound_rt"] = has_rt_comp ? compound_rt : NA_REAL,
          //     _["Compound_peak_mz"] = match_peak_mz.size() ? match_peak_mz : NumericVector::create(NA_REAL),
          //     _["Compound_peak_int"] = match_peak_int.size() ? match_peak_int : NumericVector::create(NA_REAL),
          //     _["Total_library_peaks"] = total_peaks,
          //     _["Ion_type"] = "[M+H]+",
          //     _["Matched_peaks"] = matched_peaks,
          //     _["Match_peak_score"] = score,
          //     _["Matched_features_to_peaks"] = matched_peak_features_char,
          //     _["Matched_features_mz"] = matched_feature_mz.size() ? matched_feature_mz : NumericVector::create(NA_REAL),
          //     _["Matched_features_rt"] = matched_feature_rt.size() ? matched_feature_rt : NumericVector::create(NA_REAL),
          //     _["Matched_features_int"] = matched_feature_int.size() ? matched_feature_int : NumericVector::create(NA_REAL)
          // ));
//         }
//   }
//
//   return output_list;
// }





// Rcpp::List direct_match_MS2Lib_cpp(Rcpp::DataFrame FT,
//                                    Rcpp::List compounds_list,
//                                    double diff_mz_thr,
//                                    double diff_rt_thr,
//                                    int ppm) {
//   int nrow_FT = FT.nrow();
//   int n_compounds = compounds_list.length();
//
//   Rcpp::NumericVector FT_mz = get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
//   Rcpp::CharacterVector feature_name = get_required_column<CharacterVector>(FT, "feature_name");
//   NumericVector FT_int = get_column_with_fallback<NumericVector>(FT, "intensity_mean", "mean_intensity", "intensity");
//
//
//   Rcpp::NumericVector FT_rt;
//   bool has_rt = false;
//   try {
//     FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention_time");
//     has_rt = FT_rt.size() == FT_mz.size();
//   } catch (...) {
//     has_rt = false;
//   }
//
//   if (!std::is_sorted(FT_mz.begin(), FT_mz.end())) {
//     std::vector<int> idx = compute_sorted_index(FT_mz);
//     FT_mz = reorder_vector(FT_mz, idx);
//     if (has_rt) FT_rt = reorder_vector(FT_rt, idx);
//     feature_name = reorder_vector(feature_name, idx);
//     FT_int = reorder_vector(FT_int, idx);
//   }
//
//   Rcpp::List output_list;
//   int group_index = 1;
//
//   for (int i = 0; i < n_compounds; ++i) {
//     Rcpp::List compound = compounds_list[i];
//     if (!compound.containsElementNamed("Name") || !compound.containsElementNamed("Peaks")) continue;
//
//     std::string compound_name = Rcpp::as<std::string>(compound["Name"]);
//
//     double pre_mz = NA_REAL;
//     try {
//       pre_mz = std::stod(Rcpp::as<std::string>(compound["PreIon"]));
//     } catch (...) {
//       continue;
//     }
//
//     double compound_rt = 0.0;
//     bool has_rt_comp = get_rt(compound, compound_rt);
//
//     double mz_tol = std::max(diff_mz_thr, pre_mz * ppm * 1e-6);
//
//     Rcpp::DataFrame Peaks = Rcpp::as<Rcpp::DataFrame>(compound["Peaks"]);
//     Rcpp::NumericVector peak_mz = Peaks["mz"];
//     Rcpp::NumericVector peak_int = Peaks["rel_intensity"];
//     int total_peaks = peak_mz.size();
//
//     if (!std::is_sorted(peak_mz.begin(), peak_mz.end())) {
//       std::vector<int> idx = compute_sorted_index(peak_mz);
//       peak_mz = reorder_vector(peak_mz, idx);
//       peak_int = reorder_vector(peak_int, idx);
//     }
//
//     for (int j = 0; j < nrow_FT; ++j) {
//       double fmz = FT_mz[j];
//       double frt = has_rt ? FT_rt[j] : NA_REAL;
//
//       if (std::fabs(fmz - pre_mz) > mz_tol) continue;
//       if (has_rt && has_rt_comp && std::fabs(frt - compound_rt) > diff_rt_thr) continue;
//
//       int matched_peaks = 0;
//       std::vector<std::string> matched_peak_features;
//       Rcpp::NumericVector matched_feature_mz, matched_feature_int, match_peak_mz, match_peak_int;
//
//       for (int k = 0; k < nrow_FT; ++k) {
//         double mz_k = FT_mz[k];
//         double int_k = FT_int[k];
//         double rt_k = has_rt ? FT_rt[k] : NA_REAL;
//         if (has_rt && has_rt_comp && std::fabs(rt_k - compound_rt) > diff_rt_thr) continue;
//
//         for (int p = 0; p < total_peaks; ++p) {
//           if (std::fabs(mz_k - peak_mz[p]) <= diff_mz_thr) {
//             matched_peaks++;
//             match_peak_mz.push_back(peak_mz[p]);
//             match_peak_int.push_back(peak_int[p]);
//             matched_peak_features.push_back(Rcpp::as<std::string>(feature_name[k]));
//             matched_feature_mz.push_back(mz_k);
//             matched_feature_int.push_back(int_k);
//             break;
//           }
//         }
//       }
//
//       Rcpp::CharacterVector matched_peak_features_char = Rcpp::wrap(matched_peak_features);
//       double score = matched_peaks > 0 ? matched_peaks / (double)total_peaks : NA_REAL;
//       std::string lib_group_id = "Lib_group_" + std::to_string(group_index++);
//
//       output_list.push_back(List::create(
//           _["feature_name"] = Rcpp::as<std::string>(feature_name[j]),
//           _["Lib_group"] = lib_group_id,
//           _["feature_mz"] = fmz,
//           _["feature_rt"] = frt,
//           _["Compound_name"] = compound_name,
//           _["Compound_rt"] = has_rt_comp ? compound_rt : NA_REAL,
//           _["Compound_peak_mz"] = match_peak_mz.size() ? match_peak_mz : NumericVector::create(NA_REAL),
//           _["Compound_peak_int"] = match_peak_int.size() ? match_peak_int : NumericVector::create(NA_REAL),
//           _["Compound_peak_mz"] = match_peak_mz.size() ? match_peak_mz : NumericVector::create(NA_REAL),
//           _["Compound_peak_mz_full"] = peak_mz,
//           _["Compound_peak_int_full"] = peak_int,
//           _["Precursor_mz"] = pre_mz,
//           _["Matched_peaks"] = matched_peaks,
//           _["Match_peak_score"] = score,
//           _["Total_library_peaks"] = total_peaks,
//           _["Matched_features_to_peaks"] = matched_peak_features_char.size() ? matched_peak_features_char : Rcpp::CharacterVector::create(NA_STRING),
//           _["Matched_features_mz"] = matched_feature_mz.size() ? matched_feature_mz : Rcpp::NumericVector::create(NA_REAL),
//           _["Matched_features_int"] = matched_feature_int.size() ? matched_feature_int : Rcpp::NumericVector::create(NA_REAL)
//
//       ));
//     }
//   }
//
//   return output_list;
// }
//
