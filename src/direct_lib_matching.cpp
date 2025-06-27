// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <limits>
#include "utils_check.h"
#include "utils_sort.h"
using namespace Rcpp;

// Helper: extract RT from compound
inline bool get_rt(const List& compound, double& rt_val) {
  if (!compound.containsElementNamed("RetTime")) return false;
  try {
    rt_val = std::stod(as<std::string>(compound["RetTime"]));
    return true;
  } catch (...) {
    return false;
  }
}

// Helper: extract Precursor Ion from compound
inline bool get_preion(const List& compound, double& preion_mz) {
  if (!compound.containsElementNamed("PreIon")) return false;
  try {
    preion_mz = std::stod(as<std::string>(compound["PreIon"]));
    return true;
  } catch (...) {
    return false;
  }
}

// [[Rcpp::export]]
Rcpp::List direct_match_MS1Lib_cpp(Rcpp::DataFrame FT,
                                   Rcpp::List compounds_list,
                                   double diff_mz_thr,
                                   double diff_rt_thr,
                                   int ppm) {
  int nrow_FT = FT.nrow();
  int nrow_Comp = compounds_list.length();
  NumericVector FT_mz = get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
  CharacterVector feature_name = get_required_column<CharacterVector>(FT, "feature_name");
  NumericVector FT_rt;
  bool has_rt_FT = false;

  try {
    FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention_time");
    has_rt_FT = FT_rt.size() == FT_mz.size();
  } catch (...) {
    has_rt_FT = false;
  }

  if (!std::is_sorted(FT_mz.begin(), FT_mz.end())) {
    std::vector<int> idx = compute_sorted_index(FT_mz);
    FT_mz = reorder_vector(FT_mz, idx);
    if (has_rt_FT) FT_rt = reorder_vector(FT_rt, idx);
    feature_name = reorder_vector(feature_name, idx);
  }

  Rcpp::List output_list;
  int group_index = 1;

  for (int i = 0; i < nrow_Comp; ++i) {
    List compound = compounds_list[i];
    if (!compound.containsElementNamed("Name") ||
        !compound.containsElementNamed("Peaks") ||
        !compound.containsElementNamed("ExactMass")) continue;

        std::string compound_name = compound["Name"];
        double exact_mass = compound["ExactMass"];
        double mz_pi = exact_mass + 1.007276;
        double mz_tol = std::max(diff_mz_thr, mz_pi * ppm * 1e-6);

        double compound_rt = 0.0;
        bool has_rt_comp = get_rt(compound, compound_rt);

        DataFrame Peaks = as<DataFrame>(compound["Peaks"]);

        NumericVector peak_mz = Peaks["mz"];
        NumericVector peak_int = Peaks["rel_intensity"];
        int total_peaks = peak_mz.size();

        if (!std::is_sorted(peak_mz.begin(), peak_mz.end())) {
          std::vector<int> idx = compute_sorted_index(peak_mz);
          peak_mz = reorder_vector(peak_mz, idx);
          peak_int = reorder_vector(peak_int, idx);
        }

        for (int j = 0; j < nrow_FT; ++j) {
          double fmz = FT_mz[j];
          double frt = has_rt_FT ? FT_rt[j] : NA_REAL;
          if (std::fabs(fmz - mz_pi) > mz_tol) continue;
          if (has_rt_FT && has_rt_comp && std::fabs(frt - compound_rt) > diff_rt_thr) continue;

          int matched_peaks = 0;
          std::vector<std::string> matched_peak_features;
          Rcpp::NumericVector matched_feature_mz, match_peak_mz, match_peak_int;

          for (int k = 0; k < nrow_FT; ++k) {
            double mz_k = FT_mz[k];
            double rt_k = has_rt_FT ? FT_rt[k] : NA_REAL;
            if (has_rt_FT && has_rt_comp && std::fabs(rt_k - compound_rt) > diff_rt_thr) continue;

            for (int p = 0; p < total_peaks; ++p) {
              if (std::fabs(mz_k - peak_mz[p]) <= diff_mz_thr) {
                matched_peaks++;
                match_peak_mz.push_back(peak_mz[p]);
                match_peak_int.push_back(peak_int[p]);
                matched_peak_features.push_back(as<std::string>(feature_name[k]));
                matched_feature_mz.push_back(mz_k);
                break;
              }
            }
          }
          Rcpp::CharacterVector matched_peak_features_char = Rcpp::wrap(matched_peak_features);
          double score = matched_peaks > 0 ? matched_peaks / (double)total_peaks : NA_REAL;
          std::string lib_group_id = "Lib_group_" + std::to_string(group_index++);

          output_list.push_back(List::create(
              _["feature_name"] = as<std::string>(feature_name[j]),
              _["Lib_group"] = lib_group_id,
              _["feature_mz"] = fmz,
              _["feature_rt"] = frt,
              _["Compound_name"] = compound_name,
              _["Compound_rt"] = has_rt_comp ? compound_rt : NA_REAL,
              _["Compound_peak_mz"] = match_peak_mz.size() ? match_peak_mz : NumericVector::create(NA_REAL),
              _["Compound_peak_int"] = match_peak_int.size() ? match_peak_int : NumericVector::create(NA_REAL),
              _["Ion_type"] = "[M+H]+",
              _["Matched_peaks"] = matched_peaks,
              _["Match_peak_score"] = score,
              _["Total_library_peaks"] = total_peaks,
              _["Matched_features_to_peaks"] = matched_peak_features_char.size() > 0 ? matched_peak_features_char : Rcpp::CharacterVector::create(NA_STRING),
              _["Matched_features_mz"] = matched_feature_mz.size() > 0 ? matched_feature_mz : Rcpp::NumericVector::create(NA_REAL)
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
  }

  Rcpp::List output_list;
  int group_index = 1;

  for (int i = 0; i < n_compounds; ++i) {
    Rcpp::List compound = compounds_list[i];
    if (!compound.containsElementNamed("Name") || !compound.containsElementNamed("Peaks")) continue;

    std::string compound_name = Rcpp::as<std::string>(compound["Name"]);

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

      int matched_peaks = 0;
      std::vector<std::string> matched_peak_features;
      Rcpp::NumericVector matched_feature_mz, match_peak_mz, match_peak_int;

      for (int k = 0; k < nrow_FT; ++k) {
        double mz_k = FT_mz[k];
        double rt_k = has_rt ? FT_rt[k] : NA_REAL;
        if (has_rt && has_rt_comp && std::fabs(rt_k - compound_rt) > diff_rt_thr) continue;

        for (int p = 0; p < total_peaks; ++p) {
          if (std::fabs(mz_k - peak_mz[p]) <= diff_mz_thr) {
            matched_peaks++;
            match_peak_mz.push_back(peak_mz[p]);
            match_peak_int.push_back(peak_int[p]);
            matched_peak_features.push_back(Rcpp::as<std::string>(feature_name[k]));
            matched_feature_mz.push_back(mz_k);
            break;
          }
        }
      }

      Rcpp::CharacterVector matched_peak_features_char = Rcpp::wrap(matched_peak_features);
      double score = matched_peaks > 0 ? matched_peaks / (double)total_peaks : NA_REAL;
      std::string lib_group_id = "Lib_group_" + std::to_string(group_index++);

      output_list.push_back(List::create(
          _["feature_name"] = Rcpp::as<std::string>(feature_name[j]),
          _["Lib_group"] = lib_group_id,
          _["feature_mz"] = fmz,
          _["feature_rt"] = frt,
          _["Compound_name"] = compound_name,
          _["Compound_rt"] = has_rt_comp ? compound_rt : NA_REAL,
          _["Compound_peak_mz"] = match_peak_mz.size() ? match_peak_mz : NumericVector::create(NA_REAL),
          _["Compound_peak_int"] = match_peak_int.size() ? match_peak_int : NumericVector::create(NA_REAL),
          _["Precursor_mz"] = pre_mz,
          _["Matched_peaks"] = matched_peaks,
          _["Match_peak_score"] = score,
          _["Total_library_peaks"] = total_peaks,
          _["Matched_features_to_peaks"] = matched_peak_features_char.size() ? matched_peak_features_char : Rcpp::CharacterVector::create(NA_STRING),
          _["Matched_features_mz"] = matched_feature_mz.size() ? matched_feature_mz : Rcpp::NumericVector::create(NA_REAL)
      ));
    }
  }

  return output_list;
}

