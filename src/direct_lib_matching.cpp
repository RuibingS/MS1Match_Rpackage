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
  if (!compound.containsElementNamed("RetentionTime")) return false;
  SEXP rt = compound["RetentionTime"];
  if (Rf_isString(rt)) {
    try {
      rt_val = std::stod(as<std::string>(rt));
      return !std::isnan(rt_val);
    } catch (...) {
      return false;
    }
  }
  if (Rf_isReal(rt) || Rf_isInteger(rt)) {
    rt_val = as<double>(rt);
    return !ISNA(rt_val) && !std::isnan(rt_val);
  }
  return false;
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

  if (ion_mode == "positive" || ion_mode == "pos") {
    if (adduct == "m+h" || adduct == "[m+h]+")
      return exact_mass + proton;
    if (adduct == "m+na" || adduct == "[m+na]+")
      return exact_mass + sodium;
    if (adduct == "m+k" || adduct == "[m+k]+")
      return exact_mass + potassium;
    return exact_mass + proton;
  }

  if (ion_mode == "negative" || ion_mode == "neg") {
    if (adduct == "m-h" || adduct == "[m-h]-")
      return exact_mass - proton;
    return exact_mass - proton;
  }

  return exact_mass + proton;
}


// --------------------------------------------------------
// MS1 matching
//
// Two-Step RT logic:
//   rt_thr_compound        — loose threshold, precursor gate vs compound library RT
//   rt_thr_coelution — tight threshold, co-elution pool anchored to feature j
// --------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List direct_match_MS1Lib_cpp(Rcpp::DataFrame FT,
                                   Rcpp::List compounds_list,
                                   double rt_thr_compound,
                                   double rt_thr_coelution,
                                   int ppm) {

  int nrow_FT = FT.nrow();

  NumericVector FT_mz =
    get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
  NumericVector FT_int =
    get_column_with_fallback<NumericVector>(FT, "intensity_mean",
                                            "mean_intensity", "intensity");
  CharacterVector feature_name =
    get_required_column<CharacterVector>(FT, "feature_name");

  NumericVector FT_rt;
  bool has_rt_FT = false;
  try {
    FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention_time");
    has_rt_FT = (FT_rt.size() == FT_mz.size());
  } catch (...) {
    has_rt_FT = false;
  }

  if (!std::is_sorted(FT_mz.begin(), FT_mz.end())) {
    std::vector<int> idx = compute_sorted_index(FT_mz);
    FT_mz        = reorder_vector(FT_mz, idx);
    FT_int       = reorder_vector(FT_int, idx);
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

    double mz_pi  = compute_ion_mz(exact_mass, ion_mode, adduct);
    double mz_tol = mz_pi * ppm * 1e-6;

    double compound_rt  = 0.0;
    bool   has_rt_comp  = get_rt(compound, compound_rt);

    DataFrame   Peaks    = as<DataFrame>(compound["Peaks"]);
    NumericVector peak_mz  = Peaks["mz"];
    NumericVector peak_int = Peaks["rel_intensity"];
    int total_peaks = peak_mz.size();

    if (!std::is_sorted(peak_mz.begin(), peak_mz.end())) {
      std::vector<int> idx = compute_sorted_index(peak_mz);
      peak_mz  = reorder_vector(peak_mz, idx);
      peak_int = reorder_vector(peak_int, idx);
    }

    for (int j = 0; j < nrow_FT; ++j) {

      double fmz = FT_mz[j];
      double frt = has_rt_FT ? FT_rt[j] : NA_REAL;

      // ---- Step 1: precursor m/z gate ----
      if (std::fabs(fmz - mz_pi) > mz_tol) continue;

      // ---- loose RT gate vs compound library RT ----
      if (has_rt_comp && has_rt_FT &&
          std::fabs(frt - compound_rt) > rt_thr_compound) continue;

      // ---- Step 2: co-elution pool anchored strictly to feature j ----
      // Only features that co-elute with this precursor feature
      std::vector<int> candidates;
      for (int k = 0; k < nrow_FT; ++k) {
        if (has_rt_FT && !ISNA(frt)) {
          if (std::fabs(FT_rt[k] - frt) > rt_thr_coelution) continue;
        }
        candidates.push_back(k);
      }

      NumericVector   sub_mz(candidates.size());
      NumericVector   sub_rt(candidates.size());
      NumericVector   sub_int(candidates.size());
      CharacterVector sub_name(candidates.size());

      for (size_t idx = 0; idx < candidates.size(); ++idx) {
        int k        = candidates[idx];
        sub_mz[idx]  = FT_mz[k];
        sub_rt[idx]  = has_rt_FT ? FT_rt[k] : NA_REAL;
        sub_int[idx] = FT_int[k];
        sub_name[idx]= feature_name[k];
      }

      // Peak matching: pure ppm tolerance scaled to each library peak mz.
      std::vector<std::string> matched_peak_features;
      NumericVector matched_feature_mz, matched_feature_rt, matched_feature_int;
      NumericVector match_peak_mz, match_peak_int;

      std::vector<bool> sub_used(candidates.size(), false);
      int matched_peaks = 0;

      for (int pk = 0; pk < total_peaks; ++pk) {
        double lib_mz   = peak_mz[pk];
        double peak_tol = lib_mz * ppm * 1e-6;
        double best_diff = peak_tol + 1.0;
        int    best_s    = -1;

        for (size_t s = 0; s < candidates.size(); ++s) {
          if (sub_used[s]) continue;
          double d = std::fabs(sub_mz[s] - lib_mz);
          if (d <= peak_tol && d < best_diff) {
            best_diff = d;
            best_s    = (int)s;
          }
        }
        if (best_s >= 0) {
          sub_used[best_s] = true;
          matched_peaks++;
          matched_peak_features.push_back(as<std::string>(sub_name[best_s]));
          matched_feature_mz.push_back(sub_mz[best_s]);
          matched_feature_rt.push_back(sub_rt[best_s]);
          matched_feature_int.push_back(sub_int[best_s]);
          match_peak_mz.push_back(lib_mz);
          match_peak_int.push_back(peak_int[pk]);
        }
      }

      double score = matched_peaks > 0 ?
      static_cast<double>(matched_peaks) / total_peaks : NA_REAL;

      std::string lib_group_id = "Lib_group_" + std::to_string(group_index++);

      CharacterVector matched_peak_features_char;
      if (matched_peak_features.empty())
        matched_peak_features_char = CharacterVector::create(NA_STRING);
      else
        matched_peak_features_char = wrap(matched_peak_features);

      output_list.push_back(List::create(
          _["feature_name"]   = as<std::string>(feature_name[j]),
          _["Lib_group"]      = lib_group_id,
          _["feature_mz"]     = fmz,
          _["feature_rt"]     = frt,
          _["Compound_name"]          = compound_name,
          _["Compound_M"]             = exact_mass,
          _["Compound_rt"]            = has_rt_comp ? compound_rt : NA_REAL,
          _["Compound_peak_mz"]       = match_peak_mz.size()  ? match_peak_mz  : NumericVector::create(NA_REAL),
          _["Compound_peak_int"]      = match_peak_int.size() ? match_peak_int : NumericVector::create(NA_REAL),
          _["Compound_peak_mz_full"]  = peak_mz,
          _["Compound_peak_int_full"] = peak_int,
          _["Total_library_peaks"]    = total_peaks,
          _["Ion_type"]               = adduct,
          _["Matched_peaks"]          = matched_peaks,
          _["Match_peak_score"]       = score,
          _["Matched_features_to_peaks"] = matched_peak_features_char,
          _["Matched_features_mz"]    = matched_feature_mz.size()  ? matched_feature_mz  : NumericVector::create(NA_REAL),
          _["Matched_features_rt"]    = matched_feature_rt.size()  ? matched_feature_rt  : NumericVector::create(NA_REAL),
          _["Matched_features_int"]   = matched_feature_int.size() ? matched_feature_int : NumericVector::create(NA_REAL)
      ));
    }
  }

  return output_list;
}




// --------------------------------------------------------
// Helper: extract PrecursorMZ from compound (numeric OR character safe)
// --------------------------------------------------------
inline bool get_precursor_mz(const List& compound, double& mz_val) {
  if (!compound.containsElementNamed("PrecursorMZ")) return false;
  SEXP pm = compound["PrecursorMZ"];
  if (Rf_isString(pm)) {
    try {
      mz_val = std::stod(as<std::string>(pm));
      return !std::isnan(mz_val);
    } catch (...) {
      return false;
    }
  }
  if (Rf_isReal(pm) || Rf_isInteger(pm)) {
    mz_val = as<double>(pm);
    return !ISNA(mz_val) && !std::isnan(mz_val);
  }
  return false;
}


// --------------------------------------------------------
// MS2 matching
//
// Two-Step RT logic (same as MS1):
//   rt_thr_compound        — loose threshold, precursor gate vs compound library RT
//   rt_thr_coelution — tight threshold, fragment pool anchored to feature j
//
// --------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List direct_match_MS2Lib_cpp(Rcpp::DataFrame FT,
                                   Rcpp::List compounds_list,
                                   double rt_thr_compound,
                                   double rt_thr_coelution,
                                   int ppm,
                                   bool verbose = true) {

  int nrow_FT    = FT.nrow();
  int n_compounds = compounds_list.length();

  Rcpp::NumericVector FT_mz =
    get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
  Rcpp::CharacterVector feature_name =
    get_required_column<CharacterVector>(FT, "feature_name");
  Rcpp::NumericVector FT_int =
    get_column_with_fallback<NumericVector>(FT, "intensity_mean",
                                            "mean_intensity", "intensity");

  Rcpp::NumericVector FT_rt;
  bool has_rt = false;
  try {
    FT_rt  = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention_time");
    has_rt = FT_rt.size() == FT_mz.size();
  } catch (...) {
    has_rt = false;
  }

  if (!std::is_sorted(FT_mz.begin(), FT_mz.end())) {
    std::vector<int> idx = compute_sorted_index(FT_mz);
    FT_mz        = reorder_vector(FT_mz, idx);
    if (has_rt) FT_rt = reorder_vector(FT_rt, idx);
    feature_name = reorder_vector(feature_name, idx);
    FT_int       = reorder_vector(FT_int, idx);
  }

  Rcpp::List output_list;
  int group_index = 1;
  int n_skipped = 0;
  int n_matched_precursors = 0;

  for (int i = 0; i < n_compounds; ++i) {

    Rcpp::List compound = compounds_list[i];

    std::string compound_name = compound.containsElementNamed("Name") ?
    Rcpp::as<std::string>(compound["Name"]) : ("compound_#" + std::to_string(i + 1));

    if (!compound.containsElementNamed("Name") ||
        !compound.containsElementNamed("Peaks")) {
        if (verbose) Rcpp::Rcout << "Skipping '" << compound_name
                                 << "': missing Name or Peaks\n";
        n_skipped++;
        continue;
    }

    if (!compound.containsElementNamed("ExactMass")) {
      if (verbose) Rcpp::Rcout << "Skipping '" << compound_name
                               << "': missing ExactMass\n";
      n_skipped++;
      continue;
    }
    double exact_mass = as<double>(compound["ExactMass"]);

    double pre_mz = NA_REAL;
    if (!get_precursor_mz(compound, pre_mz)) {
      if (verbose) Rcpp::Rcout << "Skipping '" << compound_name
                               << "': missing/unparseable PrecursorMZ\n";
      n_skipped++;
      continue;
    }

    double compound_rt = 0.0;
    bool   has_rt_comp = get_rt(compound, compound_rt);
    double mz_tol      = pre_mz * ppm * 1e-6;

    Rcpp::DataFrame   Peaks    = Rcpp::as<Rcpp::DataFrame>(compound["Peaks"]);
    Rcpp::NumericVector peak_mz  = Peaks["mz"];
    Rcpp::NumericVector peak_int = Peaks["rel_intensity"];
    int total_peaks = peak_mz.size();

    if (total_peaks == 0) {
      if (verbose) Rcpp::Rcout << "Skipping '" << compound_name
                               << "': Peaks table is empty\n";
      n_skipped++;
      continue;
    }

    if (!std::is_sorted(peak_mz.begin(), peak_mz.end())) {
      std::vector<int> idx = compute_sorted_index(peak_mz);
      peak_mz  = reorder_vector(peak_mz,  idx);
      peak_int = reorder_vector(peak_int, idx);
    }

    bool any_precursor_hit = false;

    for (int j = 0; j < nrow_FT; ++j) {

      double fmz = FT_mz[j];
      double frt = has_rt ? FT_rt[j] : NA_REAL;

      // ---- Step 1: precursor m/z gate ----
      if (std::fabs(fmz - pre_mz) > mz_tol) continue;

      // ---- Step 1: loose RT gate vs compound library RT ----
      if (has_rt && has_rt_comp &&
          std::fabs(frt - compound_rt) > rt_thr_compound) continue;

      any_precursor_hit = true;


      std::vector<std::string> matched_peak_features;
      Rcpp::NumericVector matched_feature_mz, matched_feature_int,
      matched_feature_rt;
      Rcpp::NumericVector match_peak_mz, match_peak_int;

      matched_peak_features.push_back(Rcpp::as<std::string>(feature_name[j]));
      matched_feature_mz.push_back(fmz);
      matched_feature_rt.push_back(frt);
      matched_feature_int.push_back(FT_int[j]);
      match_peak_mz.push_back(pre_mz);
      match_peak_int.push_back(100.0);

      // ---- Step 2: fragment matching ----
    
      int matched_peaks = 0;  // fragment count only (score denominator)

      std::vector<bool> peak_used(total_peaks, false);
      std::vector<bool> ft_used(nrow_FT, false);
      ft_used[j] = true;  // precursor already consumed above

      for (int p = 0; p < total_peaks; ++p) {
        double lib_mz    = peak_mz[p];
        // Sentinel: larger than any real mz difference.
        double best_diff = 1e9;
        int    best_k    = -1;

        for (int k = 0; k < nrow_FT; ++k) {
          if (ft_used[k]) continue;

          // ---- Step 2: strict co-elution window anchored to feature j ----
          if (has_rt && !ISNA(frt)) {
            if (std::fabs(FT_rt[k] - frt) > rt_thr_coelution) continue;
          }

          double frag_tol = lib_mz * ppm * 1e-6;
          double mz_diff = std::fabs(lib_mz - FT_mz[k]);
          if (mz_diff <= frag_tol && mz_diff < best_diff) {
            best_diff = mz_diff;
            best_k    = k;
          }
        }

        if (best_k >= 0) {
          peak_used[p]    = true;
          ft_used[best_k] = true;
          matched_peaks++;

          match_peak_mz.push_back(lib_mz);
          match_peak_int.push_back(peak_int[p]);
          matched_peak_features.push_back(
            Rcpp::as<std::string>(feature_name[best_k]));
          matched_feature_mz.push_back(FT_mz[best_k]);
          matched_feature_rt.push_back(has_rt ? FT_rt[best_k] : NA_REAL);
          matched_feature_int.push_back(FT_int[best_k]);
        }
      }

      double score = matched_peaks > 0 ?
      (double)matched_peaks / total_peaks : NA_REAL;

      std::string lib_group_id = "Lib_group_" + std::to_string(group_index++);

      // keep pre_mz is present in the full peak list for build_mirror_df
      bool prec_in_full = false;
      for (int p = 0; p < total_peaks; ++p) {
        if (std::fabs(peak_mz[p] - pre_mz) < 0.01) { prec_in_full = true; break; }
      }
      Rcpp::NumericVector peak_mz_full  = Rcpp::clone(peak_mz);
      Rcpp::NumericVector peak_int_full = Rcpp::clone(peak_int);
      if (!prec_in_full) {
        peak_mz_full.push_back(pre_mz);
        peak_int_full.push_back(100.0);
      }

      output_list.push_back(List::create(
          _["feature_name"]  = Rcpp::as<std::string>(feature_name[j]),
          _["Lib_group"]     = lib_group_id,
          _["feature_mz"]    = fmz,
          _["feature_rt"]    = frt,
          _["Compound_name"] = compound_name,
          _["Compound_M"]    = exact_mass,
          _["Compound_rt"]   = has_rt_comp ? compound_rt : NA_REAL,
          _["Precursor_mz"]  = pre_mz,
          _["Compound_peak_mz"]       = match_peak_mz,
          _["Compound_peak_int"]      = match_peak_int,
          _["Compound_peak_mz_full"]  = peak_mz_full,
          _["Compound_peak_int_full"] = peak_int_full,
          _["Matched_peaks"]          = matched_peaks,
          _["Match_peak_score"]       = score,
          _["Total_library_peaks"]    = total_peaks,
          _["Ion_type"]               = "Precursor",
          _["Matched_features_to_peaks"] = Rcpp::wrap(matched_peak_features),
          _["Matched_features_mz"]    = matched_feature_mz,
          _["Matched_features_rt"]    = matched_feature_rt,
          _["Matched_features_int"]   = matched_feature_int
      ));
    }

    if (any_precursor_hit) {
      n_matched_precursors++;
    } else if (verbose) {
      Rcpp::Rcout << "'" << compound_name
                  << "': PrecursorMZ/RT parsed OK, but no feature within "
                  << "ppm/RT window (pre_mz=" << pre_mz << ")\n";
    }
  }

  if (verbose) {
    Rcpp::Rcout << "\n[direct_match_MS2Lib_cpp] "
                << n_compounds << " compounds in library, "
                << n_skipped << " skipped (bad/missing metadata), "
                << n_matched_precursors << " had >=1 precursor-level feature match\n";
  }

  return output_list;
}
