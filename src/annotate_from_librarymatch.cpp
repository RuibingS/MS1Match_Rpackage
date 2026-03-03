// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <limits>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>  // <-- add this if you use std::iota

#include "utils_annotation.h"
#include "utils_sort.h"
#include "utils_check.h"
#include "utils_isotope.h"


using namespace Rcpp;




// -----------------------------------------------------------------------------
// helpers
// -----------------------------------------------------------------------------
static inline std::string tolower_copy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  return s;
}

static inline double mz_tol_mix(double mz, double mz_threshold, double ppm) {
  return std::max(mz_threshold, mz * ppm * 1e-6);
}

// -----------------------------------------------------------------------------
// Main: annotate from precomputed library matches (MS1 or MS2)
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List annotate_from_librarymatch_cpp(Rcpp::List library_matches,
                                          Rcpp::DataFrame adducts_df,
                                          Rcpp::DataFrame neutral_loss_table,
                                          std::string polarity,
                                          double mz_threshold = 0.005,
                                          double ppm = 20,
                                          double rt_threshold = 3.0,
                                          double iso_threshold = 1.003355,
                                          int max_charge = 3,
                                          int max_isotope = 4) {

  // ---------------- helpers ----------------
  auto tolower_copy = [](std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
  };
  auto mz_tol_mix = [](double mz, double mz_thr, double ppm_) {
    return std::max(mz_thr, mz * ppm_ * 1e-6);
  };
  auto collapse_matches = [](const std::vector<std::string>& matches) {
    if (matches.empty()) return std::string("");
    std::ostringstream oss;
    oss << matches[0];
    for (size_t i = 1; i < matches.size(); ++i) oss << "; " << matches[i];
    return oss.str();
  };


  polarity = tolower_copy(polarity);
  if (polarity == "positive" || polarity == "+" || polarity == "pos") polarity = "pos";
  else if (polarity == "negative" || polarity == "-" || polarity == "neg") polarity = "neg";

  // ---- adduct table
  Rcpp::CharacterVector adduct_ions   = adducts_df["Ion.name"];
  Rcpp::CharacterVector adduct_pol    = adducts_df["Ion.polarity"];
  Rcpp::CharacterVector adduct_mass   = adducts_df["Ion.mass"];
  Rcpp::CharacterVector adduct_charge = adducts_df.containsElementNamed("Charge")
    ? (Rcpp::CharacterVector)adducts_df["Charge"]
  : Rcpp::CharacterVector(adduct_ions.size(), NA_STRING);

  struct KeyHash {
    std::size_t operator()(const std::string& s) const noexcept {
      return std::hash<std::string>{}(s);
    }
  };

  std::unordered_set<std::string, KeyHash> adduct_seen;
  std::vector<std::string> filtered_ions;
  std::vector<std::string> filtered_formulas;

  for (int i = 0; i < adduct_ions.size(); ++i) {
    std::string pol = tolower_copy(Rcpp::as<std::string>(adduct_pol[i]));
    if (pol == "positive" || pol == "+" || pol == "pos") pol = "pos";
    else if (pol == "negative" || pol == "-" || pol == "neg") pol = "neg";
    if (pol != polarity) continue;

    std::string ion  = Rcpp::as<std::string>(adduct_ions[i]);
    std::string mass = Rcpp::as<std::string>(adduct_mass[i]);
    std::string chg  = (adduct_charge[i] == NA_STRING) ? "" : Rcpp::as<std::string>(adduct_charge[i]);

    std::string key = ion + "||" + mass + "||" + chg + "||" + pol;
    if (adduct_seen.insert(key).second) {
      filtered_ions.push_back(ion);
      filtered_formulas.push_back(mass);
    }
  }

  // ---- neutral loss table
  Rcpp::NumericVector nl_mass = neutral_loss_table["Accurate.Mass"];
  Rcpp::CharacterVector nl_name = neutral_loss_table.containsElementNamed("Formula")
    ? (Rcpp::CharacterVector)neutral_loss_table["Formula"]
  : (Rcpp::CharacterVector)neutral_loss_table["Neutral.Loss"];
  Rcpp::CharacterVector nl_pos = neutral_loss_table["Pos"];
  Rcpp::CharacterVector nl_neg = neutral_loss_table["Neg"];

  // ---- output (list-of-rows, like your other function)
  std::vector<Rcpp::List> rows;
  rows.reserve(4096);

  const int n_groups = library_matches.size();

  for (int g = 0; g < n_groups; ++g) {
    Rcpp::List grp = library_matches[g];

    // --- group metadata
    std::string compound_name = grp.containsElementNamed("Compound_name")
      ? Rcpp::as<std::string>(grp["Compound_name"]) : "";

    std::string lib_group = grp.containsElementNamed("Lib_group")
      ? Rcpp::as<std::string>(grp["Lib_group"]) : "";

    std::string ion_type = grp.containsElementNamed("Ion_type")
      ? Rcpp::as<std::string>(grp["Ion_type"]) : "";

    // --- neutral mass M:
    // prefer Compound_M (MS1 match), otherwise derive from Precursor_mz (MS2 match) assuming precursor is [M+H]+ (pos) or [M-H]- (neg)
    double M = grp.containsElementNamed("Compound_M") ? Rcpp::as<double>(grp["Compound_M"]) : NA_REAL;

    double precursor_mz = NA_REAL;
    if (Rcpp::NumericVector::is_na(M) && grp.containsElementNamed("Precursor_mz")) {
      precursor_mz = Rcpp::as<double>(grp["Precursor_mz"]);
      if (!Rcpp::NumericVector::is_na(precursor_mz)) {
        // assume precursor is primary ion
        M = (polarity == "pos") ? (precursor_mz - 1.007276) : (precursor_mz + 1.007276);
      }
    }

    // vectors per group
    Rcpp::CharacterVector fn = grp.containsElementNamed("Matched_features_to_peaks")
      ? (Rcpp::CharacterVector)grp["Matched_features_to_peaks"] : Rcpp::CharacterVector();

    Rcpp::NumericVector mz = grp.containsElementNamed("Matched_features_mz")
      ? (Rcpp::NumericVector)grp["Matched_features_mz"] : Rcpp::NumericVector();

    Rcpp::NumericVector rt = grp.containsElementNamed("Matched_features_rt")
      ? (Rcpp::NumericVector)grp["Matched_features_rt"] : Rcpp::NumericVector(mz.size(), NA_REAL);

    Rcpp::NumericVector inten = grp.containsElementNamed("Matched_features_int")
      ? (Rcpp::NumericVector)grp["Matched_features_int"] : Rcpp::NumericVector(mz.size(), NA_REAL);

    if (fn.size() == 0 || mz.size() == 0) continue;
    if (mz.size() != fn.size()) Rcpp::stop("Matched_features_mz length != Matched_features_to_peaks length.");
    if (rt.size() != mz.size()) rt = Rcpp::NumericVector(mz.size(), NA_REAL);
    if (inten.size() != mz.size()) inten = Rcpp::NumericVector(mz.size(), NA_REAL);

    // optional group scores
    double inner_cos = grp.containsElementNamed("inner_cosine") ? Rcpp::as<double>(grp["inner_cosine"]) : NA_REAL;
    double outer_cos = grp.containsElementNamed("outer_cosine") ? Rcpp::as<double>(grp["outer_cosine"]) : NA_REAL;
    double match_sc  = grp.containsElementNamed("Match_peak_score") ? Rcpp::as<double>(grp["Match_peak_score"]) : NA_REAL;

    // --- precompute theoretical adduct m/z for this M (once per group)
    std::vector<double> theo_adduct_mz;
    if (!Rcpp::NumericVector::is_na(M) && !filtered_formulas.empty()) {
      theo_adduct_mz = parse_adduct_vector(M, filtered_formulas); // from utils_annotation.h
    }

    // --- PI tagging:
    // if Ion_type says M+H / M-H -> compute PI target from M
    // else if Ion_type says Precursor and we have precursor_mz -> PI target is precursor_mz
    double PI_target = NA_REAL;
    std::string PI_tag = "";

    {
      std::string it = tolower_copy(ion_type);

      if (!Rcpp::NumericVector::is_na(precursor_mz) && it.find("precursor") != std::string::npos) {
        PI_target = precursor_mz;
        PI_tag = (polarity == "pos") ? " [M+H]+ " : " [M-H]- ";
      } else if (!Rcpp::NumericVector::is_na(M)) {
        if (it.find("m+h") != std::string::npos) {
          PI_target = M + 1.007276;
          PI_tag = " [M+H]+ ";
        } else if (it.find("m-h") != std::string::npos) {
          PI_target = M - 1.007276;
          PI_tag = " [M-H]- ";
        }
      }
    }

    int pi_idx = -1;
    if (!Rcpp::NumericVector::is_na(PI_target)) {
      double best = 1e99;
      for (int i = 0; i < mz.size(); ++i) {
        double d = std::fabs(mz[i] - PI_target);
        if (d < best) { best = d; pi_idx = i; }
      }
      if (best > mz_threshold) pi_idx = -1;
    }

    // --- isotope annotation within group (sort by m/z, run Check_Iso_Charge, map back)
    std::vector<int> ord(mz.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return mz[a] < mz[b]; });

    Rcpp::CharacterVector fn_s(fn.size());
    Rcpp::NumericVector mz_s(mz.size()), rt_s(rt.size()), in_s(inten.size());
    for (int k = 0; k < (int)ord.size(); ++k) {
      fn_s[k] = fn[ord[k]];
      mz_s[k] = mz[ord[k]];
      rt_s[k] = rt[ord[k]];
      in_s[k] = inten[ord[k]];
    }

    Rcpp::CharacterVector iso_s = Check_Iso_Charge(fn_s, mz_s, rt_s, in_s,
                                                   iso_threshold, max_charge, max_isotope);

    std::vector<std::string> iso_label(mz.size(), "");
    for (int k = 0; k < iso_s.size(); ++k) {
      int orig = ord[k];
      iso_label[orig] = Rcpp::as<std::string>(iso_s[k]);
    }

    // --- NL matching within group
    std::vector<std::string> nl_label(mz.size(), "");
    for (int i = 0; i < mz.size(); ++i) {
      std::vector<std::string> hits;

      for (int j = 0; j < mz.size(); ++j) {
        if (i == j) continue;
        if (mz[j] <= mz[i]) continue;

        // RT constraint (only if both available)
        if (!Rcpp::NumericVector::is_na(rt[i]) && !Rcpp::NumericVector::is_na(rt[j])) {
          if (std::fabs(rt[i] - rt[j]) > rt_threshold) continue;
        }

        double diff = mz[j] - mz[i];

        for (int k = 0; k < nl_mass.size(); ++k) {
          bool pol_ok = (polarity == "pos")
          ? (Rcpp::as<std::string>(nl_pos[k]) == "+")
            : (Rcpp::as<std::string>(nl_neg[k]) == "+");
          if (!pol_ok) continue;

          if (std::fabs(diff - nl_mass[k]) <= mz_threshold) {
            std::string nm = Rcpp::as<std::string>(nl_name[k]);
            std::string fj = Rcpp::as<std::string>(fn[j]);
            hits.push_back("<- " + nm + " (" + fj + ")");
          }
        }
      }
      nl_label[i] = join_annotations(hits); // from utils_annotation.h
    }

    // --- adduct match within group (single M)
    std::vector<std::string> adduct_label(mz.size(), "");
    if (!Rcpp::NumericVector::is_na(M) && !theo_adduct_mz.empty()) {
      for (int i = 0; i < mz.size(); ++i) {
        std::unordered_set<std::string, KeyHash> seen;
        std::vector<std::string> hits;

        size_t L = std::min(theo_adduct_mz.size(), filtered_ions.size());
        for (size_t k = 0; k < L; ++k) {
          double theo = theo_adduct_mz[k];
          if (std::isnan(theo)) continue;

          double d   = std::fabs(mz[i] - theo);
          double tol = mz_tol_mix(theo, mz_threshold, ppm);

          if (d <= tol) {
            std::ostringstream oss;
            oss << filtered_ions[k] << " " << std::fixed << std::setprecision(6) << M;
            std::string ann = oss.str();
            if (seen.insert(ann).second) hits.push_back(ann);
          }
        }
        adduct_label[i] = join_annotations(hits);
      }
    }

    // --- push rows (list-of-rows output)
    for (int i = 0; i < mz.size(); ++i) {
      std::string PI_str = "";
      if (i == pi_idx && !Rcpp::NumericVector::is_na(M)) {
        std::ostringstream ms;
        ms << std::fixed << std::setprecision(5) << M;
        PI_str = compound_name + PI_tag + ms.str();
      }

      rows.push_back(Rcpp::List::create(
          Rcpp::Named("feature_name")     = Rcpp::CharacterVector::create(Rcpp::as<std::string>(fn[i])),
          Rcpp::Named("mz")               = mz[i],
          Rcpp::Named("rt")               = rt[i],
          Rcpp::Named("intensity")        = inten[i],
          Rcpp::Named("Compound_name")    = compound_name,
         Rcpp::Named("Lib_group")        = lib_group,
        Rcpp::Named("Ion_type")         = ion_type,
         Rcpp::Named("Compound_M")       = M,
         Rcpp::Named("PI_match")         = PI_str,
         Rcpp::Named("adduct_match")     = adduct_label[i],
         Rcpp::Named("NL_match")         = nl_label[i],
                                                                                                                                                                                                                 Rcpp::Named("iso_annotation")   = Rcpp::CharacterVector::create(iso_label[i]),
                                                                                                                                                                                                                 Rcpp::Named("Match_peak_score") = match_sc,
                                                                                                                                                                                                                 Rcpp::Named("inner_cosine")     = inner_cos,
                                                                                                                                                                                                                 Rcpp::Named("outer_cosine")     = outer_cos
      ));
    }
  }

  // return as List
  Rcpp::List output_list(rows.size());
  for (int i = 0; i < (int)rows.size(); ++i) output_list[i] = rows[i];
  return output_list;
}




// Rcpp::DataFrame annotate_from_librarymatch_cpp(Rcpp::List library_matches,
//                                                Rcpp::DataFrame adducts_df,
//                                                Rcpp::DataFrame neutral_loss_table,
//                                                std::string polarity,
//                                                double mz_threshold = 0.005,
//                                                double ppm = 20,
//                                                double rt_threshold = 3.0,
//                                                double iso_threshold = 1.003355,
//                                                int max_charge = 3,
//                                                int max_isotope = 4) {
//   // ---- normalize polarity
//   polarity = tolower_copy(polarity);
//   if (polarity == "positive" || polarity == "+" || polarity == "pos") polarity = "pos";
//   else if (polarity == "negative" || polarity == "-" || polarity == "neg") polarity = "neg";
//
//   // ---- adduct table filtering (by polarity)
//   CharacterVector adduct_ions   = adducts_df["Ion.name"];
//   CharacterVector adduct_pol    = adducts_df["Ion.polarity"];
//   CharacterVector adduct_mass   = adducts_df["Ion.mass"];
//   CharacterVector adduct_charge = adducts_df.containsElementNamed("Charge")
//     ? (CharacterVector)adducts_df["Charge"]
//   : CharacterVector(adduct_ions.size(), NA_STRING);
//
//   struct KeyHash {
//     std::size_t operator()(const std::string& s) const noexcept {
//       return std::hash<std::string>{}(s);
//     }
//   };
//
//   std::unordered_set<std::string, KeyHash> adduct_seen;
//   std::vector<std::string> filtered_ions;
//   std::vector<std::string> filtered_formulas;
//
//   for (int i = 0; i < adduct_ions.size(); ++i) {
//     std::string pol = tolower_copy(Rcpp::as<std::string>(adduct_pol[i]));
//     if (pol == "positive" || pol == "+" || pol == "pos") pol = "pos";
//     else if (pol == "negative" || pol == "-" || pol == "neg") pol = "neg";
//
//     if (pol != polarity) continue;
//
//     std::string ion  = Rcpp::as<std::string>(adduct_ions[i]);
//     std::string mass = Rcpp::as<std::string>(adduct_mass[i]);
//     std::string chg  = (adduct_charge[i] == NA_STRING) ? "" : Rcpp::as<std::string>(adduct_charge[i]);
//
//     std::string key = ion + "||" + mass + "||" + chg + "||" + pol;
//     if (adduct_seen.insert(key).second) {
//       filtered_ions.push_back(ion);
//       filtered_formulas.push_back(mass);
//     }
//   }
//
//   // ---- neutral loss table
//   NumericVector nl_mass = neutral_loss_table["Accurate.Mass"];
//   CharacterVector nl_name = neutral_loss_table.containsElementNamed("Formula")
//     ? (CharacterVector)neutral_loss_table["Formula"]
//   : (CharacterVector)neutral_loss_table["Neutral.Loss"];
//
//   CharacterVector nl_pos = neutral_loss_table["Pos"];
//   CharacterVector nl_neg = neutral_loss_table["Neg"];
//
//   // ---- output accumulators (row per matched feature)
//   std::vector<std::string> out_feature_name;
//   std::vector<std::string> out_lib_group;
//   std::vector<std::string> out_compound_name;
//   std::vector<std::string> out_ion_type;
//
//   std::vector<double> out_M;
//   std::vector<double> out_mz;
//   std::vector<double> out_rt;
//   std::vector<double> out_int;
//
//   std::vector<std::string> out_PI;
//   std::vector<std::string> out_adduct;
//   std::vector<std::string> out_NL;
//   std::vector<std::string> out_iso;
//
//   // optional group scores if present
//   std::vector<double> out_inner_cos, out_outer_cos, out_match_score;
//
//   int n_groups = library_matches.size();
//
//   for (int g = 0; g < n_groups; ++g) {
//     List grp = library_matches[g];
//
//     // --- group metadata
//     std::string compound_name = grp.containsElementNamed("Compound_name")
//       ? Rcpp::as<std::string>(grp["Compound_name"]) : "";
//
//     std::string lib_group = grp.containsElementNamed("Lib_group")
//       ? Rcpp::as<std::string>(grp["Lib_group"]) : "";
//
//     std::string ion_type = grp.containsElementNamed("Ion_type")
//       ? Rcpp::as<std::string>(grp["Ion_type"]) : "";
//
//     double M = grp.containsElementNamed("Compound_M")
//       ? Rcpp::as<double>(grp["Compound_M"]) : NA_REAL;
//
//     // vectors per group
//     CharacterVector fn = grp.containsElementNamed("Matched_features_to_peaks")
//       ? (CharacterVector)grp["Matched_features_to_peaks"] : CharacterVector();
//
//     NumericVector mz = grp.containsElementNamed("Matched_features_mz")
//       ? (NumericVector)grp["Matched_features_mz"] : NumericVector();
//
//     NumericVector rt = grp.containsElementNamed("Matched_features_rt")
//       ? (NumericVector)grp["Matched_features_rt"] : NumericVector(mz.size(), NA_REAL);
//
//     NumericVector inten = grp.containsElementNamed("Matched_features_int")
//       ? (NumericVector)grp["Matched_features_int"] : NumericVector(mz.size(), NA_REAL);
//
//     if (fn.size() == 0 || mz.size() == 0) continue;
//     if (mz.size() != fn.size()) stop("Matched_features_mz length != Matched_features_to_peaks length.");
//     if (rt.size() != mz.size()) rt = NumericVector(mz.size(), NA_REAL);
//     if (inten.size() != mz.size()) inten = NumericVector(mz.size(), NA_REAL);
//
//     // optional group scores
//     double inner_cos = grp.containsElementNamed("inner_cosine") ? Rcpp::as<double>(grp["inner_cosine"]) : NA_REAL;
//     double outer_cos = grp.containsElementNamed("outer_cosine") ? Rcpp::as<double>(grp["outer_cosine"]) : NA_REAL;
//     double match_sc  = grp.containsElementNamed("Match_peak_score") ? Rcpp::as<double>(grp["Match_peak_score"]) : NA_REAL;
//
//     // --- precompute theoretical adduct m/z for this M (once per group)
//     std::vector<double> theo_adduct_mz;
//     if (!NumericVector::is_na(M) && !filtered_formulas.empty()) {
//       theo_adduct_mz = parse_adduct_vector(M, filtered_formulas);
//     }
//
//     // --- PI tagging (optional): mark the closest feature to expected PI if within mz_threshold
//     double PI_target = NA_REAL;
//     std::string PI_tag = "";
//
//     if (!NumericVector::is_na(M)) {
//       std::string it = tolower_copy(ion_type);
//
//       // common variants: "M+H", "M+H+", "[M+H]+", etc.
//       if (it.find("m+h") != std::string::npos) {
//         PI_target = M + 1.007276;
//         PI_tag = " [M+H]+ ";
//       } else if (it.find("m-h") != std::string::npos) {
//         PI_target = M - 1.007276;
//         PI_tag = " [M-H]- ";
//       }
//     }
//
//     int pi_idx = -1;
//     if (!NumericVector::is_na(PI_target)) {
//       double best = 1e99;
//       for (int i = 0; i < mz.size(); ++i) {
//         double d = std::fabs(mz[i] - PI_target);
//         if (d < best) { best = d; pi_idx = i; }
//       }
//       if (best > mz_threshold) pi_idx = -1;
//     }
//
//     // --- isotope annotation within group (sort by m/z, run Check_Iso_Charge, map back)
//     std::vector<int> ord(mz.size());
//     std::iota(ord.begin(), ord.end(), 0);
//     std::sort(ord.begin(), ord.end(), [&](int a, int b){ return mz[a] < mz[b]; });
//
//     CharacterVector fn_s(fn.size());
//     NumericVector mz_s(mz.size()), rt_s(rt.size()), in_s(inten.size());
//     for (int k = 0; k < (int)ord.size(); ++k) {
//       fn_s[k] = fn[ord[k]];
//       mz_s[k] = mz[ord[k]];
//       rt_s[k] = rt[ord[k]];
//       in_s[k] = inten[ord[k]];
//     }
//
//     CharacterVector iso_s = Check_Iso_Charge(fn_s, mz_s, rt_s, in_s,
//                                              iso_threshold, max_charge, max_isotope);
//
//     std::vector<std::string> iso_label(mz.size(), "");
//     for (int k = 0; k < iso_s.size(); ++k) {
//       int orig = ord[k];
//       iso_label[orig] = Rcpp::as<std::string>(iso_s[k]);
//     }
//
//     // --- NL matching within group
//     std::vector<std::string> nl_label(mz.size(), "");
//     for (int i = 0; i < mz.size(); ++i) {
//       std::vector<std::string> hits;
//
//       for (int j = 0; j < mz.size(); ++j) {
//         if (i == j) continue;
//         if (mz[j] <= mz[i]) continue;
//
//         // RT constraint (only if both are available)
//         if (!NumericVector::is_na(rt[i]) && !NumericVector::is_na(rt[j])) {
//           if (std::fabs(rt[i] - rt[j]) > rt_threshold) continue;
//         }
//
//         double diff = mz[j] - mz[i];
//
//         for (int k = 0; k < nl_mass.size(); ++k) {
//           bool pol_ok = (polarity == "pos")
//           ? (Rcpp::as<std::string>(nl_pos[k]) == "+")
//             : (Rcpp::as<std::string>(nl_neg[k]) == "+");
//
//           if (!pol_ok) continue;
//
//           if (std::fabs(diff - nl_mass[k]) <= mz_threshold) {
//             std::string nm = Rcpp::as<std::string>(nl_name[k]);
//             std::string fj = Rcpp::as<std::string>(fn[j]);
//             hits.push_back("<- " + nm + " (" + fj + ")");
//           }
//         }
//       }
//       nl_label[i] = join_annotations(hits);
//     }
//
//     // --- adduct match (within group, using single neutral mass M)
//     std::vector<std::string> adduct_label(mz.size(), "");
//     if (!NumericVector::is_na(M) && !theo_adduct_mz.empty()) {
//       for (int i = 0; i < mz.size(); ++i) {
//         std::unordered_set<std::string, KeyHash> seen;
//         std::vector<std::string> hits;
//
//         size_t L = std::min(theo_adduct_mz.size(), filtered_ions.size());
//         for (size_t k = 0; k < L; ++k) {
//           double theo = theo_adduct_mz[k];
//           if (std::isnan(theo)) continue;
//           double d = std::fabs(mz[i] - theo);
//           double tol = mz_tol_mix(theo, mz_threshold, ppm);
//
//           if (d <= tol) {
//             std::ostringstream oss;
//             oss << filtered_ions[k] << " " << std::fixed << std::setprecision(6) << M;
//             std::string ann = oss.str();
//             if (seen.insert(ann).second) hits.push_back(ann);
//           }
//         }
//         adduct_label[i] = join_annotations(hits);
//       }
//     }
//
//     // --- append rows
//     for (int i = 0; i < mz.size(); ++i) {
//       out_feature_name.push_back(Rcpp::as<std::string>(fn[i]));
//       out_lib_group.push_back(lib_group);
//       out_compound_name.push_back(compound_name);
//       out_ion_type.push_back(ion_type);
//
//       out_M.push_back(M);
//       out_mz.push_back(mz[i]);
//       out_rt.push_back(rt[i]);
//       out_int.push_back(inten[i]);
//
//       // PI annotation only for selected index
//       if (i == pi_idx && !NumericVector::is_na(M)) {
//         std::ostringstream ms;
//         ms << std::fixed << std::setprecision(5) << M;
//         out_PI.push_back(compound_name + PI_tag + ms.str());
//       } else {
//         out_PI.push_back("");
//       }
//
//       out_adduct.push_back(adduct_label[i]);
//       out_NL.push_back(nl_label[i]);
//       out_iso.push_back(iso_label[i]);
//
//       out_inner_cos.push_back(inner_cos);
//       out_outer_cos.push_back(outer_cos);
//       out_match_score.push_back(match_sc);
//     }
//   }
//
//   return DataFrame::create(
//     _["feature_name"] = out_feature_name,
//     _["Lib_group"] = out_lib_group,
//     _["Compound_name"] = out_compound_name,
//     _["Ion_type"] = out_ion_type,
//     _["Compound_M"] = out_M,
//     _["mz"] = out_mz,
//     _["rt"] = out_rt,
//     _["intensity"] = out_int,
//     _["PI_match"] = out_PI,
//     _["adduct_match"] = out_adduct,
//     _["NL_match"] = out_NL,
//     _["iso_annotation"] = out_iso,
//     _["Match_peak_score"] = out_match_score,
//     _["inner_cosine"] = out_inner_cos,
//     _["outer_cosine"] = out_outer_cos,
//     _["stringsAsFactors"] = false
//   );
// }
//
//
