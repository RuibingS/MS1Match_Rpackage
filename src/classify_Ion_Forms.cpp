#include <Rcpp.h>
#include <set>
#include <cmath>
#include <algorithm>
#include "utils_annotation.h"
using namespace Rcpp;

// =============================================================================
// classify_Ion_Forms_library
//
// Annotates each peak in a spectral library against a table of known adducts,
// using the compound's exact neutral mass (ExactMass) to compute theoretical
// m/z values for each adduct, then matching observed peaks within tolerance.
//
// Match logic:
//   1. For each peak, find the adduct whose theoretical m/z is CLOSEST to the
//      observed m/z within ± tolerance (closest-match, not first-match).
//   2. If no adduct matches, classify as:
//        "Other smaller m/z" — below the mode reference ion ([M+H]+ or [M-H]-)
//        "Other larger m/z"  — above the mode reference ion
//
// Tolerance modes (tol_in_ppm):
//   tol_in_ppm = true  (default): tol is in ppm — tolerance scales with m/z
//                                 Recommended for QTOF / Orbitrap instruments.
//                                 For Bruker maXis HD QTOF: tol = 5 (ppm)
//                                 Formula: tol_Da = observed_mz * tol / 1e6
//   tol_in_ppm = false:           tol is in Da (absolute) — fixed window
//                                 Use for low-resolution instruments
//                                 (triple quad, ion trap): tol = 0.01–0.1 Da
//
// Parameters:
//   lib              — List of compounds from import_msp() / import_msp_csu().
//                      Each element must contain: Name, ExactMass, Peaks
//                      (data.frame with column "mz"), and optionally Ion_mode
//                      or IonPolarity.
//   adduct_table     — data.frame with columns:
//                        Ion.name      — adduct label e.g. "[M+H]+"
//                        Ion.mass      — expression string parsed by
//                                        parse_adduct_vector()
//                        Ion.polarity  — "pos" / "positive" / "+" or
//                                        "neg" / "negative" / "-"
//   lib_name         — Label for this library, written into the Library column.
//   filter_names     — Optional character vector of compound names to process.
//                      If NULL (default), all compounds are processed.
//   desired_polarity — Fallback polarity ("pos" or "neg") used when a compound
//                      has no Ion_mode / IonPolarity field. Default "pos".
//   tol              — Tolerance value. Interpreted as ppm if tol_in_ppm = true,
//                      or as Da if tol_in_ppm = false. Default 5 (ppm).
//   tol_in_ppm       — If true (default), tol is in ppm (recommended for QTOF).
//                      If false, tol is in Da (for low-resolution instruments).
//
// Returns:
//   data.frame with columns:
//     Compound    — compound name
//     mz          — observed m/z value from Peaks
//     Category    — matched adduct name, "Other smaller m/z", or
//                   "Other larger m/z"
//     Library     — value of lib_name
//     tol_used_da — actual Da tolerance used for each peak (useful for QC)
// =============================================================================
static inline bool is_positive(const std::string& s) {
  return s.rfind("pos", 0) == 0 || s == "+";
}
static inline bool is_negative(const std::string& s) {
  return s.rfind("neg", 0) == 0 || s == "-";
}

// [[Rcpp::export]]
DataFrame classify_Ion_Forms_library(
    List                      lib,
    DataFrame                 adduct_table,
    std::string               lib_name,
    Nullable<CharacterVector> filter_names     = R_NilValue,
    std::string               desired_polarity = "pos",
    double                    tol              = 5.0,   // default 5 ppm for QTOF
    bool                      tol_in_ppm       = true   // ppm mode by default
) {

  // ---------------------------------------------------------------------------
  // 1. Name filter
  // ---------------------------------------------------------------------------
  bool use_filter = false;
  std::set<std::string> filter_set;
  if (filter_names.isNotNull()) {
    CharacterVector fnames(filter_names);
    if (fnames.size() > 0) {
      use_filter = true;
      for (int i = 0; i < fnames.size(); ++i)
        filter_set.insert(as<std::string>(fnames[i]));
    }
  }

  // ---------------------------------------------------------------------------
  // 2. Adduct table
  // ---------------------------------------------------------------------------
  CharacterVector ion_names = adduct_table["Ion.name"];
  CharacterVector ion_exprs = adduct_table["Ion.mass"];
  CharacterVector ion_pols  = adduct_table["Ion.polarity"];
  int n_adducts = ion_names.size();

  std::vector<std::string> ion_pols_lower(n_adducts);
  for (int k = 0; k < n_adducts; ++k) {
    std::string p = as<std::string>(ion_pols[k]);
    std::transform(p.begin(), p.end(), p.begin(), ::tolower);
    ion_pols_lower[k] = p;
  }

  // ---------------------------------------------------------------------------
  // 3. Output vectors
  // ---------------------------------------------------------------------------
  std::vector<std::string> out_compounds;
  std::vector<std::string> out_categories;
  std::vector<std::string> out_libraries;
  std::vector<double>      out_mz;
  std::vector<double>      out_tol_da;

  // ---------------------------------------------------------------------------
  // 4. Main loop over library compounds
  // ---------------------------------------------------------------------------
  for (int i = 0; i < lib.size(); ++i) {
    List comp = lib[i];

    // --- compound name ---
    if (!comp.containsElementNamed("Name")) continue;
    std::string name = as<std::string>(comp["Name"]);

    // --- name filter ---
    if (use_filter && filter_set.find(name) == filter_set.end())
      continue;

    // --- ExactMass skip missing or NA  ---

    if (!comp.containsElementNamed("ExactMass")) continue;
    SEXP em_sexp = comp["ExactMass"];
    if (Rf_isNull(em_sexp)) continue;
    NumericVector em_vec = as<NumericVector>(em_sexp);
    if (em_vec.size() == 0) continue;
    double M = em_vec[0];
    if (ISNA(M) || std::isnan(M)) continue;

    // --- Peaks data.frame ---
    if (!comp.containsElementNamed("Peaks")) continue;
    DataFrame peaks = as<DataFrame>(comp["Peaks"]);
    if (!peaks.containsElementNamed("mz")) continue;
    NumericVector mz = peaks["mz"];
    if (mz.size() == 0) continue;

    // --- polarity: prefer compound field, fall back to desired_polarity ---
    std::string comp_pol = desired_polarity;
    if (comp.containsElementNamed("Ion_mode"))
      comp_pol = as<std::string>(comp["Ion_mode"]);
    else if (comp.containsElementNamed("IonPolarity"))
      comp_pol = as<std::string>(comp["IonPolarity"]);
    std::transform(comp_pol.begin(), comp_pol.end(), comp_pol.begin(), ::tolower);

    // --- select adducts matching compound polarity ---
    std::vector<std::string> adduct_exprs;
    std::vector<std::string> adduct_labels;

    for (int k = 0; k < n_adducts; ++k) {
      const std::string& pol = ion_pols_lower[k];
      if ((is_positive(comp_pol) && is_positive(pol)) ||
          (is_negative(comp_pol) && is_negative(pol))) {
        adduct_exprs.push_back(as<std::string>(ion_exprs[k]));
        adduct_labels.push_back(as<std::string>(ion_names[k]));
      }
    }

    if (adduct_exprs.empty()) continue;

    // --- compute theoretical m/z for each adduct ---
    std::vector<double> theo_mz = parse_adduct_vector(M, adduct_exprs);

    // --- mode-appropriate reference ion for "Other" classification ---
    // Positive: [M+H]+   fallback: M + 1.007276
    // Negative: [M-H]-   fallback: M - 1.007276
    double m_ref = NAN;
    if (is_positive(comp_pol)) {
      for (size_t k = 0; k < adduct_labels.size(); ++k) {
        if (adduct_labels[k].find("[M+H]+") != std::string::npos) {
          m_ref = theo_mz[k];
          break;
        }
      }
      if (std::isnan(m_ref)) m_ref = M + 1.007276;
    } else {
      for (size_t k = 0; k < adduct_labels.size(); ++k) {
        if (adduct_labels[k].find("[M-H]-") != std::string::npos) {
          m_ref = theo_mz[k];
          break;
        }
      }
      if (std::isnan(m_ref)) m_ref = M - 1.007276;
    }

    // --- classify each peak ---
    for (int j = 0; j < mz.size(); ++j) {
      double val = mz[j];
      if (ISNA(val) || std::isnan(val)) continue;

      // Compute tolerance in Da for this peak

      double tol_da = tol_in_ppm ? (val * tol / 1e6) : tol;

      std::string cat;
      bool   matched   = false;
      double best_diff = tol_da + 1.0;  // sentinel > tol_da

      // Closest-match within tolerance (not first-match)
      for (size_t k = 0; k < theo_mz.size(); ++k) {
        double t = theo_mz[k];
        if (std::isnan(t)) continue;
        double diff = std::fabs(val - t);
        if (diff <= tol_da && diff < best_diff) {
          best_diff = diff;
          cat       = adduct_labels[k];
          matched   = true;

        }
      }

      if (!matched) {
        cat = (val < m_ref) ? "Other smaller m/z" : "Other larger m/z";
      }

      out_compounds.push_back(name);
      out_mz.push_back(val);
      out_categories.push_back(cat);
      out_libraries.push_back(lib_name);
      out_tol_da.push_back(tol_da);
    }
  }

  return DataFrame::create(
    _["Compound"]         = out_compounds,
    _["mz"]               = out_mz,
    _["Category"]         = out_categories,
    _["Library"]          = out_libraries,
    _["tol_used_da"]      = out_tol_da,
    _["stringsAsFactors"] = false
  );
}
