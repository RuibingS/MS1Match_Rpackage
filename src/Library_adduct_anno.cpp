#include <Rcpp.h>
#include <set>
#include <cmath>
#include <algorithm>
#include "utils_annotation.h"
using namespace Rcpp;

// helper: safely get string field
static inline std::string get_str_if(List x, const char* key) {
  return x.containsElementNamed(key) ? as<std::string>(x[key]) : std::string();
}

// ============================================================
// annotate_and_summarize
//
// For each compound in lib:
//   1. Annotate every peak against the adduct table (ion type)
//   2. Classify each peak relative to [M+H]+ cutoff:
//        "below_MH"  : mz <= [M+H]+ + tol  (normal fragments/precursor)
//        "above_MH"  : mz >  [M+H]+ + tol  (unexpected — adduct, contaminant, etc.)
//   3. Return two DataFrames as a List:
//        $peaks    — one row per peak: compound, mz, ion_type, above_MH flag
//        $summary  — one row per compound: counts by zone and ion type
// ============================================================

// [[Rcpp::export]]
List annotate_and_summarize(
    List        lib,
    DataFrame   adduct_table,
    std::string lib_name,
    Nullable<CharacterVector> filter_names = R_NilValue,
    std::string desired_polarity           = "pos",
    double      tol                        = 0.01,   // Da for adduct matching
    double      ppm_tol                    = 5.0     // ppm for [M+H]+ cutoff tolerance
) {
  const double proton     = 1.007276;
  const double ppm_factor = 1e-6;

  // --- name filter setup ---
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

  // --- adduct table ---
  CharacterVector ion_names = adduct_table["Ion.name"];
  CharacterVector ion_exprs = adduct_table["Ion.mass"];
  CharacterVector ion_pols  = adduct_table["Ion.polarity"];
  int n_adducts = ion_names.size();

  // --- peak-level output ---
  std::vector<std::string> pk_compound, pk_ion_type, pk_library;
  std::vector<double>      pk_mz, pk_mz_MH;
  std::vector<int>         pk_above_MH;   // 1 = above [M+H]+, 0 = at/below

  // --- summary-level output ---
  std::vector<std::string> sm_compound, sm_library;
  std::vector<double>      sm_exactmass, sm_mz_MH;
  std::vector<int>         sm_n_total, sm_n_above_MH, sm_n_below_MH;
  std::vector<int>         sm_n_above_annotated;    // above [M+H]+ AND matched adduct
  std::vector<int>         sm_n_above_undefined;    // above [M+H]+ AND unmatched

  // --- iterate over compounds ---
  for (int i = 0; i < lib.size(); ++i) {
    List comp = lib[i];

    if (!comp.containsElementNamed("Peaks") ||
        !comp.containsElementNamed("ExactMass") ||
        !comp.containsElementNamed("Name"))
        continue;

    std::string name      = as<std::string>(comp["Name"]);
    if (use_filter && filter_set.find(name) == filter_set.end()) continue;

    double        M       = as<double>(comp["ExactMass"]);
    DataFrame     peaks   = as<DataFrame>(comp["Peaks"]);
    NumericVector mz      = peaks["mz"];
    int           total   = mz.size();

    // [M+H]+ cutoff
    double base_mz = M + proton;
    double tol_abs = base_mz * ppm_tol * ppm_factor;

    // resolve polarity
    std::string comp_pol = get_str_if(comp, "IonPolarity");
    if (comp_pol.empty()) comp_pol = get_str_if(comp, "Ion_mode");
    if (comp_pol.empty()) comp_pol = desired_polarity;
    std::transform(comp_pol.begin(), comp_pol.end(), comp_pol.begin(), ::tolower);

    // subset adducts matching polarity
    std::vector<int> keep_idx;
    for (int k = 0; k < n_adducts; ++k) {
      std::string pol = as<std::string>(ion_pols[k]);
      std::transform(pol.begin(), pol.end(), pol.begin(), ::tolower);
      if ((comp_pol.rfind("pos", 0) == 0 && pol.rfind("pos", 0) == 0) ||
          (comp_pol.rfind("neg", 0) == 0 && pol.rfind("neg", 0) == 0))
        keep_idx.push_back(k);
    }

    // compute theoretical adduct m/z
    std::vector<std::string> exprs, anames;
    exprs.reserve(keep_idx.size());
    anames.reserve(keep_idx.size());
    for (int idx : keep_idx) {
      exprs.push_back(as<std::string>(ion_exprs[idx]));
      anames.push_back(as<std::string>(ion_names[idx]));
    }
    std::vector<double> theo_mz = parse_adduct_vector(M, exprs);

    // per-compound counters
    int c_above = 0, c_below = 0;
    int c_above_annotated = 0, c_above_undefined = 0;

    // --- per-peak annotation ---
    for (int j = 0; j < total; ++j) {
      double val    = mz[j];
      bool   above  = (val > base_mz + tol_abs);

      // match against adduct table
      std::string ion_type = "Undefined";
      for (size_t k = 0; k < theo_mz.size(); ++k) {
        double t = theo_mz[k];
        if (!std::isnan(t) && std::fabs(val - t) <= tol) {
          ion_type = anames[k];
          break;
        }
      }

      // peak-level record
      pk_compound.push_back(name);
      pk_mz.push_back(val);
      pk_mz_MH.push_back(base_mz);
      pk_ion_type.push_back(ion_type);
      pk_above_MH.push_back(above ? 1 : 0);
      pk_library.push_back(lib_name);

      // summary counters
      if (above) {
        c_above++;
        if (ion_type == "Undefined") c_above_undefined++;
        else                         c_above_annotated++;
      } else {
        c_below++;
      }
    }

    // summary-level record
    sm_compound.push_back(name);
    sm_library.push_back(lib_name);
    sm_exactmass.push_back(M);
    sm_mz_MH.push_back(base_mz);
    sm_n_total.push_back(total);
    sm_n_above_MH.push_back(c_above);
    sm_n_below_MH.push_back(c_below);
    sm_n_above_annotated.push_back(c_above_annotated);
    sm_n_above_undefined.push_back(c_above_undefined);
  }

  // --- build output ---
  DataFrame peaks_df = DataFrame::create(
    _["Compound"]   = pk_compound,
    _["mz"]         = pk_mz,
    _["mz_MH"]      = pk_mz_MH,       // [M+H]+ reference for this compound
    _["ion_type"]   = pk_ion_type,     // matched adduct name or "Undefined"
    _["above_MH"]   = pk_above_MH,     // 1 if peak is above [M+H]+
    _["Library"]    = pk_library,
    _["stringsAsFactors"] = false
  );

  DataFrame summary_df = DataFrame::create(
    _["Compound"]            = sm_compound,
    _["ExactMass"]           = sm_exactmass,
    _["mz_MH"]               = sm_mz_MH,
    _["n_total"]             = sm_n_total,
    _["n_below_MH"]          = sm_n_below_MH,
    _["n_above_MH"]          = sm_n_above_MH,
    _["n_above_annotated"]   = sm_n_above_annotated,  // above + matched ion type
    _["n_above_undefined"]   = sm_n_above_undefined,  // above + no match = suspicious
    _["Library"]             = sm_library,
    _["stringsAsFactors"]    = false
  );

  return List::create(
    _["peaks"]   = peaks_df,
    _["summary"] = summary_df
  );
}
