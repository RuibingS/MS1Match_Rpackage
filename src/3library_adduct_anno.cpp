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

// [[Rcpp::export]]
DataFrame annotate_from_adduct_table(
    List lib,
    DataFrame adduct_table,
    std::string lib_name,
    Nullable<CharacterVector> filter_names = R_NilValue,
    std::string desired_polarity = "pos",
    double tol = 0.01
) {
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

  // --- adduct info ---
  CharacterVector ion_names = adduct_table["Ion.name"];
  CharacterVector ion_exprs = adduct_table["Ion.mass"];
  CharacterVector ion_pols  = adduct_table["Ion.polarity"];
  int n_adducts = ion_names.size();

  // --- output containers ---
  std::vector<std::string> compounds;
  std::vector<std::string> categories;
  std::vector<std::string> libraries;
  std::vector<double> mz_out;

  // --- iterate over compounds ---
  for (int i = 0; i < lib.size(); ++i) {
    List comp = lib[i];
    std::string name = as<std::string>(comp["Name"]);
    if (use_filter && filter_set.find(name) == filter_set.end()) continue;

    double M = as<double>(comp["ExactMass"]);
    DataFrame peaks = as<DataFrame>(comp["Peaks"]);
    NumericVector mz = peaks["mz"];

    // resolve polarity from compound metadata or fallback
    std::string comp_pol = get_str_if(comp, "IonPolarity");
    if (comp_pol.empty()) comp_pol = get_str_if(comp, "Ion_mode");
    if (comp_pol.empty()) comp_pol = desired_polarity;
    std::transform(comp_pol.begin(), comp_pol.end(), comp_pol.begin(), ::tolower);

    // subset adducts for matching polarity
    std::vector<int> keep_idx;
    for (int k = 0; k < n_adducts; ++k) {
      std::string pol = as<std::string>(ion_pols[k]);
      std::transform(pol.begin(), pol.end(), pol.begin(), ::tolower);
      if ((comp_pol.rfind("pos", 0) == 0 && pol.rfind("pos", 0) == 0) ||
          (comp_pol.rfind("neg", 0) == 0 && pol.rfind("neg", 0) == 0))
        keep_idx.push_back(k);
    }

    // compute adduct m/z for kept adducts
    std::vector<std::string> exprs, names;
    exprs.reserve(keep_idx.size());
    names.reserve(keep_idx.size());
    for (int idx : keep_idx) {
      exprs.push_back(as<std::string>(ion_exprs[idx]));
      names.push_back(as<std::string>(ion_names[idx]));
    }
    std::vector<double> theo_mz = parse_adduct_vector(M, exprs);

    // --- match peaks ---
    for (int j = 0; j < mz.size(); ++j) {
      double val = mz[j];
      std::string cat = "Undefined";
      for (size_t k = 0; k < theo_mz.size(); ++k) {
        double t = theo_mz[k];
        if (!std::isnan(t) && std::fabs(val - t) <= tol) {
          cat = names[k];
          break;
        }
      }
      compounds.push_back(name);
      categories.push_back(cat);
      libraries.push_back(lib_name);
      mz_out.push_back(val);
    }
  }

  return DataFrame::create(
    _["Compound"] = compounds,
    _["mz"]       = mz_out,
    _["Category"] = categories,
    _["Library"]  = libraries,
    _["stringsAsFactors"] = false
  );
}
