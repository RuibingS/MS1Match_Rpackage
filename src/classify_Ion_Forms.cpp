#include <Rcpp.h>
#include <set>
#include <cmath>
#include <algorithm>
#include "utils_annotation.h"
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame classify_Ion_Forms_library(
    List lib,
    DataFrame adduct_table,
    std::string lib_name,
    Nullable<CharacterVector> filter_names = R_NilValue,
    std::string desired_polarity = "pos",
    double tol = 0.01 // ± tolerance
) {
  // --- name filter ---
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

  // --- output ---
  std::vector<std::string> compounds;
  std::vector<std::string> categories;
  std::vector<std::string> libraries;
  std::vector<double> mz_out;

  for (int i = 0; i < lib.size(); ++i) {
    List comp = lib[i];
    std::string name = as<std::string>(comp["Name"]);
    if (use_filter && filter_set.find(name) == filter_set.end())
      continue;

    double M = as<double>(comp["ExactMass"]);
    DataFrame peaks = as<DataFrame>(comp["Peaks"]);
    NumericVector mz = peaks["mz"];

    // --- polarity ---
    std::string comp_pol = desired_polarity;
    if (comp.containsElementNamed("IonPolarity"))
      comp_pol = as<std::string>(comp["IonPolarity"]);
    else if (comp.containsElementNamed("Ion_mode"))
      comp_pol = as<std::string>(comp["Ion_mode"]);
    std::transform(comp_pol.begin(), comp_pol.end(), comp_pol.begin(), ::tolower);

    // --- select adducts by polarity ---
    std::vector<std::string> exprs, names;
    for (int k = 0; k < n_adducts; ++k) {
      std::string pol = as<std::string>(ion_pols[k]);
      std::transform(pol.begin(), pol.end(), pol.begin(), ::tolower);
      if ((comp_pol.rfind("pos", 0) == 0 && pol.rfind("pos", 0) == 0) ||
          (comp_pol.rfind("neg", 0) == 0 && pol.rfind("neg", 0) == 0)) {
        exprs.push_back(as<std::string>(ion_exprs[k]));
        names.push_back(as<std::string>(ion_names[k]));
      }
    }

    if (exprs.empty()) continue;

    // --- theoretical m/z ---
    std::vector<double> theo_mz = parse_adduct_vector(M, exprs);

    // --- [M+Na]+ reference ---
    // double m_Na = NAN;
    // for (size_t k = 0; k < names.size(); ++k) {
    //   if (names[k].find("[M+Na]+") != std::string::npos) {
    //     m_Na = theo_mz[k];
    //     break;
    //   }
    // }
    // if (std::isnan(m_Na)) {
    //   m_Na = M + 22.9898; // compute theoretical [M+Na]+
    // }

    // --- [M+H]+ reference ---

    double m_proton = NAN;
    for (size_t k = 0; k < names.size(); ++k) {
      if (names[k].find("[M+H]+") != std::string::npos) {
        m_proton = theo_mz[k];
        break;
      }
    }
    if (std::isnan(m_proton)) {
      m_proton = M + 1.007276; // compute theoretical [M+H]+
    }

    // --- classify peaks ---
    for (int j = 0; j < mz.size(); ++j) {
      double val = mz[j];
      std::string cat;
      bool matched = false;

      for (size_t k = 0; k < theo_mz.size(); ++k) {
        double t = theo_mz[k];
        if (!std::isnan(t) && std::fabs(val - t) <= tol) {
          cat = names[k];
          matched = true;
          break;
        }
      }

      if (!matched) {
        if (val < m_proton)
          cat = "Other smaller m/z";
        else
          cat = "Other larger m/z";
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
