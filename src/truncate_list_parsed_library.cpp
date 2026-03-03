#include <Rcpp.h>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "helpers.h"

using namespace Rcpp;

// ---------------- Core truncation function ----------------
void truncate_peaks(std::vector<double>& mzs,
                    std::vector<double>& intensities,
                    double mw_or_exact,
                    const std::string& ion_mode,
                    const std::string& adduct_type,
                    int max_isotope) {
  if (mw_or_exact <= 0) return;

  const double proton = 1.007276;
  const double sodium = 22.98977;
  const double isotope_shift = 1.0034;

  // Determine base m/z based on polarity and adduct type
  double base_mz = 0.0;
  std::string ion_mode_norm = normalize_ion_mode(ion_mode);

  if (ion_mode_norm == "NEGATIVE") {
    base_mz = mw_or_exact - proton;          // [M−H]− only
  } else {
    if (adduct_type == "[M+Na]+")
      base_mz = mw_or_exact + sodium;
    else
      base_mz = mw_or_exact + proton;        // default [M+H]+
  }

  // Compute cutoff (base m/z + isotopes)
  double cutoff = base_mz + max_isotope * isotope_shift + 0.05;

  std::vector<double> keep_mzs;
  std::vector<double> keep_ints;
  for (size_t i = 0; i < mzs.size(); ++i) {
    if (mzs[i] <= cutoff) {
      keep_mzs.push_back(mzs[i]);
      keep_ints.push_back(intensities[i]);
    }
  }

  mzs.swap(keep_mzs);
  intensities.swap(keep_ints);
}


// [[Rcpp::export]]
List truncate_parsed_library_cpp(List parsed_list,
                                 int max_isotope = 2,
                                 std::string adduct_type = "[M+Na]+",
                                 Nullable<std::string> export_path = R_NilValue) {
  std::vector<List> truncated;
  truncated.reserve(parsed_list.size());

  for (int i = 0; i < parsed_list.size(); i++) {
    List cmp = parsed_list[i];
    double mw_or_exact = NA_REAL;
    bool has_mass = false;

    if (cmp.containsElementNamed("ExactMass") && !Rf_isNull(cmp["ExactMass"])) {
      mw_or_exact = as<double>(cmp["ExactMass"]);
      has_mass = true;
    } else if (cmp.containsElementNamed("MW") && !Rf_isNull(cmp["MW"])) {
      mw_or_exact = as<double>(cmp["MW"]);
      has_mass = true;
    }

    std::string ion_mode = "POSITIVE";
    if (cmp.containsElementNamed("IonPolarity"))
      ion_mode = as<std::string>(cmp["IonPolarity"]);
    else if (cmp.containsElementNamed("Ion_mode"))
      ion_mode = as<std::string>(cmp["Ion_mode"]);

    if (!has_mass) {
      Rcpp::Rcout << "Skipping compound ";
      if (cmp.containsElementNamed("Name"))
        Rcpp::Rcout << as<std::string>(cmp["Name"]);
      else
        Rcpp::Rcout << "[unknown]";
      Rcpp::Rcout << ": missing ExactMass or MW.\n";
      truncated.push_back(cmp);
      continue;
    }

    if (!cmp.containsElementNamed("Peaks")) {
      Rcpp::Rcout << "Skipping compound: missing Peaks field.\n";
      truncated.push_back(cmp);
      continue;
    }

    DataFrame peaks = as<DataFrame>(cmp["Peaks"]);
    if (!peaks.containsElementNamed("mz") ||
        !peaks.containsElementNamed("rel_intensity")) {
        Rcpp::Rcout << "Skipping compound: missing mz or rel_intensity.\n";
      truncated.push_back(cmp);
      continue;
    }

    NumericVector mz_col = peaks["mz"];
    NumericVector int_col = peaks["rel_intensity"];

    std::vector<double> mzs(mz_col.begin(), mz_col.end());
    std::vector<double> intensities(int_col.begin(), int_col.end());

    truncate_peaks(mzs, intensities, mw_or_exact, ion_mode, adduct_type, max_isotope);

    DataFrame truncated_peaks = DataFrame::create(
      _["mz"] = NumericVector(mzs.begin(), mzs.end()),
      _["rel_intensity"] = NumericVector(intensities.begin(), intensities.end())
    );
    cmp["Peaks"] = truncated_peaks;
    truncated.push_back(cmp);
  }

  // ---------------- Optional MSP export ----------------
  if (export_path.isNotNull()) {
    std::ofstream outfile(as<std::string>(export_path));
    if (!outfile.is_open()) stop("Could not open export path.");

    for (auto &comp : truncated) {
      CharacterVector nms = comp.names();
      for (int j = 0; j < nms.size(); j++) {
        std::string nm = as<std::string>(nms[j]);
        if (nm == "Peaks") continue;

        SEXP val = comp[nm];
        if (Rf_isString(val) && Rf_length(val) == 1) {
          outfile << nm << ": " << as<std::string>(val) << "\n";
        } else if (Rf_isReal(val) && Rf_length(val) == 1) {
          outfile << nm << ": " << as<double>(val) << "\n";
        }
      }

      DataFrame df = as<DataFrame>(comp["Peaks"]);
      NumericVector mz_col = df["mz"];
      NumericVector int_col = df["rel_intensity"];

      outfile << "Num Peaks: " << df.nrows() << "\n";
      for (int i = 0; i < df.nrows(); i++) {
        outfile << mz_col[i] << "\t" << int_col[i] << "\n";
      }
      outfile << "\n";
    }

    outfile.close();
    Rcpp::Rcout << "Exported truncated MSP file: "
                << as<std::string>(export_path) << "\n";
  }

  return wrap(truncated);
}
