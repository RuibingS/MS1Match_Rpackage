#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>
#include "helpers.h"

using namespace Rcpp;

// ----------- Truncate helper with adduct type ---------------
void truncate_ms1(std::vector<double>& mzs,
                  std::vector<double>& intensities,
                  double mw_or_exact,
                  const std::string& ion_mode,
                  const std::string& adduct_type,
                  int max_isotope) {
  if (mw_or_exact <= 0) return;

  const double proton = 1.007276;
  const double sodium = 22.98977;
  const double isotope_shift = 1.0034;

  // Normalize polarity
  std::string ion_mode_norm = normalize_ion_mode(ion_mode);
  double base_mz = 0.0;

  if (ion_mode_norm == "NEGATIVE") {
    base_mz = mw_or_exact - proton; // [M−H]− only
  } else {
    if (adduct_type == "[M+Na]+")
      base_mz = mw_or_exact + sodium;
    else
      base_mz = mw_or_exact + proton; // default [M+H]+
  }

  double cutoff = base_mz + max_isotope * isotope_shift + 0.1;

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

// ----------- Main Function ---------------
// [[Rcpp::export]]
Rcpp::List parse_any_library_cpp(std::string file_path,
                                 bool truncateMS1 = false,
                                 double ppm_tol = 20.0,
                                 int max_isotope = 2,
                                 std::string adduct_type = "[M+Na]+",
                                 Rcpp::Nullable<std::string> spectrum_type = R_NilValue,
                                 Rcpp::Nullable<std::string> ionPolarity   = R_NilValue,
                                 Rcpp::Nullable<std::string> export_path  = R_NilValue) {

  std::ifstream infile(file_path);
  if (!infile.is_open()) Rcpp::stop("Could not open file: " + file_path);

  std::string line;
  std::vector<Rcpp::List> compounds;
  Rcpp::List compound;
  std::vector<double> mzs, intensities;

  int num_peaks_expected = 0, num_peaks_read = 0;
  bool in_peaks_section = false;
  bool is_library_format = false;

  auto finalize_compound = [&](void) {
    if (compound.size() == 0) return;
    if (num_peaks_read > 0) {
      double mw_or_exact = NA_REAL;
      if (compound.containsElementNamed("ExactMass"))
        mw_or_exact = Rcpp::as<double>(compound["ExactMass"]);
      else if (compound.containsElementNamed("MW"))
        mw_or_exact = Rcpp::as<double>(compound["MW"]);

      std::string ion_mode = compound.containsElementNamed("Ion_mode") ?
      Rcpp::as<std::string>(compound["Ion_mode"]) : "";

      std::string spectrum = compound.containsElementNamed("Spectrum_type") ?
      Rcpp::as<std::string>(compound["Spectrum_type"]) : "MS1";

      bool accept = true;
      if (spectrum_type.isNotNull()) {
        std::string expected = normalize_spectrum_type(Rcpp::as<std::string>(spectrum_type));
        if (expected != spectrum) accept = false;
      }
      if (ionPolarity.isNotNull()) {
        std::string expected = normalize_ion_mode(Rcpp::as<std::string>(ionPolarity));
        if (expected != ion_mode) accept = false;
      }

      if (truncateMS1 && spectrum == "MS1") {
        truncate_ms1(mzs, intensities, mw_or_exact, ion_mode, adduct_type, max_isotope);
      }

      Rcpp::DataFrame peaks_df = Rcpp::DataFrame::create(
        Rcpp::Named("mz") = Rcpp::NumericVector(mzs.begin(), mzs.end()),
        Rcpp::Named("rel_intensity") = Rcpp::NumericVector(intensities.begin(), intensities.end())
      );
      compound["Peaks"] = peaks_df;
      compounds.push_back(compound);
    }
  };

  while (std::getline(infile, line)) {
    trim(line);
    if (line.empty()) continue;

    if (line.find("Name:") == 0 && compound.size() > 0) {
      finalize_compound();
      compound = Rcpp::List();
      mzs.clear(); intensities.clear();
      num_peaks_expected = num_peaks_read = 0;
      in_peaks_section = false;
      is_library_format = false;
    }

    size_t pos = line.find(":");
    if (pos != std::string::npos) {
      std::string key = trim_copy(line.substr(0, pos));
      std::string value = trim_copy(line.substr(pos + 1));
      std::string key_lower = to_lower(key);

      if (key_lower == "spectrum_type") {
        compound["Spectrum_type"] = normalize_spectrum_type(value);
      } else if (key_lower == "ion_mode" || key_lower == "ionpolarity") {
        compound["Ion_mode"] = normalize_ion_mode(value);
      } else if (key_lower == "exactmass") {
        try { compound["ExactMass"] = std::stod(value); }
        catch (...) { compound["ExactMass"] = NA_REAL; }
      } else if (key_lower == "mw") {
        try { compound["MW"] = std::stod(value); }
        catch (...) { compound["MW"] = NA_REAL; }
      } else if (key_lower == "num peaks") {
        num_peaks_expected = std::stoi(value);

        // Peek next line for format detection
        std::streampos pos_now = infile.tellg();
        std::string next_line;
        std::getline(infile, next_line);
        trim(next_line);
        std::istringstream iss(next_line);
        std::vector<double> tmp_vals;
        double v;
        while (iss >> v) tmp_vals.push_back(v);

        if (tmp_vals.size() > 2) {
          is_library_format = true; // all peaks in one line
          for (size_t i = 0; i + 1 < tmp_vals.size(); i += 2) {
            mzs.push_back(tmp_vals[i]);
            intensities.push_back(tmp_vals[i + 1]);
            ++num_peaks_read;
          }
        } else {
          is_library_format = false;
          infile.seekg(pos_now);
          in_peaks_section = true;
        }
      } else {
        compound[key] = value;
      }
    } else if (in_peaks_section) {
      std::istringstream iss(line);
      double mz, rel_intensity;
      if (iss >> mz >> rel_intensity) {
        mzs.push_back(mz);
        intensities.push_back(rel_intensity);
        ++num_peaks_read;
        if (num_peaks_read >= num_peaks_expected) in_peaks_section = false;
      }
    }
  }
  infile.close();
  finalize_compound();

  // ---------- Optional export ----------
  if (export_path.isNotNull()) {
    std::ofstream outfile(Rcpp::as<std::string>(export_path));
    for (auto &comp : compounds) {
      Rcpp::CharacterVector nms = comp.names();
      for (int j = 0; j < nms.size(); j++) {
        std::string nm = Rcpp::as<std::string>(nms[j]);
        if (nm == "Peaks") continue;

        SEXP val = comp[nm];
        if (Rf_isString(val) && Rf_length(val) == 1) {
          outfile << nm << ": " << Rcpp::as<std::string>(val) << "\n";
        } else if (Rf_isReal(val) && Rf_length(val) == 1) {
          outfile << nm << ": " << Rcpp::as<double>(val) << "\n";
        } else if (Rf_isInteger(val) && Rf_length(val) == 1) {
          outfile << nm << ": " << Rcpp::as<int>(val) << "\n";
        }
      }

      Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(comp["Peaks"]);
      Rcpp::NumericVector mz_col  = df["mz"];
      Rcpp::NumericVector int_col = df["rel_intensity"];

      outfile << "Num Peaks: " << df.nrows() << "\n";
      for (int i = 0; i < df.nrows(); i++) {
        outfile << mz_col[i] << "\t" << int_col[i] << "\n";
      }
      outfile << "\n";
    }
    outfile.close();
  }

  return Rcpp::wrap(compounds);
}
