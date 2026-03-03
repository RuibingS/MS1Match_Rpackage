#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cctype>

using namespace Rcpp;

// ----------- Helpers ---------------
static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                  [](unsigned char ch) { return !std::isspace(ch); }));
}
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}
static inline void trim(std::string &s) {
  ltrim(s); rtrim(s);
}
static inline std::string trim_copy(std::string s) {
  trim(s); return s;
}

std::string to_lower(const std::string &s) {
  std::string out = s;
  std::transform(out.begin(), out.end(), out.begin(), ::tolower);
  return out;
}

// normalize spectrum type
std::string normalize_spectrum_type(const std::string &s) {
  std::string val = to_lower(s);
  if (val == "ms" || val == "ms1") return "MS1";
  if (val == "ms2") return "MS2";
  return s;
}

// normalize ion mode
std::string normalize_ion_mode(const std::string &s) {
  std::string val = to_lower(s);
  if (val == "positive" || val == "pos") return "POSITIVE";
  if (val == "negative" || val == "neg") return "NEGATIVE";
  return s; // fallback
}

void save_compound(std::unordered_map<std::string, Rcpp::List>& best_compounds,
                   std::unordered_map<std::string, int>& best_peak_count,
                   Rcpp::List& compound,
                   std::vector<double>& mzs,
                   std::vector<double>& intensities,
                   bool compound_valid,
                   bool has_msms,
                   Rcpp::Nullable<std::string> spectrum_type,
                   int num_peaks_expected,
                   int num_peaks_read,
                   int& total_count,
                   int& filtered_count,
                   int& duplicate_count) {

  if (compound.size() == 0) return;

  total_count++;

  std::string compound_name = compound.containsElementNamed("Name") ?
  Rcpp::as<std::string>(compound["Name"]) : "";

  std::string spectrum_filter = spectrum_type.isNotNull() ?
  normalize_spectrum_type(Rcpp::as<std::string>(spectrum_type)) : "";

  bool accept = compound_valid;
  if (spectrum_type.isNotNull()) {
    if (spectrum_filter == "MS2" && !has_msms) accept = false;
    if (spectrum_filter == "MS1" && has_msms)  accept = false;
  }

  if (!accept) {
    filtered_count++;
    return;
  }

  if (num_peaks_expected == num_peaks_read && num_peaks_read > 0) {
    double max_intensity = *std::max_element(intensities.begin(), intensities.end());
    std::vector<double> rel_intensities(intensities.size(), 0.0);
    if (max_intensity > 0.0) {
      for (size_t i = 0; i < intensities.size(); ++i)
        rel_intensities[i] = (intensities[i] / max_intensity) * 100.0;
    }

    Rcpp::DataFrame peaks_df = Rcpp::DataFrame::create(
      Rcpp::Named("mz") = Rcpp::NumericVector(mzs.begin(), mzs.end()),
      Rcpp::Named("rel_intensity") = Rcpp::NumericVector(rel_intensities.begin(), rel_intensities.end())
    );
    compound["Peaks"] = peaks_df;

    if (!best_compounds.count(compound_name)) {
      best_compounds[compound_name] = compound;
      best_peak_count[compound_name] = num_peaks_read;
    } else {
      duplicate_count++;
      if (num_peaks_read > best_peak_count[compound_name]) {
        best_compounds[compound_name] = compound;
        best_peak_count[compound_name] = num_peaks_read;
      }
    }
  }
}


// ----------- Main Function ---------------
// [[Rcpp::export]]
Rcpp::List parse_msp_library_cpp(std::string file_path,
                                 Rcpp::Nullable<std::string> spectrum_type = R_NilValue,
                                 Rcpp::Nullable<std::string> ionPolarity = R_NilValue,
                                 Rcpp::Nullable<std::string> ioniMethod = R_NilValue,
                                 Rcpp::Nullable<std::string> colEnergy = R_NilValue,
                                 Rcpp::Nullable<std::string> instName = R_NilValue) {

  std::ifstream infile(file_path);
  if (!infile.is_open()) {
    Rcpp::stop("Could not open file: " + file_path);
  }

  std::string line;
  std::unordered_map<std::string, Rcpp::List> best_compounds;
  std::unordered_map<std::string, int> best_peak_count;

  Rcpp::List compound;
  std::vector<double> mzs, intensities;

  bool in_peaks_section = false;
  int num_peaks_expected = 0, num_peaks_read = 0;
  bool has_msms = false;
  bool compound_valid = true;

  int total_count = 0;
  int filtered_count = 0;
  int duplicate_count = 0;

  while (std::getline(infile, line)) {
    trim(line);
    if (line.empty()) continue;

    if (line.find("Name:") == 0 && compound.size() > 0) {
      save_compound(best_compounds, best_peak_count, compound, mzs, intensities,
                    compound_valid, has_msms, spectrum_type,
                    num_peaks_expected, num_peaks_read,
                    total_count, filtered_count, duplicate_count);

      compound = Rcpp::List();
      mzs.clear();
      intensities.clear();
      num_peaks_expected = num_peaks_read = 0;
      in_peaks_section = false;
      has_msms = false;
      compound_valid = true;
    }

    size_t pos = line.find(":");
    if (pos != std::string::npos) {
      std::string key = trim_copy(line.substr(0, pos));
      std::string value = trim_copy(line.substr(pos + 1));
      std::string key_lower = to_lower(key);

      if (key_lower == "msms") has_msms = true;

      // filter conditions
      if (key_lower == "spectrum_type") {
        std::string norm_val = normalize_spectrum_type(value);
        compound["Spectrum_type"] = norm_val;
        if (spectrum_type.isNotNull()) {
          std::string expected = normalize_spectrum_type(Rcpp::as<std::string>(spectrum_type));
          if (expected != norm_val) compound_valid = false;
        }
      } else if (key_lower == "ion_mode" || key_lower == "ionpolarity") {
        std::string norm_val = normalize_ion_mode(value);
        compound["Ion_mode"] = norm_val;
        if (ionPolarity.isNotNull()) {
          std::string expected = normalize_ion_mode(Rcpp::as<std::string>(ionPolarity));
          if (expected != norm_val) compound_valid = false;
        }
      } else if (key_lower == "instrument_type" && ioniMethod.isNotNull()) {
        if (value != Rcpp::as<std::string>(ioniMethod)) compound_valid = false;
        compound[key] = value;
      } else if (key_lower == "collision_energy" && colEnergy.isNotNull()) {
        if (value != Rcpp::as<std::string>(colEnergy)) compound_valid = false;
        compound[key] = value;
      } else if (key_lower == "instrument" && instName.isNotNull()) {
        if (value != Rcpp::as<std::string>(instName)) compound_valid = false;
        compound[key] = value;
      } else if (key_lower == "exactmass" || key_lower == "mw") {
        try {
          compound[key] = std::stod(value);
        } catch (...) {
          compound[key] = NA_REAL;
        }
      } else if (key_lower == "num peaks") {
        num_peaks_expected = std::stoi(value);
        in_peaks_section = true;
        compound[key] = num_peaks_expected;
        continue;
      } else {
        compound[key] = value;
      }
    } else if (in_peaks_section) {
      std::istringstream iss(line);
      double mz, intensity;
      if (iss >> mz >> intensity) {
        mzs.push_back(mz);
        intensities.push_back(intensity);
        ++num_peaks_read;
      }
    }
  }

  save_compound(best_compounds, best_peak_count, compound, mzs, intensities,
                compound_valid, has_msms, spectrum_type,
                num_peaks_expected, num_peaks_read,
                total_count, filtered_count, duplicate_count);

  infile.close();

  // Convert map to list
  Rcpp::List result;
  for (auto &p : best_compounds) result.push_back(p.second);

  int n_final = result.size();

  // ---- Print summary ----
  Rcpp::Rcout << "[MoNA parser] Total compounds parsed: " << total_count << std::endl;
  Rcpp::Rcout << "[MoNA parser] Filtered out (mismatch or invalid): " << filtered_count << std::endl;
  Rcpp::Rcout << "[MoNA parser] Duplicates removed: " << duplicate_count << std::endl;
  Rcpp::Rcout << "[MoNA parser] Unique compounds retained: " << n_final << std::endl;

  return result;
}
