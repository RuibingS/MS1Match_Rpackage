#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cctype>

using namespace Rcpp;

// ---------- String utilities ----------
static inline void trim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                  [](unsigned char ch) { return !std::isspace(ch); }));
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); }).base(),
                       s.end());
}

static inline std::string to_lower_copy(const std::string &s) {
  std::string res = s;
  std::transform(res.begin(), res.end(), res.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return res;
}

// ---------- Helper to save one compound ----------
void save_compound(std::vector<Rcpp::List> &compounds,
                   Rcpp::List &compound,
                   std::vector<double> &mzs,
                   std::vector<double> &intensities,
                   bool compound_valid,
                   bool has_msms,
                   int num_peaks_expected,
                   int num_peaks_read,
                   Rcpp::Nullable<std::string> spectrum_type,
                   int &total_count) {

  if (compound.size() == 0)
    return;

  total_count++;

  if (!compound_valid || num_peaks_read == 0)
    return;

  std::string spectrum_filter =
    spectrum_type.isNotNull() ? Rcpp::as<std::string>(spectrum_type) : "";
  bool is_ms2_required = (spectrum_filter == "MS2");
  bool is_ms1_required = (spectrum_filter == "MS1");

  bool accept = true;
  if (is_ms2_required && !has_msms)
    accept = false;
  if (is_ms1_required && has_msms)
    accept = false;

  if (!accept)
    return;

  Rcpp::DataFrame peaks_df = Rcpp::DataFrame::create(
    Rcpp::Named("mz") = Rcpp::NumericVector(mzs.begin(), mzs.end()),
    Rcpp::Named("rel_intensity") =
      Rcpp::NumericVector(intensities.begin(), intensities.end()));

  compound["Peaks"] = peaks_df;
  compound["NumPeaks"] = num_peaks_read;

  compounds.push_back(compound);
}

// ---------- Main Function ----------
// [[Rcpp::export]]
Rcpp::List parse_library_cpp(std::string file_path,
                             Rcpp::Nullable<std::string> spectrum_type = R_NilValue,
                             Rcpp::Nullable<std::string> ionPolarity = R_NilValue,
                             Rcpp::Nullable<std::string> ioniMethod = R_NilValue,
                             Rcpp::Nullable<std::string> colEnergy = R_NilValue,
                             Rcpp::Nullable<std::string> instName = R_NilValue) {

  std::ifstream infile(file_path);
  if (!infile.is_open())
    Rcpp::stop("Could not open file: " + file_path);

  std::string line;
  std::vector<Rcpp::List> compounds;
  Rcpp::List compound;
  std::vector<double> mzs, intensities;

  bool in_peaks_section = false;
  int num_peaks_expected = 0, num_peaks_read = 0;
  bool has_msms = false;
  bool compound_valid = true;
  int total_count = 0;

  while (std::getline(infile, line)) {
    trim(line);
    if (line.empty())
      continue;

    // --- New compound starts ---
    if (line.rfind("Name:", 0) == 0 && compound.size() > 0) {
      save_compound(compounds, compound, mzs, intensities, compound_valid, has_msms,
                    num_peaks_expected, num_peaks_read, spectrum_type, total_count);

      // reset
      compound = Rcpp::List();
      mzs.clear();
      intensities.clear();
      num_peaks_expected = 0;
      num_peaks_read = 0;
      in_peaks_section = false;
      has_msms = false;
      compound_valid = true;
    }

    // --- Parse key:value ---
    size_t pos = line.find(":");
    if (pos != std::string::npos) {
      std::string key = line.substr(0, pos);
      std::string value = line.substr(pos + 1);
      trim(key);
      trim(value);

      if (key == "MSMS")
        has_msms = true;

      // --- Case-insensitive IonPolarity filter ---
      if (key == "IonPolarity" && ionPolarity.isNotNull()) {
        std::string filter = to_lower_copy(Rcpp::as<std::string>(ionPolarity));
        std::string val = to_lower_copy(value);
        if (!(val.find(filter) != std::string::npos))
          compound_valid = false;
      }

      if (key == "IoniMethod" && ioniMethod.isNotNull() &&
          value != Rcpp::as<std::string>(ioniMethod))
        compound_valid = false;

      if (key == "ColEnergy" && colEnergy.isNotNull() &&
          value != Rcpp::as<std::string>(colEnergy))
        compound_valid = false;

      if (key == "InstName" && instName.isNotNull() &&
          value != Rcpp::as<std::string>(instName))
        compound_valid = false;

      if (key == "Num Peaks") {
        num_peaks_expected = std::stoi(value);
        in_peaks_section = true;
        continue;
      }

      compound[key] = value;
    } else if (in_peaks_section) {
      std::istringstream iss(line);
      double m, inten;
      while (iss >> m >> inten) {
        mzs.push_back(m);
        intensities.push_back(inten);
        ++num_peaks_read;
      }
    }
  }

  // --- Save last compound ---
  save_compound(compounds, compound, mzs, intensities, compound_valid, has_msms,
                num_peaks_expected, num_peaks_read, spectrum_type, total_count);

  infile.close();

  // ---------- Deduplicate by Name ----------
  std::map<std::string, int> best_idx;
  std::vector<Rcpp::List> deduped;
  int duplicate_count = 0;

  for (size_t i = 0; i < compounds.size(); ++i) {
    Rcpp::List cmp = compounds[i];
    if (!cmp.containsElementNamed("Name") ||
        !cmp.containsElementNamed("NumPeaks"))
        continue;

    std::string name = Rcpp::as<std::string>(cmp["Name"]);
    int np = Rcpp::as<int>(cmp["NumPeaks"]);

    if (best_idx.find(name) == best_idx.end()) {
      best_idx[name] = i;
    } else {
      duplicate_count++;
      int existing_idx = best_idx[name];
      Rcpp::List existing = compounds[existing_idx];
      int np_existing = existing.containsElementNamed("NumPeaks")
        ? Rcpp::as<int>(existing["NumPeaks"])
          : 0;
      if (np > np_existing) {
        best_idx[name] = i;
      }
    }
  }

  for (auto &kv : best_idx)
    deduped.push_back(compounds[kv.second]);

  int n_final = deduped.size();

  // ---------- Print summary ----------
  Rcpp::Rcout << "[In-house parser] Total compounds parsed: " << total_count << std::endl;
  Rcpp::Rcout << "[In-house parser] Duplicates removed: " << duplicate_count << std::endl;
  Rcpp::Rcout << "[In-house parser] Unique compounds retained: " << n_final << std::endl;

  return Rcpp::wrap(deduped);
}
