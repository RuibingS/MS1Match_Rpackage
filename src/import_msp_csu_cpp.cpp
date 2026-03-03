#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
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

// ---------- Save compound ----------
void save_compound(std::vector<Rcpp::List> &compounds,
                   Rcpp::List &compound,
                   std::vector<double> &mzs,
                   std::vector<double> &intensities,
                   bool compound_valid,
                   bool has_msms,
                   int num_peaks_read) {
  if (!compound_valid || num_peaks_read == 0)
    return;

  Rcpp::DataFrame peaks_df = Rcpp::DataFrame::create(
    Rcpp::Named("mz") = Rcpp::NumericVector(mzs.begin(), mzs.end()),
    Rcpp::Named("rel_intensity") = Rcpp::NumericVector(intensities.begin(), intensities.end())
  );

  compound["Peaks"] = peaks_df;
  compound["NumPeaks"] = num_peaks_read;
  compounds.push_back(compound);
}

// ---------- Main parser for CSU–PMF ----------
// [[Rcpp::export]]
Rcpp::List parse_library_CSU_cpp(std::string file_path) {

  std::ifstream infile(file_path);
  if (!infile.is_open())
    Rcpp::stop("Could not open file: " + file_path);

  std::string line;
  std::vector<Rcpp::List> compounds;
  Rcpp::List compound;
  std::vector<double> mzs, intensities;

  bool in_peaks_section = false;
  int num_peaks_read = 0;
  bool has_msms = false;
  bool compound_valid = true;
  std::string formula = "";

  int total_compounds = 0;
  int filtered_out = 0;

  while (std::getline(infile, line)) {
    trim(line);
    if (line.empty()) continue;

    // Start of a new compound
    if (line.rfind("Name:", 0) == 0 && compound.size() > 0) {
      total_compounds++;

      // filter invalid Formula
      if (to_lower_copy(formula) == "unk" || formula.empty()) {
        filtered_out++;
        compound_valid = false;
      }

      save_compound(compounds, compound, mzs, intensities,
                    compound_valid, has_msms, num_peaks_read);

      // reset
      compound = Rcpp::List();
      mzs.clear();
      intensities.clear();
      in_peaks_section = false;
      num_peaks_read = 0;
      has_msms = false;
      compound_valid = true;
      formula.clear();
    }

    // Parse key:value
    size_t pos = line.find(":");
    if (pos != std::string::npos) {
      std::string key = line.substr(0, pos);
      std::string value = line.substr(pos + 1);
      trim(key); trim(value);

      if (key == "Formula") {
        formula = value;
        if (to_lower_copy(value) == "unk" || value.empty())
          compound_valid = false;
      }

      if (key == "Num Peaks") {
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

  // Save last compound
  if (compound.size() > 0) {
    total_compounds++;
    if (to_lower_copy(formula) == "unk" || formula.empty()) {
      filtered_out++;
      compound_valid = false;
    }
    save_compound(compounds, compound, mzs, intensities,
                  compound_valid, has_msms, num_peaks_read);
  }
  infile.close();

  // ---------- Deduplicate by Name ----------
  std::map<std::string, int> best_idx;
  std::vector<Rcpp::List> deduped;
  int duplicates_removed = 0;

  for (size_t i = 0; i < compounds.size(); ++i) {
    Rcpp::List cmp = compounds[i];
    if (!cmp.containsElementNamed("Name") || !cmp.containsElementNamed("NumPeaks"))
      continue;

    std::string name = Rcpp::as<std::string>(cmp["Name"]);
    int np = Rcpp::as<int>(cmp["NumPeaks"]);

    if (best_idx.find(name) == best_idx.end()) {
      best_idx[name] = i;
    } else {
      duplicates_removed++;
      int existing_idx = best_idx[name];
      Rcpp::List existing = compounds[existing_idx];
      int np_existing = existing.containsElementNamed("NumPeaks") ?
      Rcpp::as<int>(existing["NumPeaks"]) : 0;
      if (np > np_existing) {
        best_idx[name] = i;
      }
    }
  }

  for (auto &kv : best_idx)
    deduped.push_back(compounds[kv.second]);

  // ---------- Summary print ----------
  Rcpp::Rcout << "[CSU_PMF parser] Total compounds parsed: " << total_compounds << "\n"
              << "[CSU_PMF parser] Filtered out (Formula: unk or invalid): " << filtered_out << "\n"
              << "[CSU_PMF parser] Duplicates removed: " << duplicates_removed << "\n"
              << "[CSU_PMF parser] Unique compounds retained: " << deduped.size() << "\n";

  return Rcpp::wrap(deduped);
}
