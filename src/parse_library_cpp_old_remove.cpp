#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

using namespace Rcpp;

// Helper: trim whitespace from both ends
// inline std::string trim(const std::string& s) {
//   auto start = s.begin();
//   while (start != s.end() && std::isspace(*start)) ++start;
//   auto end = s.end();
//   do {
//     --end;
//   } while (std::distance(start, end) > 0 && std::isspace(*end));
//   return std::string(start, end + 1);
// }



// ----------- String Trimming Helpers -----------
// trim from start (in place)
static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                  [](unsigned char ch) { return !std::isspace(ch); }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}

// trim and return a copy
static inline std::string trim_copy(std::string s) {
  trim(s);
  return s;
}

// ----------- Main Function ---------------
// [[Rcpp::export]]
Rcpp::List parse_library_cpp_old_remove(std::string file_path,
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
  std::vector<Rcpp::List> compounds;
  Rcpp::List compound;
  std::vector<double> mzs, intensities;

  bool in_peaks_section = false;
  int num_peaks_expected = 0, num_peaks_read = 0;
  bool has_msms = false;
  bool compound_valid = true;

  while (std::getline(infile, line)) {
    trim(line);
    if (line.empty()) continue;

    // New compound
    if (line.find("Name:") == 0 && compound.size() > 0) {
      // Handle filters
      std::string spectrum_filter = spectrum_type.isNotNull() ? Rcpp::as<std::string>(spectrum_type) : "";
      bool is_ms2_required = (spectrum_filter == "MS2");
      bool is_ms1_required = (spectrum_filter == "MS1");

      bool accept = compound_valid;
      if (spectrum_type.isNotNull()) {
        if (is_ms2_required && !has_msms) accept = false;
        if (is_ms1_required && has_msms)  accept = false;
      }

      if (accept && num_peaks_expected == num_peaks_read && num_peaks_read > 0) {
        Rcpp::DataFrame peaks_df = Rcpp::DataFrame::create(
          Rcpp::Named("mz") = Rcpp::NumericVector(mzs.begin(), mzs.end()),
          Rcpp::Named("rel_intensity") = Rcpp::NumericVector(intensities.begin(), intensities.end())
        );
        compound["Peaks"] = peaks_df;
        compounds.push_back(compound);
      }

      // Reset for next compound
      compound = Rcpp::List();
      mzs.clear();
      intensities.clear();
      num_peaks_expected = 0;
      num_peaks_read = 0;
      in_peaks_section = false;
      has_msms = false;
      compound_valid = true;
    }

    size_t pos = line.find(":");
    if (pos != std::string::npos) {
      std::string key = trim_copy(line.substr(0, pos));
      std::string value = trim_copy(line.substr(pos + 1));

      if (key == "MSMS") has_msms = true;

      // Optional filters
      if (key == "IonPolarity" && ionPolarity.isNotNull() &&
          value != Rcpp::as<std::string>(ionPolarity)) compound_valid = false;

      if (key == "IoniMethod" && ioniMethod.isNotNull() &&
          value != Rcpp::as<std::string>(ioniMethod)) compound_valid = false;

      if (key == "ColEnergy" && colEnergy.isNotNull() &&
          value != Rcpp::as<std::string>(colEnergy)) compound_valid = false;

      if (key == "InstName" && instName.isNotNull() &&
          value != Rcpp::as<std::string>(instName)) compound_valid = false;

      if (key == "Num Peaks") {
        num_peaks_expected = std::stoi(value);
        in_peaks_section = true;
        continue;
      }

      compound[key] = value;
    } else if (in_peaks_section) {
      std::istringstream iss(line);
      std::vector<double> values;
      double val;
      while (iss >> val) values.push_back(val);

      for (size_t i = 0; i + 1 < values.size(); i += 2) {
        mzs.push_back(values[i]);
        intensities.push_back(values[i + 1]);
        ++num_peaks_read;
      }
    }
  }

  // Handle last compound
  if (compound.size() > 0) {
    std::string spectrum_filter = spectrum_type.isNotNull() ? Rcpp::as<std::string>(spectrum_type) : "";
    bool is_ms2_required = (spectrum_filter == "MS2");
    bool is_ms1_required = (spectrum_filter == "MS1");

    bool accept = compound_valid;
    if (spectrum_type.isNotNull()) {
      if (is_ms2_required && !has_msms) accept = false;
      if (is_ms1_required && has_msms)  accept = false;
    }

    if (accept && num_peaks_expected == num_peaks_read && num_peaks_read > 0) {
      Rcpp::DataFrame peaks_df = Rcpp::DataFrame::create(
        Rcpp::Named("mz") = Rcpp::NumericVector(mzs.begin(), mzs.end()),
        Rcpp::Named("rel_intensity") = Rcpp::NumericVector(intensities.begin(), intensities.end())
      );
      compound["Peaks"] = peaks_df;
      compounds.push_back(compound);
    }
  }

  infile.close();
  return Rcpp::wrap(compounds);
}


// Rcpp::List parse_library_cpp(std::string file_path,
//                              Rcpp::Nullable<std::string> spectrum_type = R_NilValue,
//                              Rcpp::Nullable<std::string> ionPolarity = R_NilValue,
//                              Rcpp::Nullable<std::string> ioniMethod = R_NilValue,
//                              Rcpp::Nullable<std::string> colEnergy = R_NilValue,
//                              Rcpp::Nullable<std::string> instName = R_NilValue) {
//
//   std::ifstream infile(file_path);
//   if (!infile.is_open()) {
//     Rcpp::stop("Could not open file: " + file_path);
//   }
//
//   std::string line;
//   std::vector<Rcpp::List> compounds;
//   Rcpp::List compound;
//   std::vector<double> mzs, intensities;
//   bool in_peaks_section = false;
//   int num_peaks_expected = 0, num_peaks_read = 0;
//
//   while (std::getline(infile, line)) {
//     line = trim(line);
//     if (line.empty()) continue;
//
//     // Start of a new compound block
//     if (line.find("Name:") == 0 && compound.size() > 0) {
//
//       if (num_peaks_expected == num_peaks_read) {
//
//         Rcpp::NumericVector mz_col(mzs.begin(), mzs.end());
//         Rcpp::NumericVector intensity_col(intensities.begin(), intensities.end());
//
//         Rcpp::DataFrame peaks_df = Rcpp::DataFrame::create(
//           Rcpp::Named("mz") = mz_col,
//           Rcpp::Named("rel_intensity") = intensity_col
//         );
//
//         compound["Peaks"] = peaks_df;
//
//         compounds.push_back(compound);
//       }
//
//       // Reset for next
//       compound = Rcpp::List();
//       mzs.clear();
//       intensities.clear();
//       num_peaks_expected = 0;
//       num_peaks_read = 0;
//       in_peaks_section = false;
//     }
//
//     // Key-value pairs
//     size_t pos = line.find(":");
//     if (pos != std::string::npos) {
//       std::string key = trim(line.substr(0, pos));
//       std::string value = trim(line.substr(pos + 1));
//
//       // Filtering
//       if (key == "IonPolarity" && ionPolarity.isNotNull() && value != Rcpp::as<std::string>(ionPolarity)) {
//         compound = Rcpp::List();
//         continue;
//       }
//       if (key == "IoniMethod" && ioniMethod.isNotNull() && value != Rcpp::as<std::string>(ioniMethod)) {
//         compound = Rcpp::List();
//         continue;
//       }
//       if (key == "ColEnergy" && colEnergy.isNotNull() && value != Rcpp::as<std::string>(colEnergy)) {
//         compound = Rcpp::List();
//         continue;
//       }
//       if (key == "InstName" && instName.isNotNull() && value != Rcpp::as<std::string>(instName)) {
//         compound = Rcpp::List();
//         continue;
//       }
//       if (key == "MSMS" && spectrum_type.isNotNull() && value != Rcpp::as<std::string>(spectrum_type)) {
//         compound = Rcpp::List();
//         continue;
//       }
//
//       if (key == "Num Peaks") {
//         num_peaks_expected = std::stoi(value);
//         in_peaks_section = true;
//         continue;
//       }
//
//       // Add key-value pair
//       compound[key] = value;
//     }
//     else if (in_peaks_section) {
//       std::istringstream iss(line);
//       std::vector<double> values;
//       double val;
//       while (iss >> val) {
//         values.push_back(val);
//       }
//
//       // Convert to mz-intensity pairs
//       for (size_t i = 0; i + 1 < values.size(); i += 2) {
//         mzs.push_back(values[i]);
//         intensities.push_back(values[i + 1]);
//         ++num_peaks_read;
//       }
//     }
//
//     // else if (in_peaks_section) {
//     //   std::istringstream iss(line);
//     //   double mz, intensity;
//     //   while (iss >> mz >> intensity) {
//     //     mzs.push_back(mz);
//     //     intensities.push_back(intensity);
//     //     ++num_peaks_read;
//     //   }
//     // }
//   }
//
//   if (compound.size() > 0 && num_peaks_expected == num_peaks_read) {
//     Rcpp::NumericVector mz_col(mzs.begin(), mzs.end());
//     Rcpp::NumericVector intensity_col(intensities.begin(), intensities.end());
//
//     Rcpp::DataFrame peaks_df = Rcpp::DataFrame::create(
//       Rcpp::Named("mz") = mz_col,
//       Rcpp::Named("rel_intensity") = intensity_col
//     );
//     compound["Peaks"] = peaks_df;
//
//     compounds.push_back(compound);
//   }
//
//   infile.close();
//   return Rcpp::wrap(compounds);
// }
