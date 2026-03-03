#ifndef HELPERS_H
#define HELPERS_H

#include <string>
#include <algorithm>
#include <cctype>

static inline void trim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                  [](unsigned char ch) { return !std::isspace(ch); }));
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
}

static inline std::string trim_copy(std::string s) { trim(s); return s; }

static inline std::string to_lower(const std::string &s) {
  std::string out = s;
  std::transform(out.begin(), out.end(), out.begin(), ::tolower);
  return out;
}

static inline std::string normalize_spectrum_type(const std::string &s) {
  std::string val = to_lower(s);
  if (val == "ms" || val == "ms1") return "MS1";
  if (val == "ms2") return "MS2";
  return s;
}

static inline std::string normalize_ion_mode(const std::string &s) {
  std::string val = to_lower(s);
  if (val == "positive" || val == "pos") return "POSITIVE";
  if (val == "negative" || val == "neg") return "NEGATIVE";
  return s;
}

#endif
