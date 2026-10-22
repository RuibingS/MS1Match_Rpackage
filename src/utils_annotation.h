#ifndef UTILS_ANNOTATION_H
#define UTILS_ANNOTATION_H

#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cmath>

inline std::string join_annotations(const std::vector<std::string>& annotations,
                                    const std::string& delimiter = "; ") {
  std::string result;
  for (size_t i = 0; i < annotations.size(); ++i) {
    result += annotations[i];
    if (i < annotations.size() - 1) result += delimiter;
  }
  return result;
}


// -----------------------------------------------------------------------------
// parse_single_adduct
//
// Parses an adduct Ion.mass formula string, e.g.:
//   "M + 1.007276"       -> multiplier=1, divisor=1, offset=+1.007276
//   "2M + 22.989218"     -> multiplier=2, divisor=1, offset=+22.989218
//   "M/2 + 1.007276"     -> multiplier=1, divisor=2, offset=+1.007276
//   "4M/2 + 11.998247"   -> multiplier=4, divisor=2, offset=+11.998247
//   "3M/2 - 1.007276"    -> multiplier=3, divisor=2, offset=-1.007276
//   "M - 17.003289"      -> multiplier=1, divisor=1, offset=-17.003289
//   "4M - 1.007276"      -> multiplier=4, divisor=1, offset=-1.007276
//
// Returns: (multiplier * M / divisor) + offset
//
// Fixes vs original:
//   1. std::remove instead of unqualified remove (portability)
//   2. When string starts with '/', +/- search starts AFTER the divisor
//      digits — prevents mistaking a '-' inside "3M/2-1.007" as op_pos=1
//      (which would give substr(1,0) = "" -> stod throws)
//   3. Guard against size_t underflow when op_pos == 0
//   4. Guard against empty substrings before calling stod()
//   5. Zero-divisor check
// -----------------------------------------------------------------------------
inline double parse_single_adduct(const std::string& expr, double M) {
  std::string e = expr;

  // FIX 1: qualified std::remove
  e.erase(std::remove(e.begin(), e.end(), ' '), e.end());

  size_t m_pos = e.find('M');
  if (m_pos == std::string::npos)
    throw std::runtime_error("No 'M' found in expression: " + expr);

  std::string before_M = e.substr(0, m_pos);
  std::string after_M  = e.substr(m_pos + 1);

  // --- multiplier (everything before 'M')
  double multiplier = 1.0;
  if (!before_M.empty()) {
    try { multiplier = std::stod(before_M); }
    catch (...) {
      throw std::runtime_error("Cannot parse multiplier '" + before_M +
                               "' in: " + expr);
    }
  }

  double divisor = 1.0;
  double offset  = 0.0;

  if (!after_M.empty() && after_M[0] == '/') {
    // Form: /D  or  /D+O  or  /D-O
    // FIX 2: scan for +/- starting at index 2 (skip '/' and first digit char)
    // so we never mistake a sign that is part of the divisor number itself.
    size_t op_pos = std::string::npos;
    char   op     = '\0';
    for (size_t p = 2; p < after_M.size(); ++p) {
      if (after_M[p] == '+' || after_M[p] == '-') {
        op_pos = p;
        op     = after_M[p];
        break;
      }
    }

    if (op != '\0') {
      // /D OP O
      // FIX 3 & 4: op_pos >= 2 guaranteed by loop above, no underflow risk
      std::string div_str = after_M.substr(1, op_pos - 1);
      std::string off_str = after_M.substr(op_pos + 1);

      if (div_str.empty())
        throw std::runtime_error("Empty divisor in: " + expr);
      if (off_str.empty())
        throw std::runtime_error("Empty offset in: " + expr);

      try { divisor = std::stod(div_str); }
      catch (...) {
        throw std::runtime_error("Cannot parse divisor '" + div_str +
                                 "' in: " + expr);
      }
      try { offset = std::stod(off_str); }
      catch (...) {
        throw std::runtime_error("Cannot parse offset '" + off_str +
                                 "' in: " + expr);
      }
      if (op == '-') offset = -offset;

    } else {
      // /D only
      std::string div_str = after_M.substr(1);
      if (div_str.empty())
        throw std::runtime_error("Empty divisor in: " + expr);
      try { divisor = std::stod(div_str); }
      catch (...) {
        throw std::runtime_error("Cannot parse divisor '" + div_str +
                                 "' in: " + expr);
      }
    }

  } else if (!after_M.empty()) {
    // Form: +O  or  -O  (no divisor)
    char        op      = after_M[0];
    std::string off_str = after_M.substr(1);

    if ((op == '+' || op == '-') && !off_str.empty()) {
      try { offset = std::stod(off_str); }
      catch (...) {
        throw std::runtime_error("Cannot parse offset '" + off_str +
                                 "' in: " + expr);
      }
      if (op == '-') offset = -offset;
    }
  }

  // FIX 5: zero divisor guard
  if (divisor == 0.0)
    throw std::runtime_error("Zero divisor in: " + expr);

  return (multiplier * M / divisor) + offset;
}


// Vectorised: compute observed m/z for each adduct formula given neutral mass M
inline std::vector<double> parse_adduct_vector(double M,
                                               const std::vector<std::string>& adducts) {
  std::vector<double> result;
  result.reserve(adducts.size());
  for (const auto& expr : adducts) {
    try {
      result.push_back(parse_single_adduct(expr, M));
    } catch (...) {
      result.push_back(NAN);
    }
  }
  return result;
}

#endif
