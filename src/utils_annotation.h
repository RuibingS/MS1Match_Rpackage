#ifndef UTILS_ANNOTATION_H
#define UTILS_ANNOTATION_H

#include <string>
#include <vector>

#include <algorithm>
#include <stdexcept>
#include <sstream>

inline std::string join_annotations(const std::vector<std::string>& annotations, const std::string& delimiter = "; ") {
  std::string result;
  for (size_t i = 0; i < annotations.size(); ++i) {
    result += annotations[i];
    if (i < annotations.size() - 1) {
      result += delimiter;
    }
  }
  return result;
}



// Helper to parse and evaluate a single adduct string e.g.
inline double parse_single_adduct(const std::string& expr, double M) {
  std::string e = expr;
  e.erase(remove(e.begin(), e.end(), ' '), e.end());  // remove whitespace

  size_t m_pos = e.find('M');
  if (m_pos == std::string::npos)
    throw std::runtime_error("No 'M' found in expression: " + expr);

  std::string before_M = e.substr(0, m_pos);
  std::string after_M = e.substr(m_pos + 1);

  double multiplier = before_M.empty() ? 1.0 : std::stod(before_M);
  double divisor = 1.0;
  double offset = 0.0;

  size_t plus_pos = after_M.find('+');
  size_t minus_pos = after_M.find('-');
  size_t op_pos = (plus_pos != std::string::npos) ? plus_pos :
    (minus_pos != std::string::npos ? minus_pos : std::string::npos);
  char op = (plus_pos != std::string::npos) ? '+' :
    (minus_pos != std::string::npos ? '-' : '\0');

  if (!after_M.empty() && after_M[0] == '/') {
    if (op != '\0') {
      divisor = std::stod(after_M.substr(1, op_pos - 1));
      offset = std::stod(after_M.substr(op_pos + 1));
      if (op == '-') offset = -offset;
    } else {
      divisor = std::stod(after_M.substr(1));
    }
  } else if (op != '\0') {
    offset = std::stod(after_M.substr(op_pos + 1));
    if (op == '-') offset = -offset;
  }

  return (multiplier * M / divisor) + offset;
}

// Vectorized adduct mz
inline std::vector<double> parse_adduct_vector(double M, const std::vector<std::string>& adducts) {
  std::vector<double> result;
  result.reserve(adducts.size());

  for (const auto& expr : adducts) {
    try {
      result.push_back(parse_single_adduct(expr, M));
    } catch (const std::exception& e) {
      result.push_back(NAN);
    }
  }

  return result;
}


#endif
