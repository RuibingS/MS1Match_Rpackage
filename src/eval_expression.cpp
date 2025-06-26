#include <Rcpp.h>
#include <string>
#include <algorithm>

// [[Rcpp::export]]
double parse_and_eval_adduct(std::string expr, double M) {
  Rcpp::Rcout << "Original: [" << expr << "]" << std::endl;

  // Remove spaces
  expr.erase(remove(expr.begin(), expr.end(), ' '), expr.end());
  Rcpp::Rcout << "Cleaned: [" << expr << "]" << std::endl;

  // Find 'M'
  size_t m_pos = expr.find('M');
  if (m_pos == std::string::npos) {
    Rcpp::stop("No 'M' found in expression");
  }

  // Multiplier
  std::string before_M = expr.substr(0, m_pos);
  double multiplier = before_M.empty() ? 1.0 : std::stod(before_M);

  // After M
  std::string after_M = expr.substr(m_pos + 1);
  double divisor = 1.0;
  double offset = 0.0;

  // Find '+' or '-'
  size_t plus_pos = after_M.find('+');
  size_t minus_pos = after_M.find('-');
  size_t op_pos = (plus_pos != std::string::npos) ? plus_pos : minus_pos;
  char op = (plus_pos != std::string::npos) ? '+' : ((minus_pos != std::string::npos) ? '-' : '\0');

  if (!after_M.empty() && after_M[0] == '/') {
    if (op != '\0') {
      // Both divisor and offset
      divisor = std::stod(after_M.substr(1, op_pos - 1));
      offset = std::stod(after_M.substr(op_pos + 1));
      if (op == '-') offset = -offset;
    } else {
      // Only divisor
      divisor = std::stod(after_M.substr(1));
    }
  } else if (op != '\0') {
    // Only offset
    offset = std::stod(after_M.substr(op_pos + 1));
    if (op == '-') offset = -offset;
  }

  double result = (multiplier * M / divisor) + offset;
  Rcpp::Rcout << "Parsed: multiplier = " << multiplier
              << ", divisor = " << divisor
              << ", offset = " << offset << "\n";
  Rcpp::Rcout << "Result: " << result << "\n";

  return result;
}
