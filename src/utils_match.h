#include <tuple>
#include <vector>
#include <cmath>

inline std::vector<std::tuple<int, int, double>> find_one_to_one_matches(
    const Rcpp::NumericVector& vec1,
    const Rcpp::NumericVector& vec2,
    double tol
) {
  std::vector<std::tuple<int, int, double>> matches;
  std::vector<bool> used1(vec1.size(), false);
  std::vector<bool> used2(vec2.size(), false);

  for (int i = 0; i < vec1.size(); ++i) {
    double best_diff = tol + 1;
    int best_j = -1;

    for (int j = 0; j < vec2.size(); ++j) {
      if (used2[j]) continue;
      double diff = std::abs(vec1[i] - vec2[j]);
      if (diff <= tol && diff < best_diff) {
        best_diff = diff;
        best_j = j;
      }
    }

    if (best_j >= 0) {
      matches.emplace_back(i, best_j, best_diff);
      used1[i] = true;
      used2[best_j] = true;
    }
  }

  return matches;
}
