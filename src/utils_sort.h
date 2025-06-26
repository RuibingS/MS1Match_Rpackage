#ifndef UTILS_SORT_H
#define UTILS_SORT_H

#include <Rcpp.h>
#include <numeric>
#include <algorithm>

// Compute a sorted index vector
inline std::vector<int> compute_sorted_index(const Rcpp::NumericVector& v) {
  std::vector<int> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int i, int j) {
    return v[i] < v[j];
  });
  return idx;
}

// Reorder a numeric vector using a sorted index
inline Rcpp::NumericVector reorder_vector(const Rcpp::NumericVector& v, const std::vector<int>& idx) {
  Rcpp::NumericVector out(idx.size());
  for (int i = 0; i < idx.size(); ++i) {
    out[i] = v[idx[i]];
  }
  return out;
}

// Reorder a character vector using a sorted index
inline Rcpp::CharacterVector reorder_vector(const Rcpp::CharacterVector& v, const std::vector<int>& idx) {
  Rcpp::CharacterVector out(idx.size());
  for (int i = 0; i < idx.size(); ++i) {
    out[i] = v[idx[i]];
  }
  return out;
}

#endif
