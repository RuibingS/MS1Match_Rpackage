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
// get intensity_mean
inline Rcpp::NumericVector get_intensity_from_FT(Rcpp::DataFrame FT) {
  const char* candidates[] = {"intensity_mean", "intensity", "Inten", "I"};
  for (auto& name : candidates) {
    if (FT.containsElementNamed(name)) return FT[name];
  }
  Rcpp::stop("No intensity column found. Expected one of: intensity_mean, intensity, Inten, I.");
}

// Reorder a numeric vector using a sorted index
inline Rcpp::NumericVector reorder_vector(const Rcpp::NumericVector& v, const std::vector<int>& idx) {
  Rcpp::NumericVector out(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) {
    out[i] = v[idx[i]];
  }
  return out;
}

// Reorder a character vector using a sorted index
inline Rcpp::CharacterVector reorder_vector(const Rcpp::CharacterVector& v, const std::vector<int>& idx) {
  Rcpp::CharacterVector out(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) {
    out[i] = v[idx[i]];
  }
  return out;
}

// Reorder a list using sorted index
inline Rcpp::List reorder_list(const Rcpp::List& lst, const std::vector<int>& idx) {
  Rcpp::List out(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) {
    out[i] = lst[idx[i]];
  }
  return out;
}



#endif
