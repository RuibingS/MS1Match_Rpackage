#ifndef UTILS_CHECK_H
#define UTILS_CHECK_H

#include <numeric>
#include <algorithm>
#include <vector>
#include <string>

#include <Rcpp.h>

// Check if required column exists (or fallback)
template <typename T>
T get_column_with_fallback(Rcpp::DataFrame df, std::string primary,
                           std::string fallback = "", std::string label = "") {
  if (df.containsElementNamed(primary.c_str())) return df[primary];
  if (!fallback.empty() && df.containsElementNamed(fallback.c_str())) return df[fallback];
  if (label.empty()) label = primary + "/" + fallback;
  Rcpp::stop("None of the expected fields found: '" + label + "'");
}


// Check if required column exists (no fallback)
template <typename T>
T get_required_column(Rcpp::DataFrame df, const char* name) {
  if (!df.containsElementNamed(name)) {
    Rcpp::stop("Missing required column: '" + std::string(name) + "'");
  }
  return df[name];
}

// Check if a numeric vector is sorted in ascending order
inline void check_sorted_ascending(Rcpp::NumericVector vec, std::string label) {
  for (int i = 1; i < vec.size(); ++i) {
    if (vec[i] < vec[i - 1]) {
      Rcpp::stop("Column '%s' must be sorted in ascending order.", label);
    }
  }
}


// Sort a data frame by a numeric column ascending
inline Rcpp::DataFrame sort_dataframe_by_column(Rcpp::DataFrame df, Rcpp::NumericVector sort_col) {
  int n = sort_col.size();
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);

  std::sort(idx.begin(), idx.end(), [&](int i, int j) {
    return sort_col[i] < sort_col[j];
  });

  Rcpp::List sorted_df(df.size());
  Rcpp::CharacterVector colnames = df.names();

  for (int j = 0; j < df.size(); ++j) {
    Rcpp::RObject col = df[j];

    switch (TYPEOF(col)) {
    case INTSXP: {
      Rcpp::IntegerVector v = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::IntegerVector sorted_v(n);
      for (int i = 0; i < n; ++i) sorted_v[i] = v[idx[i]];
      sorted_df[j] = sorted_v;
      break;
    }
    case REALSXP: {
      Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(col);
      Rcpp::NumericVector sorted_v(n);
      for (int i = 0; i < n; ++i) sorted_v[i] = v[idx[i]];
      sorted_df[j] = sorted_v;
      break;
    }
    case STRSXP: {
      Rcpp::CharacterVector v = Rcpp::as<Rcpp::CharacterVector>(col);
      Rcpp::CharacterVector sorted_v(n);
      for (int i = 0; i < n; ++i) sorted_v[i] = v[idx[i]];
      sorted_df[j] = sorted_v;
      break;
    }
    case LGLSXP: {
      Rcpp::LogicalVector v = Rcpp::as<Rcpp::LogicalVector>(col);
      Rcpp::LogicalVector sorted_v(n);
      for (int i = 0; i < n; ++i) sorted_v[i] = v[idx[i]];
      sorted_df[j] = sorted_v;
      break;
    }
    default:
      Rcpp::stop("Unsupported column type in FT: " + std::string(colnames[j]));
    }
  }

  sorted_df.attr("names") = colnames;
  sorted_df.attr("class") = df.attr("class");
  sorted_df.attr("row.names") = df.attr("row.names");

  return Rcpp::DataFrame(sorted_df);
}


#endif
