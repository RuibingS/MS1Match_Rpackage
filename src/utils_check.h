#ifndef UTILS_CHECK_H
#define UTILS_CHECK_H

#include <numeric>
#include <algorithm>
#include <vector>
#include <string>
#include <Rcpp.h>

// -----------------------------------------------------------------------------
// get_column_with_fallback
// Try primary column name first, then fallback, then stop with clear message.
// -----------------------------------------------------------------------------
template <typename T>
T get_column_with_fallback(Rcpp::DataFrame df,
                           std::string primary,
                           std::string fallback = "",
                           std::string label    = "") {
  if (df.containsElementNamed(primary.c_str())) return df[primary];
  if (!fallback.empty() && df.containsElementNamed(fallback.c_str())) return df[fallback];
  // FIX: build error label from whatever was actually provided
  if (label.empty()) {
    label = "'" + primary + "'";
    if (!fallback.empty()) label += " or '" + fallback + "'";
  }
  Rcpp::stop("None of the expected columns found: " + label);
}


// -----------------------------------------------------------------------------
// get_required_column — no fallback variant
// -----------------------------------------------------------------------------
template <typename T>
T get_required_column(Rcpp::DataFrame df, const char* name) {
  if (!df.containsElementNamed(name))
    Rcpp::stop("Missing required column: '" + std::string(name) + "'");
  return df[name];
}


// -----------------------------------------------------------------------------
// check_sorted_ascending — diagnostic only, does not sort
// -----------------------------------------------------------------------------
inline void check_sorted_ascending(Rcpp::NumericVector vec, std::string label) {
  for (int i = 1; i < vec.size(); ++i) {
    if (vec[i] < vec[i - 1])
      Rcpp::stop("Column '%s' must be sorted in ascending order.", label.c_str());
  }
}


// -----------------------------------------------------------------------------
// sort_dataframe_by_column
//
// Sorts all columns of df by the values in sort_col (ascending).
//
// FIX: row.names are reset to 1..n after sorting.
//      The original code copied the old row.names unchanged, so after sorting
//      row 3 would still be labelled "3" even if it was originally row 7.
// -----------------------------------------------------------------------------
inline Rcpp::DataFrame sort_dataframe_by_column(Rcpp::DataFrame df,
                                                Rcpp::NumericVector sort_col) {
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
      Rcpp::IntegerVector sv(n);
      for (int i = 0; i < n; ++i) sv[i] = v[idx[i]];
      sorted_df[j] = sv;
      break;
    }
    case REALSXP: {
      Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(col);
      Rcpp::NumericVector sv(n);
      for (int i = 0; i < n; ++i) sv[i] = v[idx[i]];
      sorted_df[j] = sv;
      break;
    }
    case STRSXP: {
      Rcpp::CharacterVector v = Rcpp::as<Rcpp::CharacterVector>(col);
      Rcpp::CharacterVector sv(n);
      for (int i = 0; i < n; ++i) sv[i] = v[idx[i]];
      sorted_df[j] = sv;
      break;
    }
    case LGLSXP: {
      Rcpp::LogicalVector v = Rcpp::as<Rcpp::LogicalVector>(col);
      Rcpp::LogicalVector sv(n);
      for (int i = 0; i < n; ++i) sv[i] = v[idx[i]];
      sorted_df[j] = sv;
      break;
    }
    default:
      Rcpp::stop("Unsupported column type in sort_dataframe_by_column: " +
                 std::string(colnames[j]));
    }
  }

  sorted_df.attr("names")     = colnames;
  sorted_df.attr("class")     = df.attr("class");
  // FIX: reset row names to sequential 1..n after sort
  sorted_df.attr("row.names") = Rcpp::seq(1, n);

  return Rcpp::DataFrame(sorted_df);
}

#endif
