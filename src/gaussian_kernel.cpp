#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP gaussian_kernel_dist(NumericVector rt_values,
                          double sigma,
                          Nullable<DataFrame> feature_pairs = R_NilValue) {
  double denom = 2 * sigma * sigma;

  // full matrix
  if (feature_pairs.isNull()) {
    int n = rt_values.size();
    NumericMatrix dist_matrix(n, n);
    NumericMatrix sim_matrix(n, n);

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double diff = rt_values[i] - rt_values[j];
        double dist = std::fabs(diff);
        dist_matrix(i, j) = dist;
        sim_matrix(i, j)  = std::exp(-(dist * dist) / denom);
      }
    }

    return List::create(
      _["dist_matrix"] = dist_matrix,
      _["sim_matrix"]  = sim_matrix
    );
  }

  // pairwise features
  DataFrame pairs = feature_pairs.get();
  CharacterVector f1 = pairs["feature_1"];
  CharacterVector f2 = pairs["feature_2"];
  CharacterVector feature_names = pairs.attr("feature_names");

  int n_pairs = f1.size();

  // map name -> index map
  std::unordered_map<std::string, int> feature_index;
  for (int i = 0; i < feature_names.size(); ++i) {
    feature_index[Rcpp::as<std::string>(feature_names[i])] = i;
  }

  std::vector<double> sim_vec(n_pairs);

  for (int i = 0; i < n_pairs; ++i) {
    auto it1 = feature_index.find(Rcpp::as<std::string>(f1[i]));
    auto it2 = feature_index.find(Rcpp::as<std::string>(f2[i]));

    if (it1 == feature_index.end() || it2 == feature_index.end()) {
      sim_vec[i] = NA_REAL;
      continue;
    }

    double diff = rt_values[it1->second] - rt_values[it2->second];
    double dist = std::fabs(diff);
    sim_vec[i]  = std::exp(-(dist * dist) / denom);
  }

  return DataFrame::create(
    _["feature_1"]     = f1,
    _["feature_2"]     = f2,
    _["rt_similarity"] = sim_vec
  );
}
