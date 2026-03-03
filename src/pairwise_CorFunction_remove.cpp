#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double calCorrelation(NumericVector col1, NumericVector col2, std::string method) {
  LogicalVector completeIndices = !is_na(col1) & !is_na(col2);

  if (sum(completeIndices) > 1) {
    arma::vec v1 = as<arma::vec>(col1[completeIndices]);
    arma::vec v2 = as<arma::vec>(col2[completeIndices]);

    if (method == "pearson") {
      return arma::as_scalar(arma::cor(v1, v2));
    } else if (method == "spearman" || method == "kendall") {

      Function rank("rank");
      NumericVector rank1 = rank(wrap(v1));
      NumericVector rank2 = rank(wrap(v2));
      Function cor("cor");
      return as<double>(cor(rank1, rank2, Named("method") = "pearson"));
    } else {
      stop("Unsupported correlation method. Use 'pearson', 'spearman', or 'kendall'.");
    }
  } else {
    return NA_REAL;
  }
}


// [[Rcpp::export]]
NumericMatrix pairwiseCor_mat(NumericMatrix x, std::string method) {
  int nCols = x.ncol();
  NumericMatrix result(nCols, nCols);

  for (int i = 0; i < nCols - 1; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      NumericVector col1 = x(_, i);
      NumericVector col2 = x(_, j);
      result(i, j) = calCorrelation(col1, col2, method);
      result(j, i) = result(i, j);  // Symmetric
    }
  }

  for (int i = 0; i < nCols; ++i) {
    result(i, i) = 1.0;
  }

  return result;
}


// [[Rcpp::export]]
DataFrame pairwiseCor_long(NumericMatrix x, CharacterVector feature_names, std::string method) {
  int nCols = x.ncol();

  std::vector<std::string> feature1;
  std::vector<std::string> feature2;
  std::vector<double> correlations;

  for (int i = 0; i < nCols - 1; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      NumericVector col1 = x(_, i);
      NumericVector col2 = x(_, j);
      double corr = calCorrelation(col1, col2, method);

      feature1.push_back(as<std::string>(feature_names[i]));
      feature2.push_back(as<std::string>(feature_names[j]));
      correlations.push_back(corr);
    }
  }

  return DataFrame::create(
    _["feature_1"] = feature1,
    _["feature_2"] = feature2,
    _["cor_intensity"] = correlations
  );
}
