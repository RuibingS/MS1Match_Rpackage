// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double pearson_cpp(NumericVector x, NumericVector y) {
  arma::vec a = as<arma::vec>(x);
  arma::vec b = as<arma::vec>(y);
  return arma::as_scalar(arma::cor(a, b));
}

// [[Rcpp::export]]
double spearman_via_r(NumericVector x, NumericVector y) {
  Function cor("cor");
  return as<double>(cor(x, y, Named("method") = "spearman"));
}

// [[Rcpp::export]]
double kendall_via_r(NumericVector x, NumericVector y) {
  Function cor("cor");
  return as<double>(cor(x, y, Named("method") = "kendall"));
}
