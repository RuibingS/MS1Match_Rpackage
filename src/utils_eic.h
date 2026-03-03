#ifndef UTILS_EIC_H
#define UTILS_EIC_H

#include <Rcpp.h>

Rcpp::DataFrame closest_rt_fun(Rcpp::NumericVector rt1, Rcpp::NumericVector int1,
                               Rcpp::NumericVector rt2, Rcpp::NumericVector int2,
                               double rt_tol);

double calCor(Rcpp::NumericVector x, Rcpp::NumericVector y, std::string method);


Rcpp::List feature_EIC_fun(Rcpp::List input_list, int nrows, int ncols);

#endif
