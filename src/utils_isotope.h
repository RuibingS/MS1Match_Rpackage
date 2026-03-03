#ifndef UTILS_ISOTOPE_H
#define UTILS_ISOTOPE_H

#include <Rcpp.h>

Rcpp::CharacterVector Check_Iso_Charge(Rcpp::CharacterVector name,
                                       Rcpp::NumericVector mz,
                                       Rcpp::NumericVector rt,
                                       Rcpp::NumericVector intensity,
                                       double iso_threshold,
                                       int max_charge,
                                       int max_isotope);

#endif
