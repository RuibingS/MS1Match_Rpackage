#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double test_eic_similarity_df(DataFrame eic_df1, DataFrame eic_df2,
                              double rt_thr = 2.0) {
  // Extract columns
  NumericVector rt1 = eic_df1["rt"];
  NumericVector int1 = eic_df1["intensity"];
  NumericVector rt2 = eic_df2["rt"];
  NumericVector int2 = eic_df2["intensity"];

  std::vector<double> aligned1;
  std::vector<double> aligned2;

  // For each point in rt1, find nearest point in rt2
  for (int i = 0; i < rt1.size(); ++i) {
    double best_diff = rt_thr + 1.0;
    int best_j = -1;

    for (int j = 0; j < rt2.size(); ++j) {
      double diff = std::fabs(rt1[i] - rt2[j]);
      if (diff < best_diff && diff <= rt_thr &&
          !NumericVector::is_na(int1[i]) && !NumericVector::is_na(int2[j]) &&
          int1[i] > 0 && int2[j] > 0) {
        best_diff = diff;
        best_j = j;
      }
    }

    if (best_j >= 0) {
      aligned1.push_back(std::log10(int1[i]));
      aligned2.push_back(std::log10(int2[best_j]));
    }
  }

  if (aligned1.size() < 3) {
    Rcout << "Too few matched points: " << aligned1.size() << std::endl;
    return NA_REAL;
  }

  // Compute Pearson correlation
  NumericVector x = wrap(aligned1);
  NumericVector y = wrap(aligned2);

  double mean_x = mean(x);
  double mean_y = mean(y);

  double num = 0, den_x = 0, den_y = 0;
  for (int i = 0; i < x.size(); ++i) {
    double dx = x[i] - mean_x;
    double dy = y[i] - mean_y;
    num += dx * dy;
    den_x += dx * dx;
    den_y += dy * dy;
  }

  return num / std::sqrt(den_x * den_y);
}
