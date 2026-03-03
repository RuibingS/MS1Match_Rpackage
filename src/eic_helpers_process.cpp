#include <RcppArmadillo.h>
#include <unordered_map>
#include <cmath>
using namespace Rcpp;
#include "utils_eic.h"
using namespace Rcpp;



// double calCor(NumericVector col1, NumericVector col2, std::string method = "pearson") {
//   // Treat zero as missing
//   for (int i = 0; i < col1.size(); ++i) {
//     if (col1[i] <= 0 || col2[i] <= 0) {
//       col1[i] = NA_REAL;
//       col2[i] = NA_REAL;
//     }
//   }
//
//   LogicalVector complete = !is_na(col1) & !is_na(col2);
//   if (sum(complete) > 1) {
//     arma::vec v1 = as<arma::vec>(col1[complete]);
//     arma::vec v2 = as<arma::vec>(col2[complete]);
//
//     // Apply log10 safely
//     v1 = arma::log10(v1);
//     v2 = arma::log10(v2);
//
//     if (method == "pearson") {
//       return arma::as_scalar(arma::cor(v1, v2));
//     } else if (method == "spearman" || method == "kendall") {
//       Function rank("rank");
//       NumericVector rank1 = rank(wrap(v1));
//       NumericVector rank2 = rank(wrap(v2));
//       Function cor("cor");
//       return as<double>(cor(rank1, rank2, Named("method") = method));
//     } else {
//       stop("Unsupported correlation method. Use 'pearson', 'spearman', or 'kendall'.");
//     }
//   }
//   return NA_REAL;
// }


// [[Rcpp::export]]
DataFrame closest_rt_fun(NumericVector rt1, NumericVector int1,
                         NumericVector rt2, NumericVector int2, double rt_EIC_thr) {
  std::vector<double> m_rt1, m_int1, m_rt2, m_int2;

  for (int i = 0; i < rt1.size(); ++i) {
    double best_diff = rt_EIC_thr + 1.0;
    double best_rt2 = NA_REAL, best_int2 = NA_REAL;

    for (int j = 0; j < rt2.size(); ++j) {
      double diff = std::abs(rt1[i] - rt2[j]);
      if (diff <= rt_EIC_thr && diff < best_diff) {
        best_diff = diff;
        best_rt2 = rt2[j];
        best_int2 = int2[j];
      }
    }

    if (!NumericVector::is_na(best_rt2)) {
      m_rt1.push_back(rt1[i]);
      m_int1.push_back(int1[i]);
      m_rt2.push_back(best_rt2);
      m_int2.push_back(best_int2);
    }
  }

  return DataFrame::create(
    Named("rt1") = m_rt1,
    Named("int1") = m_int1,
    Named("rt2") = m_rt2,
    Named("int2") = m_int2
  );
}

// DataFrame closest_rt_fun(NumericVector rt1, NumericVector int1,
//                          NumericVector rt2, NumericVector int2, double rt_EIC_thr) {
//   int n1 = rt1.size(), n2 = rt2.size();
//   std::vector<double> m_rt1, m_int1, m_rt2, m_int2;
//   int j = 0;
//   for (int i = 0; i < n1; ++i) {
//     while (j < n2 && rt2[j] < rt1[i] - rt_EIC_thr) j++;
//     if (j < n2 && std::abs(rt1[i] - rt2[j]) <= rt_EIC_thr) {
//       m_rt1.push_back(rt1[i]);
//       m_int1.push_back(int1[i]);
//       m_rt2.push_back(rt2[j]);
//       m_int2.push_back(int2[j]);
//     }
//   }
//   return DataFrame::create(Named("rt1") = m_rt1,
//                            Named("int1") = m_int1,
//                            Named("rt2") = m_rt2,
//                            Named("int2") = m_int2);
// }



// Helper function: Get maximum row
// [[Rcpp::export]]
NumericVector get_max_fun_cpp(DataFrame chromPeaks) {
  NumericVector into = chromPeaks["into"];
  int max_idx = which_max(into);
  NumericVector rtmin = chromPeaks["rtmin"];
  NumericVector rtmax = chromPeaks["rtmax"];
  NumericVector sample = chromPeaks["sample"];

  return NumericVector::create(rtmin[max_idx], rtmax[max_idx], sample[max_idx]);
}

// Main processing function
// [[Rcpp::export]]
List feature_EIC_process_fun_cpp(List input_list, int nrows, int ncols) {
  List result(nrows);

  for (int i = 0; i < nrows; ++i) {
    List chrom_peaks_list(ncols);

    for (int j = 0; j < ncols; ++j) {
      S4 xChrom = as<S4>(input_list(i, j)); // get one XChromatogram
      NumericMatrix chrom_peaks = as<NumericMatrix>(xChrom.slot("chromPeaks"));
      chrom_peaks_list[j] = chrom_peaks;
    }

    // Combine all chrom_peaks
    DataFrame chrom_peaks_df = Rcpp::as<DataFrame>(
      Rcpp::wrap(Rcpp::Language("do.call", "rbind", chrom_peaks_list).eval())
    );

    // Extract RT window from maximum peak
    NumericVector FT_max = get_max_fun_cpp(chrom_peaks_df);
    double rt_min = FT_max[0];
    double rt_max = FT_max[1];
    int sample = static_cast<int>(FT_max[2]) - 1; // adjust for 0-based index

    // Get RT and intensity from sample
    S4 xChrom_sample = as<S4>(input_list(i, sample));
    NumericVector rt_raw = xChrom_sample.slot("rtime");
    NumericVector intensity_raw = xChrom_sample.slot("intensity");

    LogicalVector valid_rt = (rt_raw >= rt_min) & (rt_raw <= rt_max);
    NumericVector rt = rt_raw[valid_rt];
    NumericVector intensity = intensity_raw[valid_rt];

    result[i] = List::create(
      Named("chrom_peaks_df") = chrom_peaks_df,
      Named("EIC_df") = DataFrame::create(
        Named("rt") = rt,
        Named("intensity") = intensity
      )
    );
  }

  return result;
}



