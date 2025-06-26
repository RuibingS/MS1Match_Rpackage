#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <omp.h>  // parallel computing

//prototypes
List feature_EIC_fun(List input_list, int nrows, int ncols);
DataFrame closest_rt_fun(NumericVector rt1, NumericVector int1,
                             NumericVector rt2, NumericVector int2, double rt_EIC_thr);
NumericVector get_max_fun_cpp(DataFrame chromPeaks);


double calCor(NumericVector col1, NumericVector col2, std::string method);





//' Function to process the feature EIC
//' @param input_list A list is the output from 'runXCMS'.
//' @param nrows A numeric value. the number of rows in the FT
//' @param ncols  A numeric value. The number of samples
//'  for each feature, the Chromatogram with the maximum 'into' will remain.
//'  @return The list with two elements: chrom_peaks_df (data frame of xChromatogram.slot"chromPeaks" XCMS and EIC_df (dataframe with two column 'rt' and 'intensity').
// [[Rcpp::export]]
List feature_EIC_process_fun(List input_list, int nrows, int ncols) {
  List result(nrows);

  for (int i = 0; i < nrows; ++i) {
    List chrom_peaks_list = List(ncols);

    NumericVector rt_combined;
    NumericVector intensity_combined;



    for (int j = 0; j < ncols; ++j) {

      S4 xChromatogram = as<S4>(input_list(i, j));

      NumericMatrix chrom_peaks = as<NumericMatrix>(xChromatogram.slot("chromPeaks"));
      chrom_peaks_list[j] = chrom_peaks;
    }


    DataFrame chrom_peaks_df = Rcpp::as<DataFrame>(
      Rcpp::wrap(Rcpp::Language("do.call", "rbind", chrom_peaks_list).eval())
    );


    int FT_index_adj = i;
    S4 FT_obj = as<S4>(input_list(FT_index_adj, 0));
    NumericVector rt, intensity;

    NumericVector FT_max = get_max_fun_cpp(chrom_peaks_df);
    double rt_min = FT_max[0];
    double rt_max = FT_max[1];
    int sample = static_cast<int>(FT_max[2]) - 1; // 0-based index


    S4 FT_sample = as<S4>(input_list(FT_index_adj, sample));
    NumericVector rt_raw = FT_sample.slot("rtime");
    NumericVector intensity_raw = FT_sample.slot("intensity");


    LogicalVector valid_ft = (rt_raw >= rt_min) & (rt_raw <= rt_max);
    rt = rt_raw[valid_ft];
    intensity = intensity_raw[valid_ft];


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






//' Function to reorder the list
//'  @param input_list A list is the output from 'runXCMS'.
//'  @param order_indices  A numeric value
//'  @param ncols  A numeric value. The number of samples
//' for each feature, the Chromatogram with the maximum 'into' will remain.
//' @return The list with two elements: chrom_peaks_df (data frame of xChromatogram.slot"chromPeaks" XCMS and EIC_df (dataframe with two column 'rt' and 'intensity').
// [[Rcpp::export]]
List reorder_input_list(List input_list, IntegerVector order_indices, int ncols) {
  int nrows = order_indices.size();
  List reordered_list(nrows * ncols);

  for (int i = 0; i < nrows; ++i) {
    int new_index = order_indices[i] - 1; // Convert to 0-based indexing

    for (int j = 0; j < ncols; ++j) {
      reordered_list(i * ncols + j) = input_list(new_index * ncols + j);
    }
  }

  return reordered_list;
}



//' Helper function: Extract feature EICs
//' @param input_list A list is the output from 'runXCMS'.
//' @param nrows  A numeic value. The number of rows in the feature table
//' @param ncols  A numeic value. The number of samples in the feature table
//' @return The list of the each feature chrom_peaks (data frame of xChromatogram.slot"chromPeaks" XCMS .
// [[Rcpp::export]]
List feature_EIC_fun(List input_list, int nrows, int ncols) {

  List result(nrows);

  // Loop over rows
  for (int i = 0; i < nrows; ++i) {

    List chrom_peaks_list(ncols);

    for (int j = 0; j < ncols; ++j) {

      S4 xChromatogram = as<S4>(input_list(i, j));

      if (!xChromatogram.hasSlot("chromPeaks")) {
        Rcpp::stop("XChromatogram at row " + std::to_string(i + 1) +
          ", column " + std::to_string(j + 1) +
          " does not have a 'chromPeaks' slot.");
      }

      NumericMatrix chrom_peaks = as<NumericMatrix>(xChromatogram.slot("chromPeaks"));

      chrom_peaks_list[j] = chrom_peaks;
    }
    result[i] = Rcpp::as<DataFrame>(Rcpp::wrap(Rcpp::Rcpp_eval(
      Rcpp::Language("do.call", "rbind", chrom_peaks_list)
    )));
  }

  return result;
}



// Helper function: Extract feature all EICs from all samples

// [[Rcpp::export]]
List feature_all_EIC_process(List input_list, DataFrame FT, int ncols) {
  int nrows = input_list.size() / ncols; // number of features
  List result(nrows);

  for (int i = 0; i < nrows; ++i) {
    List chrom_peaks_list(ncols);
    std::vector<double> all_rt, all_intensity;
    std::vector<int> all_sample_idx;

    for (int j = 0; j < ncols; ++j) {
      S4 xChromatogram = as<S4>(input_list(i, j));

      if (!xChromatogram.hasSlot("chromPeaks")) {
        Rcpp::stop("XChromatogram at row " + std::to_string(i + 1) +
          ", column " + std::to_string(j + 1) +
          " does not have a 'chromPeaks' slot.");
      }

      DataFrame chrom_peaks_df = as<DataFrame>(xChromatogram.slot("chromPeaks"));
      chrom_peaks_list[j] = chrom_peaks_df;

      if (Rf_isNull(chrom_peaks_df) || chrom_peaks_df.nrows() == 0)
        continue;

      NumericVector rtmins = chrom_peaks_df["rtmin"];
      NumericVector rtmaxs = chrom_peaks_df["rtmax"];

      NumericVector rt = xChromatogram.slot("rtime");
      NumericVector intensity = xChromatogram.slot("intensity");

      for (int k = 0; k < rtmins.size(); ++k) {
        double rtmin = rtmins[k];
        double rtmax = rtmaxs[k];

        LogicalVector valid = (rt >= rtmin) & (rt <= rtmax);
        for (int m = 0; m < rt.size(); ++m) {
          if (valid[m]) {
            all_rt.push_back(rt[m]);
            all_intensity.push_back(intensity[m]);
            all_sample_idx.push_back(j + 1);
          }
        }
      }
    }

    DataFrame EIC_all_df;

    if (all_rt.size() > 0 &&
        all_rt.size() == all_intensity.size() &&
        all_rt.size() == all_sample_idx.size()) {
      EIC_all_df = DataFrame::create(
        Named("rt") = all_rt,
        Named("intensity") = all_intensity,
        Named("sample") = all_sample_idx
      );
    } else {
      // return NA row
      EIC_all_df = DataFrame::create(
        Named("rt") = NumericVector::create(NA_REAL),
        Named("intensity") = NumericVector::create(NA_REAL),
        Named("sample") = IntegerVector::create(NA_INTEGER)
      );
    }

    result[i] = List::create(
      Named("chrom_peaks_df") = chrom_peaks_list,
      Named("EIC_all_df") = EIC_all_df
    );
  }

  return result;
}


// Helper function: Extract feature EICs which has maximum into
// [[Rcpp::export]]
List feature_EIC_process(List input_list, DataFrame FT,
                         int ncols) {
  int nrows = input_list.size() / ncols; // number of rows from the list dimensions
  List result(nrows);

  // Loop over rows
  for (int i = 0; i < nrows; ++i) {

    std::string fea_name = FT["feature_name"];
    List chrom_peaks_list = List(ncols);
    NumericVector rt, intensity; //  store RT and intensity for this feature

    for (int j = 0; j < ncols; ++j) {

      if (i < 0 || i >= nrows || j < 0 || j >= ncols) {
        Rcpp::stop("Index out of bounds: [row index = " + std::to_string(i) +
          "; column index = " + std::to_string(j) + "]");
      }

      S4 xChromatogram = as<S4>(input_list(i, j));

      if (!xChromatogram.hasSlot("chromPeaks")) {
        Rcpp::stop("XChromatogram at row " + std::to_string(i + 1) +
          ", column " + std::to_string(j + 1) +
          " does not have a 'chromPeaks' slot.");
      }

      NumericMatrix chrom_peaks = as<NumericMatrix>(xChromatogram.slot("chromPeaks"));
      chrom_peaks_list[j] = chrom_peaks;
    }

    int FT_index_adj = i;
    S4 FT_obj = as<S4>(input_list(FT_index_adj, 0)); // XChromatograms
    DataFrame FT_chromPeaks = as<DataFrame>(chrom_peaks_list[0]);

    NumericVector FT_max = get_max_fun_cpp(FT_chromPeaks);
    double rt_min = FT_max[0];
    double rt_max = FT_max[1];
    int sample = static_cast<int>(FT_max[2]) -1 ; // 0-based index

    if (sample < 0 || sample >= ncols) {
      Rcpp::stop("Sample index out of bounds: sample = " + std::to_string(sample));
    }

    S4 FT_sample = as<S4>(input_list(FT_index_adj, sample)); //specific XChromatogram
    rt = FT_sample.slot("rtime");
    intensity = FT_sample.slot("intensity");

    LogicalVector valid_ft = (rt >= rt_min) & (rt <= rt_max);
    rt = rt[valid_ft];
    intensity = intensity[valid_ft];

    result[i] = List::create(
      Named("chrom_peaks_df") = DataFrame::create(Rcpp::Named("chrom_peaks", chrom_peaks_list)),
      Named("EIC_df") = DataFrame::create(
        Named("rt") = rt,
        Named("intensity") = intensity
      )
    );
  }

  return result;
}







DataFrame pairwise_EIC_similarity(
    DataFrame FT,
    List feature_EIC_combined_result,
    double rt_threshold,
    double eic_threshold,
    double rt_EIC_thr
) {
  int n_features = FT.nrow();
  NumericVector rt = FT["rtmed"];
  std::vector<int> feature_1, feature_2;
  std::vector<double> similarity_scores;

  for (int i = 0; i < n_features; ++i) {
    for (int j = i + 1; j < n_features; ++j) {
      if (std::abs(rt[i] - rt[j]) > rt_threshold) {
        Rcpp::Rcout << "Skipping: " << i + 1 << " - " << j + 1 << " (RT threshold)\n";
        continue;
      }

      List feature1_EIC = feature_EIC_combined_result[i];
      List feature2_EIC = feature_EIC_combined_result[j];

      Rcpp::Rcout << "Checking feature pair: " << i + 1 << " - " << j + 1 << "\n";

      if (std::abs(rt[i] - rt[j]) > rt_threshold) {
        Rcpp::Rcout << "Skipping (RT threshold exceeded): " << i + 1 << " - " << j + 1 << "\n";
        continue;
      }


      if (feature1_EIC.size() == 0 || feature2_EIC.size() == 0) {
        Rcpp::Rcout << "Skipping: " << i + 1 << " - " << j + 1 << " (Empty EIC)\n";
        continue;
      }

      DataFrame EIC1 = as<DataFrame>(feature1_EIC["EIC_df"]);
      DataFrame EIC2 = as<DataFrame>(feature2_EIC["EIC_df"]);

      NumericVector rt1 = EIC1["rt"];
      NumericVector int1 = EIC1["intensity"];
      NumericVector rt2 = EIC2["rt"];
      NumericVector int2 = EIC2["intensity"];

      DataFrame merged = closest_rt_fun(rt1, int1, rt2, int2, rt_EIC_thr);
      NumericVector merged_int1 = merged["int1"];
      NumericVector merged_int2 = merged["int2"];

      LogicalVector valid = (!is_na(merged_int1)) & (!is_na(merged_int2));
      NumericVector filtered_int1 = merged_int1[valid];
      NumericVector filtered_int2 = merged_int2[valid];

      if (filtered_int1.size() < 4) {
        Rcpp::Rcout << "Skipping: " << i + 1 << " - " << j + 1 << " (Not enough points)\n";
        continue;
      }

      NumericVector log_int1 = log10(filtered_int1);
      NumericVector log_int2 = log10(filtered_int2);
      // arma::vec v1 = as<arma::vec>(log_int1);
      // arma::vec v2 = as<arma::vec>(log_int2);
      double cor_val = calCor(log_int1, log_int2,"pearson"); //arma::as_scalar(arma::cor(v1, v2));

      if (cor_val >= eic_threshold) {
        feature_1.push_back(i + 1);
        feature_2.push_back(j + 1);
        similarity_scores.push_back(cor_val);
      }
    }
  }

  return DataFrame::create(
    Named("feature_1") = feature_1,
    Named("feature_2") = feature_2,
    Named("similarity_score") = similarity_scores
  );
}





// main function: EIC matrix

// [[Rcpp::export]]
List EIC_cor_mat_fun(
    DataFrame FT,
    List feature_EIC_combined_result,
    double rt_threshold = 3,
    double rt_EIC_thr = 1
) {
  int n_features = FT.nrow();
  CharacterVector feature_names = FT["feature_name"];
  NumericVector rt = FT["rtmed"];

  std::unordered_map<std::string, int> feature_index_map;
  for (int i = 0; i < n_features; ++i) {
    feature_index_map[std::string(feature_names[i])] = i;
  }

  std::vector<std::string> feature_1, feature_2;
  std::vector<double> similarity_scores;

  for (int i = 0; i < n_features; ++i) {
    for (int j = i + 1; j < n_features; ++j) {
      if (std::abs(rt[i] - rt[j]) <= rt_threshold) {
        feature_1.push_back(std::string(feature_names[i]));
        feature_2.push_back(std::string(feature_names[j]));
      }
    }
  }

  // Rcpp::Rcout << "Total valid feature pairs: " << feature_1.size() << std::endl;

  for (size_t k = 0; k < feature_1.size(); ++k) {
    std::string name1 = feature_1[k];
    std::string name2 = feature_2[k];

    if (feature_index_map.find(name1) == feature_index_map.end() ||
        feature_index_map.find(name2) == feature_index_map.end()) {
      // Rcpp::Rcout << "Skipping: Feature not found in index map (" << name1 << ", " << name2 << ")\n";
      continue;
    }

    int i = feature_index_map[name1];
    int j = feature_index_map[name2];

    if (i >= feature_EIC_combined_result.size() || j >= feature_EIC_combined_result.size()) {
     // Rcpp::Rcout << "Skipping: Index out of bounds (" << i << ", " << j << ")\n";
      continue;
    }

    List feature1_EIC = feature_EIC_combined_result[i];
    List feature2_EIC = feature_EIC_combined_result[j];

    if (!feature1_EIC.containsElementNamed("EIC_df") ||
        !feature2_EIC.containsElementNamed("EIC_df")) {
       // Rcpp::Rcout << "Skipping: Missing EIC_df for (" << name1 << ", " << name2 << ")\n";
      continue;
    }

    DataFrame EIC1 = as<DataFrame>(feature1_EIC["EIC_df"]);
    DataFrame EIC2 = as<DataFrame>(feature2_EIC["EIC_df"]);

    NumericVector rt1 = EIC1["rt"];
    NumericVector int1 = EIC1["intensity"];
    NumericVector rt2 = EIC2["rt"];
    NumericVector int2 = EIC2["intensity"];

    DataFrame merged = closest_rt_fun(rt1, int1, rt2, int2, rt_EIC_thr);
    NumericVector merged_int1 = merged["int1"];
    NumericVector merged_int2 = merged["int2"];

    LogicalVector valid = (!is_na(merged_int1)) & (!is_na(merged_int2));
    NumericVector filtered_int1 = merged_int1[valid];
    NumericVector filtered_int2 = merged_int2[valid];

    if (filtered_int1.size() < 4) {
      // Rcpp::Rcout << "Skipping: Not enough matched RT points for (" << name1 << ", " << name2 << ")\n";
      similarity_scores.push_back(NA_REAL);
      continue;
    }

    NumericVector log_int1 = log10(filtered_int1);
    NumericVector log_int2 = log10(filtered_int2);

    double cor_val = calCor(log_int1, log_int2,"pearson");

    similarity_scores.push_back(cor_val);
  }

  if (feature_1.size() == 0) {
    return List::create(
      Named("feature_1") = CharacterVector(1, NA_STRING),
      Named("feature_2") = CharacterVector(1, NA_STRING),
      Named("similarity_score") = NumericVector(1, NA_REAL)
    );
  }

  while (similarity_scores.size() < feature_1.size()) {
    similarity_scores.push_back(NA_REAL);
  }

  return List::create(
    Named("feature_1") = feature_1,
    Named("feature_2") = feature_2,
    Named("similarity_score") = similarity_scores
  );
}


// [[Rcpp::export]]
List EIC_cor_fun(
    DataFrame FT,
    List feature_EIC_combined_result,
    double rt_threshold = 3,
    double eic_threshold = 0.8,
    double rt_EIC_thr = 1
) {
  int n_features = FT.nrow();
  CharacterVector feature_names = FT["feature_name"];
  NumericVector rt = FT["rtmed"];

  // Step 1: map from feature_name to its index
  std::unordered_map<std::string, int> feature_index_map;
  for (int i = 0; i < n_features; ++i) {
    feature_index_map[std::string(feature_names[i])] = i;
  }

  std::vector<std::string> feature_1, feature_2;
  std::vector<double> similarity_scores;

  // Step 2: feature pairs
  for (int i = 0; i < n_features; ++i) {
    for (int j = i + 1; j < n_features; ++j) {
      if (std::abs(rt[i] - rt[j]) <= rt_threshold) {
        feature_1.push_back(std::string(feature_names[i]));
        feature_2.push_back(std::string(feature_names[j]));
      }
    }
  }

  Rcpp::Rcout << "Total valid feature pairs: " << feature_1.size() << std::endl;

  // Step 3: Compute Similarity for Each Pair
  for (size_t k = 0; k < feature_1.size(); ++k) {
    std::string name1 = feature_1[k];
    std::string name2 = feature_2[k];

    // Lookup indices
    if (feature_index_map.find(name1) == feature_index_map.end() ||
        feature_index_map.find(name2) == feature_index_map.end()) {
      Rcpp::Rcout << "Skipping: Feature not found in index map (" << name1 << ", " << name2 << ")\n";
      continue;
    }

    int i = feature_index_map[name1];
    int j = feature_index_map[name2];

    if (i >= feature_EIC_combined_result.size() || j >= feature_EIC_combined_result.size()) {
      Rcpp::Rcout << "Skipping: Index out of bounds (" << i << ", " << j << ")\n";
      continue;
    }

    List feature1_EIC = feature_EIC_combined_result[i];
    List feature2_EIC = feature_EIC_combined_result[j];

    if (!feature1_EIC.containsElementNamed("EIC_df") ||
        !feature2_EIC.containsElementNamed("EIC_df")) {
        Rcpp::Rcout << "Skipping: Missing EIC_df for (" << name1 << ", " << name2 << ")\n";
      continue;
    }

    DataFrame EIC1 = as<DataFrame>(feature1_EIC["EIC_df"]);
    DataFrame EIC2 = as<DataFrame>(feature2_EIC["EIC_df"]);

    NumericVector rt1 = EIC1["rt"];
    NumericVector int1 = EIC1["intensity"];
    NumericVector rt2 = EIC2["rt"];
    NumericVector int2 = EIC2["intensity"];

    // Step 4: Match Closest RT
    DataFrame merged = closest_rt_fun(rt1, int1, rt2, int2, rt_EIC_thr);
    NumericVector merged_int1 = merged["int1"];
    NumericVector merged_int2 = merged["int2"];

    LogicalVector valid = (!is_na(merged_int1)) & (!is_na(merged_int2));
    NumericVector filtered_int1 = merged_int1[valid];
    NumericVector filtered_int2 = merged_int2[valid];

    if (filtered_int1.size() < 4) {
      Rcpp::Rcout << "Skipping: Not enough matched RT points for (" << name1 << ", " << name2 << ")\n";
      continue;
    }

    // Step 5: Log-Transformed Correlation
    NumericVector log_int1 = log10(filtered_int1);
    NumericVector log_int2 = log10(filtered_int2);



    double cor_val = calCor(log_int1,log_int2,  "pearson");
    if (cor_val >= eic_threshold) {
      similarity_scores.push_back(cor_val);}
    else {
      similarity_scores.push_back(NA_REAL);
      }
  }

  if (feature_1.size() == 0) {
    return List::create(
      Named("feature_1") = CharacterVector(1, NA_STRING),
      Named("feature_2") = CharacterVector(1, NA_STRING),
      Named("similarity_score") = NumericVector(1, NA_REAL)
    );
  }

  while (similarity_scores.size() < feature_1.size()) {
    similarity_scores.push_back(NA_REAL);
  }

  return List::create(
    Named("feature_1") = feature_1,
    Named("feature_2") = feature_2,
    Named("similarity_score") = similarity_scores
  );
}





// Helper function: get the maximum row
NumericVector get_max_fun_cpp(DataFrame chromPeaks) {
  NumericVector into = chromPeaks["into"];
  int max_idx = which_max(into);
  NumericVector rtmin = chromPeaks["rtmin"];
  NumericVector rtmax = chromPeaks["rtmax"];
  NumericVector sample = chromPeaks["sample"];

  return NumericVector::create(rtmin[max_idx], rtmax[max_idx], sample[max_idx]);
}


// Helper function: Function to get the closed rt


// [[Rcpp::export]]
DataFrame closest_rt_fun(NumericVector rt1, NumericVector int1,
                             NumericVector rt2, NumericVector int2, double rt_EIC_thr) {
  int n1 = rt1.size();
  int n2 = rt2.size();

  std::vector<double> matched_rt1, matched_int1, matched_rt2, matched_int2;

  int j = 0;
  for (int i = 0; i < n1; ++i) {
    while (j < n2 && rt2[j] < rt1[i] - rt_EIC_thr) {
      j++;
    }
    if (j < n2 && std::abs(rt1[i] - rt2[j]) <= rt_EIC_thr) {
      matched_rt1.push_back(rt1[i]);
      matched_int1.push_back(int1[i]);
      matched_rt2.push_back(rt2[j]);
      matched_int2.push_back(int2[j]);
    }
  }

  return DataFrame::create(
    Named("rt1") = matched_rt1,
    Named("int1") = matched_int1,
    Named("rt2") = matched_rt2,
    Named("int2") = matched_int2
  );
}




// [[Rcpp::export]]
double calCor(NumericVector col1, NumericVector col2, std::string method) {
  LogicalVector completeIndices = (!is_na(col1)) & (!is_na(col2));

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
NumericMatrix pairwiseCor(NumericMatrix x, std::string method) {
  int nCols = x.ncol();
  NumericMatrix result(nCols, nCols);

  for (int i = 0; i < nCols - 1; ++i) {
    for (int j = i + 1; j < nCols; ++j) {
      NumericVector col1 = x(_, i);
      NumericVector col2 = x(_, j);
      result(i, j) = calCor(col1, col2, method);
      result(j, i) = result(i, j);
    }
  }

  for (int i = 0; i < nCols; ++i) {
    result(i, i) = 1.0;
  }

  return result;
}








