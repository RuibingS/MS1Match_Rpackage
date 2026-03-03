// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>


using namespace Rcpp;


// [[Rcpp::export]]
DataFrame CorPairs_cpp(NumericMatrix log_intensity_t,
                       CharacterVector feature_names,
                       DataFrame feature_pairs,
                       std::string method = "pearson") {
  arma::mat intensity_mat = as<arma::mat>(log_intensity_t); // samples x features (columns are features)
  std::unordered_map<std::string, int> name_to_col;

  for (int i = 0; i < feature_names.size(); ++i)
    name_to_col[as<std::string>(feature_names[i])] = i;

  CharacterVector f1 = feature_pairs["feature_1"];
  CharacterVector f2 = feature_pairs["feature_2"];
  int n = f1.size();

  std::vector<double> cor_vals(n);

  for (int i = 0; i < n; ++i) {
    auto it1 = name_to_col.find(as<std::string>(f1[i]));
    auto it2 = name_to_col.find(as<std::string>(f2[i]));

    if (it1 == name_to_col.end() || it2 == name_to_col.end()) {
      cor_vals[i] = NA_REAL;
      continue;
    }

    arma::vec x = intensity_mat.col(it1->second);
    arma::vec y = intensity_mat.col(it2->second);

    arma::uvec valid = arma::find_finite(x % y);
    if (valid.n_elem < 3) {
      cor_vals[i] = NA_REAL;
      continue;
    }

    arma::vec x_filt = x.elem(valid);
    arma::vec y_filt = y.elem(valid);

    if (method == "pearson") {
      cor_vals[i] = arma::as_scalar(arma::cor(x_filt, y_filt));
    } else if (method == "spearman") {
      arma::vec x_rank = arma::conv_to<arma::vec>::from(arma::sort_index(arma::sort_index(x_filt))) + 1;
      arma::vec y_rank = arma::conv_to<arma::vec>::from(arma::sort_index(arma::sort_index(y_filt))) + 1;
      cor_vals[i] = arma::as_scalar(arma::cor(x_rank, y_rank));
    } else {
      Function cor("cor");
      cor_vals[i] = as<double>(cor(wrap(x_filt), wrap(y_filt), Named("method") = method));
    }
    // negative correlations 0
    if (!NumericVector::is_na(cor_vals[i]) && cor_vals[i] < 0) {
      cor_vals[i] = 0.0;
      }

  }

  return DataFrame::create(
    _["feature_1"] = f1,
    _["feature_2"] = f2,
    _["cor_value"] = cor_vals
  );
}

