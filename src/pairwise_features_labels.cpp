#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
using namespace Rcpp;

// Helper: collapse vector into string
std::string collapse_cpp(const std::vector<std::string>& vec, std::string sep = ";") {
  std::string result;
  for (size_t i = 0; i < vec.size(); ++i) {
    result += vec[i];
    if (i != vec.size() - 1) result += sep;
  }
  return result;
}

// [[Rcpp::export]]
DataFrame pairwise_features_labels(DataFrame feature_table,
                                               List feature_compounds,
                                               double rt_threshold = 5.0) {

  // Extract feature names and retention times
  CharacterVector feature_names = feature_table["feature_name"];
  NumericVector rtmed = feature_table["rtmed"];
  int n = feature_names.size();

  // Build feature-to-compound map
  std::unordered_map<std::string, std::vector<std::string>> feature2compounds;

  for (int i = 0; i < feature_compounds.size(); ++i) {
    List compound_entry = feature_compounds[i];
    CharacterVector peaks = compound_entry["Matched_features_to_peaks"];
    std::string cname_str = Rcpp::as<std::string>(compound_entry["Compound_name"]);

    for (auto& peak : peaks) {
      feature2compounds[Rcpp::as<std::string>(peak)].push_back(cname_str);
    }
  }

  std::vector<std::string> feature_1_vec, feature_2_vec, shared_compounds_vec;
  std::vector<int> label_vec;

  // RT filtering
  for (int i = 0; i < n; ++i) {
    std::string f1 = Rcpp::as<std::string>(feature_names[i]);
    double rt1 = rtmed[i];

    // Skip features with no compound assignments
    if (feature2compounds.find(f1) == feature2compounds.end()) continue;

    for (int j = i + 1; j < n; ++j) {
      double rt2 = rtmed[j];

      if (rt2 - rt1 > rt_threshold) break;

      std::string f2 = Rcpp::as<std::string>(feature_names[j]);

      if (feature2compounds.find(f2) == feature2compounds.end()) continue;

      // Get compounds for features
      const std::vector<std::string>& compounds1 = feature2compounds[f1];
      const std::vector<std::string>& compounds2 = feature2compounds[f2];

      // Find shared compounds
      std::unordered_set<std::string> set1(compounds1.begin(), compounds1.end());
      std::vector<std::string> shared;
      for (const auto& c : compounds2) {
        if (set1.find(c) != set1.end()) {
          shared.push_back(c);
        }
      }

      // Label = 1 if shared compounds, else 0
      int label = shared.empty() ? 0 : 1;
      std::string shared_str = collapse_cpp(shared);

      feature_1_vec.push_back(f1);
      feature_2_vec.push_back(f2);
      shared_compounds_vec.push_back(shared_str);
      label_vec.push_back(label);
    }
  }

  return DataFrame::create(
    _["feature_1"] = feature_1_vec,
    _["feature_2"] = feature_2_vec,
    _["shared_compounds"] = shared_compounds_vec,
    _["label"] = label_vec
  );
}
