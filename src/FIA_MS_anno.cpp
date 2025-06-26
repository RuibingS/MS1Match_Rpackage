// This is the annotation functions for FIA-MS data.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iomanip>
#include <regex>
#include "utils_annotation.h"

using namespace Rcpp;
//prototypes
//
// std::string join_annotations(const std::vector<std::string>& annotations, const std::string& delimiter);

Rcpp::List PImatch_fia_fun(Rcpp::DataFrame FT,
                       Rcpp::DataFrame Comp_data,
                       Rcpp::CharacterVector PIon,
                       double diff_mz_thr,
                       int ppm);

Rcpp::List findNL_fia_fun(
    Rcpp::DataFrame NL_data,
    Rcpp::DataFrame adj_long,
    double diff_mz_thr,
    int ppm);



CharacterVector Check_Iso_Charge_fia(CharacterVector name,
                                 NumericVector mz,
                                 NumericVector intensity,
                                 double iso_threshold ,
                                 int max_charge ,
                                 int max_isotope);



// double compute_adduct_mass(double M, std::string adduct_formula);


// std::string join_annotations(const std::vector<std::string>& annotations, const std::string& delimiter = "; ") {
//   std::string result;
//   for (size_t i = 0; i < annotations.size(); ++i) {
//     result += annotations[i];
//     if (i < annotations.size() - 1) result += delimiter;
//   }
//   return result;
// }

// double compute_adduct_mass(double M, std::string adduct_formula) {
//   adduct_formula.erase(std::remove(adduct_formula.begin(), adduct_formula.end(), ' '), adduct_formula.end());
//   std::regex expr_regex("([+-]?[0-9]*\\.?[0-9]*)M(/([0-9]+))?");
//   std::smatch match;
//   double computed_mass = 0.0;
//   std::string formula = adduct_formula;
//   while (std::regex_search(formula, match, expr_regex)) {
//     double coefficient = match[1].str().empty() ? 1.0 : std::stod(match[1].str());
//     double divisor = match[3].str().empty() ? 1.0 : std::stod(match[3].str());
//     computed_mass += (coefficient * M) / divisor;
//     formula = match.suffix();
//   }
//   if (!formula.empty()) {
//     try {
//       computed_mass += std::stod(formula);
//     } catch (...) {
//       Rcpp::Rcerr << "Warning: Could not parse: " << formula << std::endl;
//     }
//   }
//   return computed_mass;
// }
//



// [[Rcpp::export]]
List annotate_peaks_fia(List compounds_list, DataFrame FT,
                        DataFrame adducts_df, DataFrame neutral_loss_table,
                        std::string polarity, NumericVector intensity,
                        double threshold = 0.005,
                        double iso_threshold = 1.003355, int max_charge = 3, int max_isotope = 4) {

  CharacterVector feature_name = FT["feature_name"];
  NumericVector FT_mz = FT["mz"];

  int n_features = FT.nrows();
  int n_compounds = compounds_list.size();

  NumericVector nl_mz_difference = neutral_loss_table["Accurate.Mass"];
  CharacterVector nl_name = neutral_loss_table["Formula"];
  CharacterVector pos_polarity = neutral_loss_table["Pos"];
  CharacterVector neg_polarity = neutral_loss_table["Neg"];

  CharacterVector adduct_ions = adducts_df["Ion.name"];
  CharacterVector adduct_polarities = adducts_df["Ion.polarity"];
  CharacterVector adduct_formulas = adducts_df["Ion.mass"];

  List output_list(n_features);
  std::vector<std::string> PI_annotations(n_features, "");
  std::vector<std::string> adduct_annotations(n_features, "");
  std::vector<std::string> NL_annotations(n_features, "");

  for (int f = 0; f < n_features; f++) {
    double feature_mz = FT_mz[f];
    std::vector<std::string> PI_matches;

    for (int i = 0; i < n_compounds; i++) {
      List compound = compounds_list[i];
      if (!compound.containsElementNamed("ExactMass")) continue;
      double M = compound["ExactMass"];

      std::vector<std::pair<double, std::string>> PI_ions;
      if (polarity == "pos") {
        PI_ions = {{M + 1.007276, "[M+H]+"}, {M + 22.989218, "[M+Na]+"}};
      } else {
        PI_ions = {{M - 1.007276, "[M-H]-"}};
      }

      for (const auto& ion : PI_ions) {
        if (std::abs(feature_mz - ion.first) <= threshold) {
          PI_matches.push_back(as<std::string>(compound["Name"]) + " " + ion.second);
        }
      }
    }
    PI_annotations[f] = join_annotations(PI_matches);
  }

  // Adduct match
  // for (int i = 0; i < n_compounds; ++i) {
  //   List compound = compounds_list[i];
  //   if (!compound.containsElementNamed("ExactMass")) continue;
  //   double M = compound["ExactMass"];
  //
  //   for (int k = 0; k < adduct_formulas.size(); ++k) {
  //     if (as<std::string>(adduct_polarities[k]) != polarity) continue;
  //     double mass = compute_adduct_mass(M, as<std::string>(adduct_formulas[k]));
  //
  //     for (int g = 0; g < n_features; ++g) {
  //       if (std::abs(FT_mz[g] - mass) <= threshold) {
  //         std::string annotation = as<std::string>(adduct_ions[k]) + " " + std::to_string(M);
  //         if (adduct_annotations[g].empty()) adduct_annotations[g] = annotation;
  //         else adduct_annotations[g] += "; " + annotation;
  //       }
  //     }
  //   }
  // }

  // Neutral Loss
  for (int f = 0; f < n_features; f++) {
    double mz_f = FT_mz[f];
    std::vector<std::string> matches;

    for (int g = 0; g < n_features; g++) {
      if (g == f) continue;
      double mz_g = FT_mz[g];
      if (mz_g <= mz_f) continue;

      double diff = mz_g - mz_f;
      for (int m = 0; m < nl_mz_difference.size(); m++) {
        bool match_polarity = ((polarity == "pos" && as<std::string>(pos_polarity[m]) == "+") ||
                               (polarity == "neg" && as<std::string>(neg_polarity[m]) == "+"));

        if (match_polarity && std::abs(diff - nl_mz_difference[m]) <= threshold) {
          matches.push_back("-" + as<std::string>(nl_name[m]) + " (" + as<std::string>(feature_name[g]) + ")");
        }
      }
    }
    NL_annotations[f] = join_annotations(matches);
  }

  CharacterVector isotope_labels = Check_Iso_Charge_fia(feature_name, FT_mz, intensity, iso_threshold, max_charge, max_isotope);

  for (int f = 0; f < n_features; f++) {
    output_list[f] = List::create(
      Named("feature_name") = as<std::string>(feature_name[f]),
      Named("mz") = FT_mz[f],
      Named("PI_match") = PI_annotations[f],
      // Named("adduct_match") = adduct_annotations[f],
      Named("NL_match") = NL_annotations[f],
      Named("iso_annotation") = as<std::string>(isotope_labels[f])
    );
  }

  return output_list;
}

// [[Rcpp::export]]
List PImatch_fia_fun(DataFrame FT, DataFrame Comp_data, CharacterVector PIon, double diff_mz_thr, int ppm) {
  int nrow_FT = FT.nrow();
  int nrow_Comp = Comp_data.nrow();
  int length_PIon = PIon.length();

  CharacterVector feature_name = FT["feature_name"];
  NumericVector FT_mz = FT["mz"];

  CharacterVector Comp_name = Comp_data["Name"];
  NumericVector Comp_M = Comp_data["exactMass"];

  List output(nrow_FT);
  for (int i = 0; i < nrow_FT; i++) {
    List df(5);
    df.names() = CharacterVector::create("Feature_name", "Feature_mz", "Comp_name", "PI_name", "M");

             CharacterVector feaNames;
             CharacterVector CompNames;
             CharacterVector PIname;
             NumericVector fea_mz;
             NumericVector M_vec;

             for (int p = 0; p < length_PIon; p++) {
               String PIon_name = PIon[p];
               NumericVector comp_mz = Comp_data[PIon_name];

               for (int j = 0; j < nrow_Comp; j++) {
                 float diff_mz = std::fabs(FT_mz[i] - comp_mz[j]);

                 if ((diff_mz <= comp_mz[j] * ppm * 1e-6) || (diff_mz <= diff_mz_thr)) {
                   feaNames.push_back(feature_name[i]);
                   fea_mz.push_back(FT_mz[i]);
                   CompNames.push_back(Comp_name[j]);
                   PIname.push_back(PIon[p]);
                   M_vec.push_back(Comp_M[j]);
                 }
               }
             }

             df["Feature_name"] = feaNames;
             df["Feature_mz"] = fea_mz;
             df["Comp_name"] = CompNames;
             df["PI_name"] = PIname;
             df["M"] = M_vec;
             output[i] = df;
  }
  return output;
}

// [[Rcpp::export]]
List findNL_fia_fun(DataFrame NL_data, DataFrame adj_long, double diff_mz_thr, int ppm) {
  int nrow_adj_long = adj_long.nrow();
  int nrow_NL = NL_data.nrow();

  NumericVector mz_diff = adj_long["mz_diff"];
  CharacterVector fea_name1 = adj_long["Var1"];
  CharacterVector fea_name2 = adj_long["Var2"];
  NumericVector mz1 = adj_long["mz_x"];
  NumericVector mz2 = adj_long["mz_y"];

  NumericVector NL_mz = NL_data["Accurate.Mass"];
  CharacterVector NL_name = NL_data["Neutral.Loss"];

  List out(nrow_adj_long);

  for (int i = 0; i < nrow_adj_long; i++) {
    List df(3);
    df.names() = CharacterVector::create("Feature_name1", "Feature_name2", "NL_Names");
    CharacterVector feaNames1;
    CharacterVector feaNames2;
    CharacterVector NLNames;

    for (int j = 0; j < nrow_NL; j++) {
      float diff_mz1 = std::fabs(NL_mz[j] - mz_diff[i]);
      if (diff_mz1 <= diff_mz_thr || diff_mz1 <= NL_mz[j] * ppm * 1e-6) {
        if (mz1[i] < mz2[i]) {
          feaNames1.push_back(fea_name2[i]);
          feaNames2.push_back(fea_name1[i]);
        } else {
          feaNames1.push_back(fea_name1[i]);
          feaNames2.push_back(fea_name2[i]);
        }
        NLNames.push_back(NL_name[j]);
      }
    }

    df["Feature_name1"] = feaNames1;
    df["Feature_name2"] = feaNames2;
    df["NL_Names"] = NLNames;
    out[i] = df;
  }

  return out;
}
// [[Rcpp::export]]
Rcpp::CharacterVector Check_Iso_Charge_fia(Rcpp::CharacterVector name,
                                       Rcpp::NumericVector mz,
                                       Rcpp::NumericVector intensity,
                                       double iso_threshold = 1.003355,
                                       int max_charge = 3,
                                       int max_isotope = 4) {
  bool is_sorted = true;
  for (int i = 1; i < mz.size(); ++i) {
    if (mz[i] < mz[i - 1]) {
      is_sorted = false;
      break;
    }
  }

  if (!is_sorted) {
    Rcpp::Rcout << "Warning: input m/z is not sorted in ascending order. Please make sure the feature table is sorted by m/z.\n";
  }

  Rcpp::CharacterVector labels(name.size());
  int groupIndex = 1;

  for (int i = 0; i < name.size(); i++) {
    if (labels[i] != "") continue;

    std::vector<int> maxIsotopeForCharge(max_charge + 1, 0);

    for (int j = i + 1; j < name.size(); j++) {
      if (labels[j] != "") continue;

      double mz_diff = mz[j] - mz[i];

      for (int c = 1; c <= max_charge; c++) {
        for (int m = 1; m <= max_isotope; m++) {
          double expected_mz_diff = m * iso_threshold / c;
          double max_int_ratio = (mz[i] * c / 12) * 0.011;

          if (m == 1) {
            if (std::abs(mz_diff - expected_mz_diff) <= 0.002 &&
                intensity[i] * max_int_ratio >= intensity[j]) {

              labels[i] = Rcpp::String("[" + std::to_string(groupIndex) + "][M]" +
                (c > 1 ? std::to_string(c) : "") + "+");

              labels[j] = Rcpp::String("[" + std::to_string(groupIndex) + "][M+" +
                std::to_string(m) + "]" +
                (c > 1 ? std::to_string(c) : "") + "+");

              maxIsotopeForCharge[c] = m;
            }
          } else {
            if (std::abs(mz_diff - expected_mz_diff) <= 0.01 &&
                intensity[j] <= intensity[i]) {
              if (maxIsotopeForCharge[c] == m - 1) {
                labels[j] = Rcpp::String("[" + std::to_string(groupIndex) + "][M+" +
                  std::to_string(m) + "]" +
                  (c > 1 ? std::to_string(c) : "") + "+");

                maxIsotopeForCharge[c] = m;
              }
            }
          }
        }
      }
    }

    if (!labels[i].empty()) {
      groupIndex++;
    }
  }

  return labels;
}
