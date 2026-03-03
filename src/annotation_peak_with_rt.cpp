// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <limits>
#include <sstream>
#include <iomanip>
#include <algorithm> // transform, min
#include "utils_annotation.h"
#include "utils_sort.h"
#include "utils_check.h"
#include "utils_isotope.h"

using namespace Rcpp;


//prototypes

Rcpp::List PImatch_fun(Rcpp::DataFrame FT,
                       Rcpp::DataFrame Comp_data,
                       Rcpp::CharacterVector PIon,
                       double diff_mz_thr,
                       double diff_rt_thr,
                       int ppm);

Rcpp::List findNL_fun(
    Rcpp::DataFrame NL_data,
    Rcpp::DataFrame adj_long,
    double diff_mz_thr,
    double diff_rt_thr,
    int ppm);



CharacterVector Check_Iso_Charge(CharacterVector name,
                                 NumericVector mz,
                                 NumericVector rt,
                                 NumericVector intensity,
                                 double iso_threshold ,
                                 int max_charge ,
                                 int max_isotope);


double compute_adduct_mass(double M, std::string adduct_formula);


std::string collapse_matches(const std::vector<std::string>& matches) {
  if (matches.empty()) return "";
  std::ostringstream oss;
  oss << matches[0];
  for (size_t i = 1; i < matches.size(); ++i) {
    oss << "; " << matches[i];
  }
  return oss.str();
}


// [[Rcpp::export]]
List annotate_with_rt_cpp(List compounds_list,
                          DataFrame FT,
                          DataFrame adducts_df,
                          DataFrame neutral_loss_table,
                          std::string polarity,

                          double mz_threshold = 0.005,
                          double ppm = 20,                // ppm for adduct matching
                          double PI_rt_threshold = 5.0,
                          double rt_threshold = 3.0,
                          double iso_threshold = 1.003355,
                          int max_charge = 3,
                          int max_isotope = 4,
                          bool run_PI_match = true,
                          bool run_adduct_match = true,
                          bool run_NL_match = true) {

  NumericVector FT_mz = get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
  CharacterVector feature_name = get_required_column<CharacterVector>(FT, "feature_name");
  NumericVector FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention time");
  NumericVector FT_intensity = get_intensity_from_FT(FT);

  // ensure FT is sorted by m/z (stable)
  bool is_sorted = true;
  for (int i = 1; i < FT_mz.size(); ++i) if (FT_mz[i] < FT_mz[i-1]) { is_sorted = false; break; }
  if (!is_sorted) {
    Rcpp::Rcout << "Sorting FT by m/z...\n";
    std::vector<int> idx = compute_sorted_index(FT_mz);
    FT_mz = reorder_vector(FT_mz, idx);
    FT_rt = reorder_vector(FT_rt, idx);
    feature_name = reorder_vector(feature_name, idx);
    FT_intensity = reorder_vector(FT_intensity, idx);
  }

  const int n_features  = FT.nrows();
  const int n_compounds = compounds_list.size();
  Rcpp::Rcout << "Features: " << n_features << ", Compounds: " << n_compounds << "\n";

  std::transform(polarity.begin(), polarity.end(), polarity.begin(), ::tolower);
  if (polarity == "positive" || polarity == "+" || polarity == "pos") polarity = "pos";
  else if (polarity == "negative" || polarity == "-" || polarity == "neg") polarity = "neg";

  std::vector<double> compound_mass(n_compounds, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> compound_rt (n_compounds, 0.0);
  std::vector<std::string> compound_name(n_compounds);

  for (int i = 0; i < n_compounds; ++i) {
    List compound = compounds_list[i];
    if (compound.containsElementNamed("ExactMass")) {
      try { compound_mass[i] = Rcpp::as<double>(compound["ExactMass"]); } catch (...) {}
    }
    if (compound.containsElementNamed("RetTime")) {
      try { compound_rt[i] = std::atof(Rcpp::as<std::string>(compound["RetTime"]).c_str()); } catch (...) {}
    }
    if (compound.containsElementNamed("Name")) {
      try { compound_name[i] = Rcpp::as<std::string>(compound["Name"]); } catch (...) {}
    }
  }

  // adduct table columns
  CharacterVector adduct_ions      = adducts_df["Ion.name"];
  CharacterVector adduct_polarity  = adducts_df["Ion.polarity"];
  CharacterVector adduct_formulas  = adducts_df["Ion.mass"];
  CharacterVector adduct_charge    = adducts_df.containsElementNamed("Charge")
    ? (CharacterVector)adducts_df["Charge"]
  : CharacterVector(adduct_ions.size(), NA_STRING);

  // Check duplicate adduct table
  struct KeyHash {
    std::size_t operator()(const std::string& s) const noexcept {
      return std::hash<std::string>{}(s);
    }
  };
  std::unordered_set<std::string, KeyHash> adduct_seen;
  std::vector<std::string> filtered_ions;
  std::vector<std::string> filtered_formulas;

  for (int i = 0; i < adduct_ions.size(); ++i) {
    std::string pol = Rcpp::as<std::string>(adduct_polarity[i]);
    std::transform(pol.begin(), pol.end(), pol.begin(), ::tolower);
    if (pol == "positive" || pol == "+" || pol == "pos") pol = "pos";
    else if (pol == "negative" || pol == "-" || pol == "neg") pol = "neg";

    if (pol != polarity) continue;

    std::string ion  = Rcpp::as<std::string>(adduct_ions[i]);
    std::string mass = Rcpp::as<std::string>(adduct_formulas[i]);
    std::string chg  = adduct_charge[i] == NA_STRING ? "" : Rcpp::as<std::string>(adduct_charge[i]);

    std::string key = ion + "||" + mass + "||" + chg + "||" + pol;
    if (adduct_seen.insert(key).second) {
      filtered_ions.push_back(ion);
      filtered_formulas.push_back(mass);
    }
  }

  // neutral loss
  NumericVector nl_mz_difference = neutral_loss_table["Accurate.Mass"];
  CharacterVector nl_name = neutral_loss_table["Formula"];
  CharacterVector pos_pol = neutral_loss_table["Pos"];
  CharacterVector neg_pol = neutral_loss_table["Neg"];

  std::vector<std::string> PI_annotations(n_features);
  std::vector<std::string> adduct_annotations(n_features);
  std::vector<std::string> NL_annotations(n_features);

  // ----------------- PI match -----------------
  if (run_PI_match) {
    Rcpp::Rcout << "Starting PI match...\n";
    int skipped_no_mass = 0, skipped_no_rt = 0, skipped_rt_tol = 0, skipped_mz_tol = 0, matched_count = 0;

    for (int f = 0; f < n_features; ++f) {
      double mz = FT_mz[f];
      double rt = FT_rt[f];

      for (int i = 0; i < n_compounds; ++i) {
        double M = compound_mass[i];
        if (std::isnan(M)) { skipped_no_mass++; continue; }
        if (compound_rt[i] == 0.0) { skipped_no_rt++; continue; }

        double theo_mz = (polarity == "pos") ? (M + 1.007276) : (M - 1.007276);
        double rt_diff = std::fabs(rt - compound_rt[i]);
        if (rt_diff > PI_rt_threshold) { skipped_rt_tol++; continue; }

        double mz_diff = std::fabs(mz - theo_mz);
        if (mz_diff > mz_threshold) { skipped_mz_tol++; continue; }

        // M to 5 decimal places:
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(5) << M;
        std::string M_str = oss.str();

        PI_annotations[f] = compound_name[i] + ((polarity == "pos") ? " [M+H]+ " : " [M-H]- ") + M_str;
        matched_count++;
        break;
      }
    }
    Rcpp::Rcout << "PI match done.\n";
  }

  // ----------------- Adduct match (with de-dup) -----------------
  std::vector<std::vector<std::string>> all_adduct_matches(n_features);

  if (run_adduct_match) {
    Rcpp::Rcout << "Starting adduct match...\n";

    for (int f = 0; f < n_features; ++f) {
      double mz_feature = FT_mz[f];

      std::unordered_set<std::string, KeyHash> seen;
      std::vector<std::string> matches;
      matches.reserve(16);

      for (int c = 0; c < n_compounds; ++c) {
        double M = compound_mass[c];
        if (std::isnan(M)) continue;

        // Compute adduct m/z vector from formulas
        const std::vector<std::string>& forms = filtered_formulas;
        std::vector<double> mz_adducts = parse_adduct_vector(M, forms);

        size_t L = std::min(mz_adducts.size(), filtered_ions.size());
        for (size_t k = 0; k < L; ++k) {
          double mz_adduct = mz_adducts[k];
          double diff = std::fabs(mz_feature - mz_adduct);
          double mz_tol = std::max(mz_threshold, mz_adduct * ppm * 1e-6);

          if (diff <= mz_tol) {
            std::ostringstream oss;
            oss << filtered_ions[k] << " " << std::fixed << std::setprecision(6) << M;
            std::string annotation = oss.str();

            if (seen.insert(annotation).second) {
              matches.push_back(annotation);
            }
          }
        }
      }

      // Save matches for this feature
      all_adduct_matches[f] = matches;
    }

    Rcpp::Rcout << "Adduct annotation completed.\n";
  }


  // ----------------- Neutral loss match -----------------
  if (run_NL_match) {
    Rcpp::Rcout << "Starting neutral loss match...\n";
    for (int f1 = 0; f1 < n_features; ++f1) {
      double mz1 = FT_mz[f1];
      double rt1 = FT_rt[f1];
      std::vector<std::string> nl_matches;

      for (int f2 = 0; f2 < n_features; ++f2) {
        if (f1 == f2) continue;
        double mz2 = FT_mz[f2];
        double rt2 = FT_rt[f2];
        if (std::fabs(rt1 - rt2) > rt_threshold || mz2 <= mz1) continue;

        double diff = mz2 - mz1;
        for (int m = 0; m < nl_mz_difference.size(); ++m) {
          bool match_ok = ((polarity == "pos" && as<std::string>(pos_pol[m]) == "+") ||
                           (polarity == "neg" && as<std::string>(neg_pol[m]) == "+"));
          if (match_ok && std::fabs(diff - nl_mz_difference[m]) <= mz_threshold) {
            nl_matches.push_back("<- " + as<std::string>(nl_name[m]) + " (" + as<std::string>(feature_name[f2]) + ")");
          }
        }
      }
      NL_annotations[f1] = join_annotations(nl_matches);
    }
  }

  // ----------------- Isotope / charge annotations -----------------
  Rcpp::Rcout << "Starting isotope annotation...\n";
  CharacterVector isotope_labels = Check_Iso_Charge(feature_name, FT_mz, FT_rt, FT_intensity,
                                                    iso_threshold, max_charge, max_isotope);
  Rcpp::Rcout << "Isotope annotation completed.\n";

  // ----------------- output -----------------
  List output_list(n_features);
  for (int i = 0; i < n_features; ++i) {
    output_list[i] = List::create(
      Named("feature_name")   =  CharacterVector::create(feature_name[i]),
      Named("mz")             = FT_mz[i],
      Named("rt")             = FT_rt[i],
      Named("PI_match")       = PI_annotations[i],
      Named("adduct_match")   = collapse_matches(all_adduct_matches[i]),
      Named("NL_match")       = NL_annotations[i],
      Named("iso_annotation") =  CharacterVector::create(isotope_labels[i])
    );
  }

  Rcpp::Rcout << "Annotation completed.\n";
  return output_list;
}









// Isotope annotation

Rcpp::CharacterVector Check_Iso_Charge(Rcpp::CharacterVector name,
                                       Rcpp::NumericVector mz,
                                       Rcpp::NumericVector rt,
                                       Rcpp::NumericVector intensity,
                                       double iso_threshold = 1.003355,
                                       int max_charge = 3,
                                       int max_isotope = 4) {
  // Just warn if mz is not sorted
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
      double rt_diff = std::abs(rt[j] - rt[i]);

      for (int c = 1; c <= max_charge; c++) {
        for (int m = 1; m <= max_isotope; m++) {
          double expected_mz_diff = m * iso_threshold / c;
          double max_int_ratio = (mz[i] * c / 12) * 0.011;

          if (m == 1) {
            if (std::abs(mz_diff - expected_mz_diff) <= 0.002 &&
                rt_diff < 3 &&
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
                rt_diff < 3 &&
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




Rcpp::List PImatch_fun(Rcpp::DataFrame FT,
                       Rcpp::DataFrame Comp_data,
                       Rcpp::CharacterVector PIon,
                       double diff_mz_thr,
                       double diff_rt_thr,
                       int ppm) {
  int nrow_FT = FT.nrow();
  int nrow_Comp = Comp_data.nrow();
  int length_PIon = PIon.length();

  Rcpp::CharacterVector feature_name = FT["feature_name"];
  Rcpp::NumericVector FT_mz = FT["mzmed"];
  Rcpp::NumericVector FT_rt = FT["rtmed"];

  Rcpp::CharacterVector Comp_name = Comp_data["Name"];
  Rcpp::NumericVector Comp_rt = Comp_data["RT"];
  Rcpp::NumericVector Comp_M = Comp_data["exactMass"];

  List output(nrow_FT);
  for(int i = 0; i < nrow_FT ; i++)
  {
    Rcpp::List df(6);
    df.names() = Rcpp::CharacterVector::create("Feature_name", "Feature_mz","Feature_rt","Comp_name","PI_name","M");

    CharacterVector feaNames;
    CharacterVector CompNames;
    CharacterVector PIname;
    Rcpp::NumericVector  fea_rt;
    Rcpp::NumericVector  fea_mz;
    Rcpp::NumericVector  M_vec;

    for(int p=0; p<length_PIon ; p++){
      String PIon_name = PIon[p];
      Rcpp::NumericVector comp_mz = Comp_data[PIon_name];

      for(int j = 0; j < nrow_Comp; j++){

        float diff_mz= std::fabs(FT_mz[i]-comp_mz[j]);

        if (((Comp_rt[j] != 0 || !NumericVector::is_na(Comp_rt[j])) &&
            (diff_mz <= comp_mz[j] * ppm * 1e-6) &&
            (std::abs(FT_rt[i] - Comp_rt[j]) <= diff_rt_thr)) ||
            ((diff_mz <= diff_mz_thr) &&
            (std::abs(FT_rt[i] - Comp_rt[j]) <= diff_rt_thr)))
        {
          feaNames.push_back(feature_name[i]);
          fea_rt.push_back(FT_rt[i]);
          fea_mz.push_back(FT_mz[i]);
          CompNames.push_back(Comp_name[j]);
          PIname.push_back(PIon[p]);
          M_vec.push_back(Comp_M[j]);


        }

        else{
          if(diff_mz<= comp_mz[j]*ppm* 1e-6||diff_mz<= diff_mz_thr )
          {
            feaNames.push_back(feature_name[i]);
            fea_rt.push_back(FT_rt[i]);
            fea_mz.push_back(FT_mz[i]);
            CompNames.push_back(Comp_name[j]);
            PIname.push_back(PIon[p]);
            M_vec.push_back(Comp_M[j]);

          }
        }

      }

    }

    df["Feature_name"] = feaNames;
    df["Feature_mz"] = fea_mz;
    df["Feature_rt"] = fea_rt;
    df["Comp_name"] = CompNames;
    df["PI_name"] = PIname;
    df["M"] = M_vec;
    output[i] = df;
  }
  return output;

}



// NL match
// [[Rcpp::export]]
Rcpp::List findNL_fun(
    Rcpp::DataFrame NL_data,
    Rcpp::DataFrame adj_long,
    double diff_mz_thr,
    double diff_rt_thr,
    int ppm) {

  int nrow_adj_long = adj_long.nrow();
  int nrow_NL = NL_data.nrow();

  Rcpp::NumericVector mz_diff = adj_long["mz_diff"];
  Rcpp::NumericVector rt_diff = adj_long["rt_diff"];
  Rcpp::NumericVector NL_mz = NL_data["Accurate.Mass"];
  Rcpp::CharacterVector NL_name = NL_data["Neutral.Loss"];
  Rcpp::CharacterVector fea_name1 = adj_long["Var1"];
  Rcpp::CharacterVector fea_name2 = adj_long["Var2"];

  Rcpp::NumericVector mz1 = adj_long["mz_x"];
  Rcpp::NumericVector mz2 = adj_long["mz_y"];

  List out(nrow_adj_long);

  for (int i = 0; i < nrow_adj_long - 1; i++) {
    Rcpp::List df(3);
    df.names() = Rcpp::CharacterVector::create("Feature_name1", "Feature_name2", "NL_Names");
    CharacterVector feaNames1;
    CharacterVector feaNames2;
    CharacterVector NLNames;

    for (int j = 0; j < nrow_NL - 1; j++) {
      float diff_mz1 = std::fabs(NL_mz[j] - mz_diff[i]);

      if ((diff_mz1 <= diff_mz_thr || diff_mz1 <= NL_mz[j] * ppm * 1e-6) && rt_diff[i] <= diff_rt_thr) {


        if (mz1[i] < mz2[i]) {
          feaNames1.push_back(fea_name2[i]); // swap names
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


